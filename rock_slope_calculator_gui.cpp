// ======================================================================================
// PROGRAM DESCRIPTION & METHODOLOGY
// ======================================================================================
//
// 1. PURPOSE:
//    This software calculates the required size and weight of rock armor units for 
//    rubble mound breakwaters and revetments. It performs a comprehensive analysis 
//    of hydraulic stability across different water depth zones (Deep, Shallow, 
//    Very Shallow, and Swash zones), selecting the most appropriate empirical formula 
//    based on the hydraulic regime.
//
// 2. METHODOLOGY & LOGIC:
//    The calculator implements a "Multi-Model Consensus" approach. It computes stability 
//    using several state-of-the-art empirical formulas derived from extensive hydraulic 
//    model testing. It then evaluates the hydraulic context (Relative Depth h/Hm0) 
//    to recommend the scientifically most accurate formula.
//
//    The logic follows a 4-Step Process:
//    a. Input Processing: Parsing wave data, geometry, and material properties.
//    b. Hydraulic Analysis: Computing local wavelength, celerity, wave steepness, 
//       surf similarity (Iribarren number), and determining the breaker type.
//    c. Stability Calculation: Running multiple empirical models (Hudson, Van der Meer, 
//       Van Gent, Etemad-Shahidi, Eldrup & Andersen, Scaravaglione).
//    d. Intelligent Selection: Automatically detecting the "Hydraulic Zone" 
//       (Deep vs. Shallow vs. Swash) to recommend the most valid result.
//
// 3. HYDRAULIC ZONES & FORMULA SELECTION STRATEGY:
//    - ZONE 1: Deep/Intermediate (h/Hm0 > 3.0) -> Recommends Van der Meer (2021).
//    - ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0) -> Recommends Van der Meer (2021) 
//      or Van Gent (2003), depending on spectral shape.
//    - ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5) -> Recommends Scaravaglione (2025) 
//      (Modified ES) to account for wave breaking and horizontal stability trends.
//    - ZONE 4: Extremely Shallow/Swash (h/Hm0 <= 0.5) -> Recommends Scaravaglione (2025) 
//      (Modified VG) to handle infragravity dominance and high turbulence.
//
// 4. BIBLIOGRAPHY & SCIENTIFIC REFERENCES:
//    (Identical to CLI version - omitted for brevity in GUI source but logic is same)
//
// 5. COMPILATION INSTRUCTIONS (MinGW on Windows):
//
//    g++ -O3 -std=c++17 -municode rock_slope_calculator_gui.cpp \
//    -o rock_slope_calculator_gui.exe -mwindows -static -static-libgcc -static-libstdc++
//
// 6. EXECUTION:
//    1. Launch the application (rock_slope_calculator_gui.exe).
//    2. Enter hydraulic parameters.
//    3. Select "Automatic" or a specific formula.
//    4. Click "Calculate".
//
// ======================================================================================

// ----------------------------------------------------------------------
// MINGW STATIC LINKING SHIM (Fixes "undefined reference" on Windows/GCC)
// ----------------------------------------------------------------------
#if defined(_WIN32) && defined(__GNUC__)
    #include <cstdio>
    extern "C" {
        int (*__imp_fseeko64)(FILE*, long long, int) = reinterpret_cast<int(*)(FILE*, long long, int)>(&_fseeki64);
        long long (*__imp_ftello64)(FILE*) = reinterpret_cast<long long(*)(FILE*)>(&_ftelli64);
    }
#endif

// GUI Headers
#define _USE_MATH_DEFINES 
#include <windows.h>
#include <codecvt> 

// Standard Headers
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>

// Define Pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ----------------------------------------------------------------------
// DATA STRUCTURES
// ----------------------------------------------------------------------

struct GradingDef {
    std::string name;
    double min_mass; // kg
    double max_mass; // kg
    double M50;      // kg
};

struct Inputs {
    double Hs;          // Significant Wave Height (m)
    double Tm10;        // Spectral Period Tm-1,0 (s)
    double h_toe;       // Water Depth at Toe (m)
    double slope_m;     // Structure Slope cot(alpha)
    double foreshore_m; // Foreshore Slope cot(beta)
    double rho_r;       // Rock Density (kg/m3)
    double rho_w;       // Water Density (kg/m3) - Default 1025
    double P;           // Notional Permeability
    double D_ratio;     // Core/Armor Diameter Ratio
    double Sd;          // Design Damage Level
    double duration;    // Storm Duration (hours)
    bool use_en13383;   // Use Standard Grading
    
    // GUI Override Logic (0 = Auto, >0 = Specific Index)
    int formula_override_index; 
};

struct DerivedParams {
    double cot_alpha;
    double alpha_rad;
    double alpha_deg;
    double Delta;       // Relative buoyant density
    double Cp;          // Physical permeability coefficient
    double N_waves;     // Number of waves
};

struct Hydraulics {
    double L0;
    double L_toe;
    double C;
    double Cg;
    double s_m10;       // Deep water steepness
    double s_local;     // Local steepness
    double xi_m10;      // Surf similarity parameter
    double rel_depth;   // h / Hm0
    std::string breaker_type;
    std::string zone_desc;
};

struct FormulaResult {
    std::string name;
    double Dn50;
    double Ns;
    double Kd;
    std::string note;
    bool valid;         // If formula calculation was successful (Dn > 0)
};

struct LayerDesign {
    std::string layer_name; // "Primary Armor" or "Underlayer"
    std::string grading_name;
    double target_W_kN;
    double target_M50_kg;
    double target_Dn_m;
    
    double w_min_kn;
    double w_max_kn;
    double w_min_kg;
    double w_max_kg;
    
    double m_mean_kg;
    double w_mean_kn;
    double actual_dn;
    double thickness;
    double packing_density; // rocks/100m2
    
    bool design_valid; // If a valid grading/design was found
};

struct FullReport {
    Inputs inputs;
    DerivedParams derived;
    Hydraulics hydro;
    std::vector<FormulaResult> comparison;
    FormulaResult recommended;
    std::vector<std::string> justification;
    LayerDesign armor_layer;
    LayerDesign underlayer;
};

// ----------------------------------------------------------------------
// CLASS DEFINITION
// ----------------------------------------------------------------------

class RockSlopeCalculator {
private:
    const double G = 9.80665;
    std::vector<GradingDef> standard_gradings;
    std::vector<std::string> log_buffer;
    
    // GUI stream capture
    std::stringstream gui_ss;

public:
    Inputs defaults;

    RockSlopeCalculator() {
        // Default Inputs
        defaults = {
            2.0,    // Hs
            12.0,   // Tm10
            7.0,    // h_toe
            2.0,    // slope_m
            30.0,   // foreshore_m
            2650.0, // rho_r
            1025.0, // rho_w
            0.4,    // P
            0.3,    // D_ratio
            2.0,    // Sd
            6.0,    // duration
            true,   // use_en13383
            0       // formula_override_index (0 = Auto)
        };

        // Initialize EN 13383 Database (Mass in kg)
        standard_gradings = {
            {"CP 45/125",       0.4,   1.2,    0.8},
            {"CP 63/180",       1.2,   3.8,    2.5},
            {"CP 90/250",       3.1,   9.3,    6.2},
            {"CP 45/180",       0.4,   1.2,    0.8},
            {"CP 90/180",       2.1,   2.8,    2.45},
            {"LMA 5-40",        10,    20,     15},
            {"LMA 10-60",       20,    35,     27.5},
            {"LMA 15-120",      35,    60,     47.5},
            {"LMA 40-200",      80,    120,    100},
            {"LMA 60-300",      120,   190,    155},
            {"LMA 15-300",      45,    135,    90},
            {"HMA 300-1000",    540,   690,    615},
            {"HMA 1000-3000",   1700,  2100,   1900},
            {"HMA 3000-6000",   4200,  4800,   4500},
            {"HMA 6000-10000",  7500,  8500,   8000},
            {"HMA 10000-15000", 12000, 13000,  12500}
        };

        std::sort(standard_gradings.begin(), standard_gradings.end(), 
            [](const GradingDef& a, const GradingDef& b) {
                return a.M50 < b.M50;
            });
    }

    std::string get_gui_output() const {
        return gui_ss.str();
    }

    void log(std::string message) {
        // Echo to stdout for debug/CLI consistency
        std::cout << message << std::endl;
        // Save to buffer
        log_buffer.push_back(message);
        // Save to GUI buffer
        gui_ss << message << "\n";
    }

    // --- Wave Mechanics ---

    double solve_wavelength(double T, double h) {
        if (h <= 0) return 0.0;
        double L0 = G * T * T / (2 * M_PI);
        double k0h = 2 * M_PI * h / L0;
        
        // Initial Guess
        double term = std::pow(6.0/5.0, k0h) * std::sqrt(k0h);
        double L = L0 * std::tanh(term);
        
        // Newton-Raphson
        double dL = 1.0;
        double delta = 1e-8;
        int iter = 0;
        while (std::abs(dL / L) >= delta && iter < 100) {
            double val = 2 * M_PI * h / L;
            double f1 = L - L0 * std::tanh(val);
            double val_delta = 2 * M_PI * h / (L + delta);
            double f2 = (L + delta) - L0 * std::tanh(val_delta);
            double denom = f2 - f1;
            if (denom == 0) break;
            dL = delta * f1 / denom;
            L = L - dL;
            iter++;
        }
        return L;
    }

    Hydraulics analyze_hydraulics(const Inputs& in, const DerivedParams& dp) {
        Hydraulics h;
        h.L0 = (G * std::pow(in.Tm10, 2)) / (2 * M_PI);
        h.L_toe = solve_wavelength(in.Tm10, in.h_toe);
        
        if (in.Tm10 > 0) h.C = h.L_toe / in.Tm10;
        else h.C = 0;
        
        double k = (h.L_toe > 0) ? (2 * M_PI / h.L_toe) : 0;
        double kh = k * in.h_toe;
        double n;
        if (kh > 20) n = 0.5;
        else if (kh <= 0) n = 1.0;
        else n = 0.5 * (1 + (2 * kh) / std::sinh(2 * kh));
        h.Cg = n * h.C;

        h.s_m10 = in.Hs / h.L0;
        h.s_local = (h.L_toe > 0) ? (in.Hs / h.L_toe) : 0;
        
        if (h.s_m10 > 0) 
            h.xi_m10 = std::tan(dp.alpha_rad) / std::sqrt(h.s_m10);
        else 
            h.xi_m10 = 0;

        h.rel_depth = (in.Hs > 0) ? (in.h_toe / in.Hs) : 999.0;

        if (h.xi_m10 < 0.5) h.breaker_type = "Spilling";
        else if (h.xi_m10 < 1.8) h.breaker_type = "Plunging";
        else if (h.xi_m10 < 3.0) h.breaker_type = "Surging";
        else h.breaker_type = "Collapsing/Surging";

        if (h.rel_depth > 3.0) h.zone_desc = "ZONE 1: Deep to Intermediate (h/Hm0 > 3.0)";
        else if (h.rel_depth > 1.5) h.zone_desc = "ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0)";
        else if (h.rel_depth > 0.5) h.zone_desc = "ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5)";
        else h.zone_desc = "ZONE 4: Extremely Shallow Water (h/Hm0 <= 0.5)";

        return h;
    }

    // --- Stability Formulas ---

    FormulaResult calc_hudson(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        double Kd;
        if (h.rel_depth > 3.0) Kd = 4.0;
        else if (h.rel_depth > 1.5) Kd = 3.5;
        else if (h.rel_depth > 0.5) Kd = 3.0;
        else Kd = 2.0;

        double Dn = (1.27 * p.Hs) / (dp.Delta * std::pow(Kd * dp.cot_alpha, 1.0/3.0));
        double Ns = (1.27 * p.Hs) / (dp.Delta * Dn);
        return {"Hudson (1959)", Dn, Ns, Kd, "Hudson (1959) - Legacy Ref", true};
    }

    FormulaResult calc_vdm_2021(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        double c_pl = 6.49, c_su = 0.97;
        double term_crit = (c_pl / c_su) * std::pow(p.P, 0.31) * std::sqrt(std::tan(dp.alpha_rad));
        double xi_cr = std::pow(term_crit, 1.0 / (p.P + 0.5));
        double damage_term = std::pow(p.Sd / std::sqrt(dp.N_waves), 0.2);
        
        double Ns;
        std::string note;
        std::stringstream ss; ss << std::fixed << std::setprecision(2) << xi_cr;

        if (h.xi_m10 < xi_cr) {
            Ns = c_pl * std::pow(p.P, 0.18) * damage_term * std::pow(h.xi_m10, -0.5);
            note = "Van der Meer 2021 (Plunging: xi < " + ss.str() + ")";
        } else {
            Ns = c_su * std::pow(p.P, -0.13) * damage_term * std::sqrt(dp.cot_alpha) * std::pow(h.xi_m10, p.P);
            note = "Van der Meer 2021 (Surging: xi > " + ss.str() + ")";
        }
        
        double Dn = p.Hs / (dp.Delta * Ns);
        double Kd = std::pow(Ns, 3) / dp.cot_alpha;
        return {"Van der Meer (2021)", Dn, Ns, Kd, note, true};
    }

    FormulaResult calc_van_gent_mod(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        double h2_ratio;
        if (h.rel_depth >= 3.0) {
            h2_ratio = 1.4;
        } else if (h.rel_depth < 1.5) {
            h2_ratio = 1.2 + 0.2 * (1.5 - h.rel_depth) / 1.5;
            h2_ratio = std::min(h2_ratio, 1.4);
        } else {
            h2_ratio = 1.2 + 0.2 * (h.rel_depth - 1.5) / 1.5;
        }

        double c_pl = 8.4;
        double c_su = 1.3;
        
        double term_crit = (c_pl / c_su) * std::pow(p.P, 0.31) * std::sqrt(std::tan(dp.alpha_rad));
        double xi_cr = std::pow(term_crit, 1.0 / (p.P + 0.5));
        double damage_term = std::pow(p.Sd / std::sqrt(dp.N_waves), 0.2);
        double ratio_term = std::pow(h2_ratio, -1.0);

        double Ns;
        std::string note;
        std::stringstream ss; ss << std::fixed << std::setprecision(2) << h2_ratio;

        if (h.xi_m10 < xi_cr) {
            Ns = 8.4 * std::pow(p.P, 0.18) * damage_term * std::pow(h.xi_m10, -0.5) * ratio_term;
            note = "Van Gent Modified (2003) Plunging (H2%/Hs=" + ss.str() + ")";
        } else {
            Ns = 1.3 * std::pow(p.P, -0.13) * damage_term * std::sqrt(dp.cot_alpha) * std::pow(h.xi_m10, p.P) * ratio_term;
            note = "Van Gent Modified (2003) Surging (H2%/Hs=" + ss.str() + ")";
        }

        double Dn = p.Hs / (dp.Delta * Ns);
        double Kd = std::pow(Ns, 3) / dp.cot_alpha;
        return {"Van Gent Modified (2003)", Dn, Ns, Kd, note, true};
    }

    FormulaResult calc_van_gent_simp(const Inputs& p, const DerivedParams& dp, const Hydraulics&) {
        double c_VG = 1.75;
        double damage_term = std::pow(p.Sd / std::sqrt(dp.N_waves), 0.2);
        double perm_term = 1.0 + p.D_ratio;
        
        double Ns = c_VG * std::sqrt(dp.cot_alpha) * perm_term * damage_term;
        double Dn = p.Hs / (dp.Delta * Ns);
        double Kd = std::pow(Ns, 3) / dp.cot_alpha;
        return {"Van Gent Simplified (2003)", Dn, Ns, Kd, "Van Gent Simplified (2003)", true};
    }

    FormulaResult calc_eldrup(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        double c_EA1 = 4.5, c_EA2 = 3.1;
        double damage_term = std::pow(p.Sd / std::sqrt(dp.N_waves), 0.2);
        
        double Ns_pl = c_EA1 * damage_term * std::pow(1.6, p.P) * std::pow(h.xi_m10, 0.4*p.P - 0.67);
        double cot_term = std::pow(std::min(dp.cot_alpha, 2.0), 0.23);
        double Ns_su = c_EA2 * damage_term * std::pow(p.P, 0.17) * cot_term;
        
        double Ns;
        std::string note;
        if (h.xi_m10 < 2.5) { Ns = Ns_pl; note = "Eldrup (Plunging)"; }
        else { Ns = Ns_su; note = "Eldrup (Surging)"; }
        
        double Dn = p.Hs / (dp.Delta * Ns);
        double Kd = std::pow(Ns, 3) / dp.cot_alpha;
        return {"Eldrup & Andersen (2019)", Dn, Ns, Kd, note, true};
    }

    FormulaResult calc_es_2020(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        double c_ES1 = 4.5, c_ES2 = 3.9;
        double m = 1.0 / p.foreshore_m;
        double term_Nw = std::pow(dp.N_waves, -0.1);
        double term_Sd = std::pow(p.Sd, 1.0/6.0);
        
        double fs_factor; 
        std::string suffix;

        if (h.rel_depth < 3.0) {
            fs_factor = std::max(0.1, 1.0 - 3.0*m);
            suffix = " (Depth-Limited)";
        } else {
            fs_factor = 1.0;
            suffix = " (Deep)";
        }
        
        double Ns;
        std::string note;
        if (h.xi_m10 < 1.8) {
            Ns = c_ES1 * dp.Cp * term_Nw * term_Sd * std::pow(h.xi_m10, -7.0/12.0) * fs_factor;
            note = "Etemad-Shahidi (Plunging)" + suffix;
        } else {
            Ns = c_ES2 * dp.Cp * term_Nw * term_Sd * std::pow(h.xi_m10, -1.0/3.0) * fs_factor;
            note = "Etemad-Shahidi (Surging)" + suffix;
        }
        
        double Dn = p.Hs / (dp.Delta * Ns);
        double Kd = std::pow(Ns, 3) / dp.cot_alpha;
        return {"Etemad-Shahidi (2020)", Dn, Ns, Kd, note, true};
    }

    FormulaResult calc_mod_vg(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        double c_VG_new = 3.3;
        double perm_term = 1.0 + p.D_ratio;
        double s_safe = std::max(h.s_m10, 0.005);
        double steep_term = std::pow(s_safe, 0.1);
        double damage_term = std::pow(p.Sd / std::sqrt(dp.N_waves), 0.2);
        
        double Ns = c_VG_new * std::sqrt(dp.cot_alpha) * perm_term * steep_term * damage_term;
        
        // Safety Cap
        double Kd_calc = std::pow(Ns, 3) / dp.cot_alpha;
        std::string note = "Scaravaglione (Mod. VG)";
        if (Kd_calc > 5.0) {
            Ns = std::pow(5.0 * dp.cot_alpha, 1.0/3.0);
            Kd_calc = 5.0; 
            note += " [Capped Kd=5.0]";
        }
        double Dn = p.Hs / (dp.Delta * Ns);
        double Kd = std::pow(Ns, 3) / dp.cot_alpha;
        return {"Scaravaglione (Mod. VG 2025)", Dn, Ns, Kd, note, true};
    }

    FormulaResult calc_mod_es(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        if (h.xi_m10 < 1.8) return {"Scaravaglione (Mod. ES 2025)", 0, 0, 0, "N/A (Req. xi > 1.8)", false};
        
        double c_ES_new = 3.55;
        double term_N = std::pow(dp.N_waves, -0.1);
        double term_Sd = std::pow(p.Sd, 1.0/6.0);
        double term_cot = std::pow(dp.cot_alpha, 1.0/3.0);
        double s_safe = std::max(h.s_m10, 0.005);
        double term_s = std::pow(s_safe, 1.0/20.0);
        
        double Ns = c_ES_new * dp.Cp * term_N * term_cot * term_Sd * term_s;
        
        // Safety Cap
        double Kd_calc = std::pow(Ns, 3) / dp.cot_alpha;
        std::string note = "Scaravaglione (Mod. ES)";
        if (Kd_calc > 5.0) {
            Ns = std::pow(5.0 * dp.cot_alpha, 1.0/3.0);
            note += " [Capped Kd=5.0]";
        }
        double Dn = p.Hs / (dp.Delta * Ns);
        double Kd = std::pow(Ns, 3) / dp.cot_alpha;
        return {"Scaravaglione (Mod. ES 2025)", Dn, Ns, Kd, note, true};
    }

    // --- Recommendation Logic ---

    void select_recommendation(FullReport& report) {
        double depth_ratio = report.hydro.rel_depth;
        std::vector<std::string>& log = report.justification;

        auto get_res = [&](std::string key) -> FormulaResult {
            for (auto& r : report.comparison) {
                if (r.name.find(key) != std::string::npos && r.valid) return r;
            }
            return {"", 0,0,0,"", false};
        };

        FormulaResult vdm = get_res("Van der Meer (2021)");
        FormulaResult vg_mod = get_res("Van Gent Modified");
        FormulaResult vg_simp = get_res("Van Gent Simplified");
        FormulaResult es = get_res("Etemad-Shahidi");
        FormulaResult mod_vg = get_res("Mod. VG");
        FormulaResult mod_es = get_res("Mod. ES");

        // Logic Tree replicating Python's intelligent selector
        if (depth_ratio > 3.0) {
            // ZONE 1: Deep
            report.recommended = vdm;
            log.push_back("### 1. HYDRAULIC CONTEXT: Deep / Intermediate Water");
            std::stringstream ss; ss << std::fixed << std::setprecision(2) << depth_ratio;
            log.push_back("   The structure is located in deep water relative to the wave height (h/Hm0 = " + ss.str() + " > 3.0).");
            log.push_back("   In this regime, the wave height distribution strictly follows the **Rayleigh distribution**.");
            log.push_back("   Key characteristics:");
            log.push_back("     * The ratio H2%/Hs is constant at approximately 1.4.");
            log.push_back("     * Wave breaking is limited to whitecapping or direct interaction with the armor layer.");
            log.push_back("     * The spectral shape is standard (JONSWAP/Pierson-Moskowitz), and energy transfer to low-frequencies is minimal.");
            log.push_back("");
            log.push_back("### 2. FORMULA COMPARISON & ANALYSIS");
            
            log.push_back("   **A. Van der Meer (2021 Rewritten) [RECOMMENDED]**");
            log.push_back("      * **Advantages:** This formula is the modernized industry standard.");
            log.push_back("        Van der Meer (2021) rewrote the original formula to use the spectral period (Tm-1,0),");
            log.push_back("        eliminating the influence of spectral shape.");
            log.push_back("        Van der Meer et al. (2024) confirmed its validity for h/Hm0 > 1.5, preferring Hm0 over H1/3 for nonlinear waves.");
            log.push_back("      * **Physics:** It correctly assumes a Rayleigh distribution of wave heights, aligning with the actual deep-water statistics.");

            log.push_back("   **B. Van Gent Modified (2003)**");
            log.push_back("      * **Context:** This formula incorporates the ratio H2%/Hs.");
            log.push_back("        In deep water, with H2%/Hs = 1.4, this formula essentially converges closely with the Van der Meer predictions.");
            log.push_back("        However, its specific calibration was focused on the effects of shallow foreshores.");

            log.push_back("   **C. Etemad-Shahidi (2020)**");
            log.push_back("      * **Comparison:** Etemad-Shahidi (2020) provides a robust formula validated for both deep and shallow water.");
            log.push_back("        It introduces a physical permeability parameter (D_core/D_armor) to replace the nominal P factor,");
            log.push_back("        reducing uncertainty. However, Van der Meer remains the primary standard for deep water.");

            if (es.valid && vdm.valid && vdm.Dn50 > 0) {
                double diff_pct = std::abs(vdm.Dn50 - es.Dn50) / vdm.Dn50;
                if (diff_pct > 0.10) {
                    log.push_back("      * **Note on Divergence:** The result deviates from Van der Meer here. This typically occurs in the 'transition zone'");
                    log.push_back("        of the surf similarity parameter (xi approx 2.0 - 4.5). Etemad-Shahidi transitions to Surging physics");
                    log.push_back("        earlier (xi > 1.8), predicting higher stability, whereas Van der Meer maintains Plunging physics");
                    log.push_back("        (lower stability) until a higher critical threshold. Van der Meer is more conservative here.");
                } else {
                    log.push_back("        It typically converges with Van der Meer here.");
                }
            }

            log.push_back("");
            log.push_back("### 3. FINAL JUSTIFICATION");
            log.push_back("   **Use [Van der Meer (2021 Rewritten)]**.");
            log.push_back("   It provides the most theoretically consistent result for non-depth-limited waves.");
            
            if (es.valid && vdm.valid) {
                 std::stringstream ss_diff; ss_diff << std::fixed << std::setprecision(3) << es.Dn50;
                 std::stringstream ss_delta; ss_delta << std::fixed << std::setprecision(3) << std::abs(vdm.Dn50 - es.Dn50);
                 log.push_back("   *Verification:* Etemad-Shahidi yields Dn50 = " + ss_diff.str() + "m (Difference: " + ss_delta.str() + "m).");
            }
        } 
        else if (depth_ratio > 1.5) {
            // ZONE 2: Shallow
            report.recommended = vdm;
            std::stringstream ss; ss << std::fixed << std::setprecision(2) << depth_ratio;
            log.push_back("### 1. HYDRAULIC CONTEXT: Shallow Water (Transition Zone)");
            log.push_back("   The structure is in the transition zone (1.5 < h/Hm0 = " + ss.str() + " <= 3.0).");
            log.push_back("   Key characteristics:");
            log.push_back("     * **Spectral Truncation:** The largest waves in the spectrum break on the foreshore.");
            log.push_back("     * **Distribution Shift:** The wave height distribution deviates from Rayleigh; H2%/Hm0 drops below 1.4.");
            log.push_back("     * **Shoaling:** Significant shoaling modifies the wave shape before impact, creating peaked crests and flat troughs.");
            log.push_back("");
            log.push_back("### 2. FORMULA COMPARISON & ANALYSIS");
            
            log.push_back("   **A. Van der Meer (2021)**");
            log.push_back("      * **Advantages:** Van der Meer et al. (2024) extensively re-analyzed shallow water data and concluded");
            log.push_back("        that the rewritten Van der Meer formula (using Tm-1,0) is valid down to h/Hm0 = 1.5.");
            log.push_back("        It performs reasonably well, with slightly less reliability in the 1.0 < h/Hm0 < 1.5 range.");
            log.push_back("      * **Note:** For nonlinear waves in this zone, using Hm0 is preferred over H1/3 for nonlinear waves to avoid deviations.");

            log.push_back("   **B. Van Gent Modified (2003)**");
            log.push_back("      * **Constraint:** This formula explicitly relies on the ratio H2%/Hs. Research by Van der Meer et al. (2024)");
            log.push_back("        highlights that predicting H2% accurately in this transition zone (where the ratio dips to ~1.2)");
            log.push_back("        is notoriously inaccurate without physical modeling. The formula is valid, but the input uncertainty is high.");

            log.push_back("   **C. Van Gent et al. (2003) Simplified**");
            log.push_back("      * **Context:** This formula was specifically derived for shallow foreshores.");
            log.push_back("        However, Van der Meer et al. (2024) found that the simplified formula often does not match the data");
            log.push_back("        in the surging domain as well as the rewritten Van der Meer formula.");
            
            log.push_back("");
            log.push_back("### 3. FINAL JUSTIFICATION");
            log.push_back("   **Use [Van der Meer (2021 Rewritten)]**.");
            log.push_back("   Recent research (2024) confirms its validity in this depth range (h/Hm0 > 1.5), favoring it over simplified methods");
            log.push_back("   due to the uncertainties in predicting H2% required for the Van Gent Modified formula.");
        }
        else if (depth_ratio > 0.5) {
            // ZONE 3: Very Shallow
            std::stringstream ss; ss << std::fixed << std::setprecision(2) << depth_ratio;
            
            if (mod_es.valid) {
                report.recommended = mod_es;
                log.push_back("### 1. HYDRAULIC CONTEXT: Very Shallow Water (Surf Zone)");
                log.push_back("   The structure is in the surf zone (0.5 < h/Hm0 = " + ss.str() + " <= 1.5).");
                log.push_back("   Key characteristics:");
                log.push_back("     * **Severe Breaking:** Waves are constantly breaking.");
                log.push_back("     * **Saturation:** Wave height is depth-limited (H ~ 0.5h). Increasing offshore energy does not increase load.");
                log.push_back("     * **Infragravity Dominance:** Scaravaglione et al. (2025) and VdM (2024) note that infragravity waves");
                log.push_back("       begin to dominate the spectrum, causing Tm-1,0 to increase massively (up to 4x).");
                log.push_back("     * **Formula Deviation:** Standard formulas fail here because the stability curves flatten out (Horizontal Trend).");
                log.push_back("");
                log.push_back("### 2. FORMULA COMPARISON & ANALYSIS");
                
                log.push_back("   **A. Standard Formulas (VdM, VG, Standard ES)**");
                log.push_back("      * **Failure Mode:** Scaravaglione et al. (2025) demonstrated that these formulas fail to converge here.");
                log.push_back("        Using the inflated spectral period results in over-predicted stability (for surging) or under-predicted (for plunging).");
                
                log.push_back("   **B. Scaravaglione (Modified ES 2025) [RECOMMENDED]**");
                log.push_back("      * **Advantages:** This formula explicitly **decouples** the wave steepness term (s^0.05) from the structure slope term.");
                log.push_back("      * **Physics:** It is calibrated specifically for surging/bore-like waves in the surf zone using new coefficients (c_ES,new=3.55).");
                log.push_back("      * **Safety Cap Applied:** The formula has a weak dependence on steepness. For very long swells, it may predict");
                log.push_back("        theoretically high stability (Kd > 10). This system has capped Kd at 5.0 to ensure physical realism.");
                
                log.push_back("");
                log.push_back("### 3. FINAL JUSTIFICATION");
                log.push_back("   **Use [Scaravaglione (Modified ES 2025)]**.");
                log.push_back("   This represents the state-of-the-art for broken waves in very shallow water, correcting the overestimation of damage.");
            } else {
                report.recommended = vg_simp;
                log.push_back("### 1. HYDRAULIC CONTEXT: Very Shallow Water (Surf Zone) (Plunging)");
                log.push_back("   The structure is in the surf zone, but conditions are **Plunging** (xi < 1.8).");
                log.push_back("   The Modified ES formula is only calibrated for surging/bore conditions.");
                log.push_back("");
                log.push_back("### 3. FINAL JUSTIFICATION");
                log.push_back("   **Use [Van Gent Simplified (2003)]**.");
                log.push_back("   It acts as a robust fallback. Caution is advised as damage may be underpredicted for impermeable structures.");
            }
        }
        else {
            // ZONE 4: Swash
            report.recommended = mod_vg;
            std::stringstream ss; ss << std::fixed << std::setprecision(2) << depth_ratio;
            log.push_back("### 1. HYDRAULIC CONTEXT: Extremely Shallow Water (Swash Zone)");
            log.push_back("   The structure is located effectively in the **Swash Zone** (h/Hm0 = " + ss.str() + " <= 0.5).");
            log.push_back("   Key characteristics:");
            log.push_back("     * **Aeration:** High air entrainment reduces the effective fluid density and buoyancy of the rocks.");
            log.push_back("     * **Impact:** Wave impact is characterized by a high-velocity turbulent bore.");
            log.push_back("     * **Hydrostatics:** The hydrostatic cushioning effect is negligible.");
            log.push_back("");
            log.push_back("### 2. FORMULA COMPARISON & ANALYSIS");
            
            log.push_back("   **A. Standard Van Gent (2003)**");
            log.push_back("      * **Disadvantages:** Using the standard coefficient (1.75) here is **Unsafe**.");
            log.push_back("        Scaravaglione et al. (2025) showed that stability is significantly lower than predicted by intermediate-depth formulas");
            log.push_back("        due to the lack of buoyancy and intense turbulence.");
            
            log.push_back("   **B. Scaravaglione (Modified VG 2025) [RECOMMENDED]**");
            log.push_back("      * **Advantages:** This formula uses a recalibrated coefficient (C_VG = 3.3 instead of 1.75).");
            log.push_back("      * **Physics:** It explicitly accounts for the increased instability in the swash zone, correcting the underestimation");
            log.push_back("        of damage by the original VG formula in this specific regime.");
            
            log.push_back("");
            log.push_back("### 3. FINAL JUSTIFICATION");
            log.push_back("   **Use [Scaravaglione (Modified VG 2025)]**.");
            log.push_back("   It provides the necessary safety margin for swash zone instability where standard formulas fail.");
        }
    }

    // --- Layer Design ---

    LayerDesign design_layer(double target_mass, double target_dn, bool is_armor) {
        LayerDesign ld;
        ld.layer_name = is_armor ? "Primary Armor" : "Underlayer";
        ld.target_M50_kg = target_mass;
        ld.target_Dn_m = target_dn;
        ld.target_W_kN = target_mass * G / 1000.0;
        
        double gamma_r = defaults.rho_r * G / 1000.0;
        bool grading_EN13383 = defaults.use_en13383;
        
        // --- EN 13383 Standard Grading Logic ---
        if (grading_EN13383) {
            // Find standard grading
            GradingDef selected = {"", 0,0,0};
            bool found = false;
            
            if (is_armor) {
                // For Armor: Select lightest where M50 >= Target
                for (const auto& g : standard_gradings) {
                    if (g.M50 >= target_mass) { selected = g; found = true; break; }
                }
            } else {
                // For Underlayer: Select closest M50
                double min_diff = std::numeric_limits<double>::infinity();
                for (const auto& g : standard_gradings) {
                    double diff = std::abs(g.M50 - target_mass);
                    if (diff < min_diff) { min_diff = diff; selected = g; found = true; }
                }
            }

            if (found) {
                ld.grading_name = selected.name;
                ld.m_mean_kg = selected.M50;
                ld.w_min_kn = selected.min_mass * G / 1000.0;
                ld.w_max_kn = selected.max_mass * G / 1000.0;
                ld.w_min_kg = selected.min_mass;
                ld.w_max_kg = selected.max_mass;
                ld.w_mean_kn = ld.m_mean_kg * G / 1000.0;
                ld.actual_dn = std::pow(ld.w_mean_kn / gamma_r, 1.0/3.0);
                ld.design_valid = true;
            } else {
                ld.grading_name = "No Standard Fit";
                ld.design_valid = false;
            }
        } 
        
        // --- Custom Power Law Calculation ---
        if (!grading_EN13383 || !ld.design_valid) {
             ld.grading_name = is_armor ? "Custom Grading" : "Custom Grading Underlayer";
             
             // Using 'x' as Theoretical Required M50 (in kg)
             double x_val = ld.target_M50_kg;
             
             // --- MINIMUM LIMIT ---
             // Grading Min Params
             double a_min = 1.056832014477894E+00;
             double b_min = 1.482769823574055E+00;
             double c_min = -2.476127406338004E-01;
             
             // Calculate Scaling Factor (Dimensionless Ratio)
             double factor_min = a_min / (1.0 + std::pow(x_val / b_min, c_min));
             
             // Calculate Mass first (kg), then Weight (kN)
             // This ensures the 'requires weights in kgs' rule is applied correctly
             double w_min_kg_calc = target_mass * factor_min;
             ld.w_min_kn = (w_min_kg_calc * G) / 1000.0;
             
             // --- MAXIMUM LIMIT ---
             // Grading Max Params
             double a_max = 1.713085676568561E+00;
             double b_max = 2.460481255856126E+05;
             double c_max = 1.327263214034671E-01;
             
             // Calculate Scaling Factor (Dimensionless Ratio)
             double factor_max = a_max / (1.0 + std::pow(x_val / b_max, c_max));
             
             // Calculate Mass first (kg), then Weight (kN)
             double w_max_kg_calc = target_mass * factor_max;
             ld.w_max_kn = (w_max_kg_calc * G) / 1000.0;
             
             // --- DESIGN VALUES ---
             ld.w_mean_kn = ld.target_W_kN;
             ld.m_mean_kg = target_mass;
             ld.actual_dn = target_dn;
             ld.w_min_kg = w_min_kg_calc;
             ld.w_max_kg = w_max_kg_calc;
             ld.design_valid = true;
        }

        // Geometry
        ld.thickness = 2.0 * 1.0 * ld.actual_dn; // 2 layers, kt=1.0
        double porosity = 0.30;
        ld.packing_density = 100.0 * 2.0 * 1.0 * (1.0 - porosity) / std::pow(ld.actual_dn, 2);
        
        return ld;
    }

    // --- Main Solve Routine ---

    FullReport solve(Inputs in) {
        FullReport report;
        report.inputs = in;
        defaults = in; // Update internal defaults for access in sub-funcs

        // Derived Params
        DerivedParams dp;
        dp.cot_alpha = in.slope_m;
        dp.alpha_rad = std::atan(1.0 / in.slope_m);
        dp.alpha_deg = dp.alpha_rad * 180.0 / M_PI;
        dp.Delta = (in.rho_r / in.rho_w) - 1.0;
        dp.Cp = std::pow(1.0 + std::pow(in.D_ratio, 0.3), 0.6);
        dp.N_waves = (in.duration * 3600.0) / in.Tm10;
        report.derived = dp;

        // Logging Inputs
        log_buffer.clear(); 
        gui_ss.str(""); // Clear GUI buffer
        
        log("ROCK SLOPE STABILITY CALCULATOR");
        log("\n" + std::string(95, '='));
        log("   1. DESIGN INPUT PARAMETERS");
        log(std::string(95, '='));
        
        // Helper to format table rows exactly like output.txt
        auto fmt_row = [](std::string lab, double val, std::string unit, int prec) {
            std::stringstream ss;
            ss << std::left << std::setw(35) << lab << " | " 
               << std::left << std::setw(10) << std::fixed << std::setprecision(prec) << val << " | " << unit;
            return ss.str();
        };

        std::stringstream ss_head_inputs;
        ss_head_inputs << std::left << std::setw(35) << "PARAMETER" << " | " << std::setw(10) << "VALUE" << " | UNIT";
        log(ss_head_inputs.str());
        
        log(std::string(95, '-'));
        log(fmt_row("Significant Wave Height (Hs)", in.Hs, "m", 2));
        log(fmt_row("Spectral Period (Tm-1,0)", in.Tm10, "s", 2));
        log(fmt_row("Water Depth (h_toe)", in.h_toe, "m", 2));
        
        std::stringstream ss_slope; ss_slope << "1:" << std::fixed << std::setprecision(2) << in.slope_m;
        std::stringstream ss_slope_out;
        ss_slope_out << std::left << std::setw(35) << "Structure Slope" << " | " << std::left << std::setw(10) << ss_slope.str() << " | (V:H)";
        log(ss_slope_out.str());

        std::stringstream ss_fore; ss_fore << "1:" << std::fixed << std::setprecision(2) << in.foreshore_m;
        std::stringstream ss_fore_out;
        ss_fore_out << std::left << std::setw(35) << "Foreshore Slope" << " | " << std::left << std::setw(10) << ss_fore.str() << " | (V:H)";
        log(ss_fore_out.str());
        
        log(fmt_row("Permeability (P_notional)", in.P, "(-)", 2));
        log(fmt_row("Physical Permeability (Cp)", dp.Cp, "(-)", 2));
        log(fmt_row("Damage Level (S)", in.Sd, "(-)", 1));
        log(std::string(95, '-'));

        // Hydraulics
        report.hydro = analyze_hydraulics(in, dp);
        
        log("\n" + std::string(95, '='));
        log("   2. CALCULATED HYDRAULIC PARAMETERS");
        log(std::string(95, '='));
        
        std::stringstream ss_head_hydro;
        ss_head_hydro << std::left << std::setw(40) << "PARAMETER" << " | " << std::setw(10) << "VALUE" << " | UNIT";
        log(ss_head_hydro.str());
        
        log(std::string(95, '-'));
        
        auto fmt_h = [](std::string lab, double val, std::string unit, int prec) {
             std::stringstream ss;
             ss << std::left << std::setw(40) << lab << " | " << std::left << std::setw(10) << std::fixed << std::setprecision(prec) << val << " | " << unit;
             return ss.str();
        };

        log(fmt_h("Deep Water Wavelength (L0)", report.hydro.L0, "m", 2));
        log(fmt_h("Wavelength at Toe (L_toe)", report.hydro.L_toe, "m", 2));
        log(fmt_h("Wave Celerity at Toe (C)", report.hydro.C, "m/s", 2));
        log(fmt_h("Group Celerity at Toe (Cg)", report.hydro.Cg, "m/s", 2));
        log(fmt_h("Deep Water Steepness (s_m-1,0)", report.hydro.s_m10, "(-)", 5));
        log(fmt_h("Local Wave Steepness (s_local)", report.hydro.s_local, "(-)", 5));
        log(fmt_h("Surf Similarity (xi_m-1,0)", report.hydro.xi_m10, "(-)", 2));
        
        std::stringstream ss_break; 
        ss_break << std::left << std::setw(40) << "Breaker Type (Visual/Physical)" << " | " << std::left << std::setw(10) << report.hydro.breaker_type << " | (-)";
        log(ss_break.str());
        
        log(fmt_h("Relative Depth (h/Hm0)", report.hydro.rel_depth, "(-)", 2));
        log(fmt_h("Relative Depth (h/L0)", report.hydro.rel_depth * report.hydro.s_m10, "(-)", 3));
        
        std::stringstream ss_zone;
        ss_zone << std::left << std::setw(40) << "Hydraulic Zone" << " | " << std::left << std::setw(40) << report.hydro.zone_desc;
        log(ss_zone.str());
        
        log(std::string(95, '-'));

        // Run Formulas
        report.comparison.push_back(calc_hudson(in, dp, report.hydro));
        report.comparison.push_back(calc_vdm_2021(in, dp, report.hydro));
        report.comparison.push_back(calc_van_gent_mod(in, dp, report.hydro));
        report.comparison.push_back(calc_van_gent_simp(in, dp, report.hydro));
        report.comparison.push_back(calc_eldrup(in, dp, report.hydro));
        report.comparison.push_back(calc_es_2020(in, dp, report.hydro));
        report.comparison.push_back(calc_mod_vg(in, dp, report.hydro));
        report.comparison.push_back(calc_mod_es(in, dp, report.hydro));

        log("\n" + std::string(95, '='));
        log("   3. FORMULA SELECTION & JUSTIFICATION");
        log(std::string(95, '='));
        
        // Recommend
        select_recommendation(report);

        // Comparison Table
        log("COMPARISON OF RESULTS:");
        std::stringstream header;
        header << std::left << std::setw(30) << "Method" 
               << " | " << std::left << std::setw(8) << "Ns (-)" 
               << " | " << std::left << std::setw(10) << "Dn50 (m)" 
               << " | " << std::left << std::setw(10) << "M50 (kg)" 
               << " | " << std::left << std::setw(8) << "Kd_eq (-)" 
               << " | " << "NOTES";
        log(header.str());
        log(std::string(115, '-'));
        
        for (const auto& res : report.comparison) {
            if (!res.valid) continue;
            double mass = in.rho_r * std::pow(res.Dn50, 3);
            std::stringstream ss;
            ss << std::left << std::setw(30) << res.name 
               << " | " << std::left << std::setw(8) << std::fixed << std::setprecision(4) << res.Ns 
               << " | " << std::left << std::setw(10) << std::fixed << std::setprecision(3) << res.Dn50 
               << " | " << std::left << std::setw(10) << std::fixed << std::setprecision(0) << mass 
               << " | " << std::left << std::setw(8) << std::fixed << std::setprecision(2) << res.Kd 
               << " | " << res.note;
            log(ss.str());
        }
        
        log("\nJUSTIFICATION & ANALYSIS:");
        for (const auto& line : report.justification) log(line);
        log(std::string(95, '-'));

        return report;
    }

    void finalize_design(FullReport& report) {
        double target_dn = report.recommended.Dn50;
        double target_mass = defaults.rho_r * std::pow(target_dn, 3);
        
        report.armor_layer = design_layer(target_mass, target_dn, true);
        
        // --- Report Armor Layer ---
        log("\n" + std::string(95, '='));
        if (defaults.use_en13383) log("   4. ROCK ARMOUR LAYER DESIGN (EN 13383 Standard)");
        else log("   4. ROCK ARMOUR LAYER DESIGN (Custom Grading)");
        log(std::string(95, '='));
        log("PRIMARY ARMOR LAYER");
        
        auto fmt_line = [](std::string lab, std::string val) {
            std::stringstream ss;
            ss << "   " << std::left << std::setw(36) << lab << ": " << val;
            return ss.str();
        };

        std::stringstream ss_tw, ss_tm, ss_td;
        ss_tw << std::fixed << std::setprecision(2) << report.armor_layer.target_W_kN << " kN";
        log("   Theoretical Required W    : " + ss_tw.str());
        
        ss_tm << std::fixed << std::setprecision(0) << report.armor_layer.target_M50_kg << " kg";
        log("   Theoretical Required M50  : " + ss_tm.str());
        
        ss_td << std::fixed << std::setprecision(3) << report.armor_layer.target_Dn_m << " m";
        log("   Theoretical Required Dn50 : " + ss_td.str());
        
        log(std::string(40, '-'));

        if (!report.armor_layer.design_valid && defaults.use_en13383) {
             log("   [WARNING] No standard EN13383 grading found for this mass.");
             // Don't return, allow printing empty/0 values as per python logic if desired, 
             // but python returns early. Let's return to match python.
             return; 
        } else {
             log(fmt_line("Adopted rock grading", report.armor_layer.grading_name));

             std::stringstream ss_min; 
             ss_min << std::fixed << std::setprecision(2) << report.armor_layer.w_min_kn << " kN (" 
                    << std::fixed << std::setprecision(0) << report.armor_layer.w_min_kg << " kg)";
             log(fmt_line("Grading Min (Lower Limit)", ss_min.str()));

             std::stringstream ss_max; 
             ss_max << std::fixed << std::setprecision(2) << report.armor_layer.w_max_kn << " kN (" 
                    << std::fixed << std::setprecision(0) << report.armor_layer.w_max_kg << " kg)";
             log(fmt_line("Grading Max (Upper Limit)", ss_max.str()));

             std::stringstream ss_rep;
             ss_rep << std::fixed << std::setprecision(0) << report.armor_layer.m_mean_kg << " kg";
             log(fmt_line("Representative M50", ss_rep.str()));

             std::stringstream ss_dn;
             ss_dn << std::fixed << std::setprecision(3) << report.armor_layer.actual_dn << " m";
             log(fmt_line("Nominal Diameter (Dn_rock)", ss_dn.str()));

             std::stringstream ss_thick;
             ss_thick << std::fixed << std::setprecision(2) << report.armor_layer.thickness << " m";
             log(fmt_line("Double Layer Thickness", ss_thick.str()));

             std::stringstream ss_pack; 
             ss_pack << std::fixed << std::setprecision(2) << report.armor_layer.packing_density;
             log(fmt_line("Packing Density [rocks/100m2]", ss_pack.str()));
        }
        log(std::string(95, '-'));

        // --- Underlayer ---
        double target_mass_ul = report.armor_layer.m_mean_kg / 10.0;
        double target_dn_ul = std::pow(target_mass_ul / defaults.rho_r, 1.0/3.0);
        report.underlayer = design_layer(target_mass_ul, target_dn_ul, false);

        log("UNDERLAYER (FILTER LAYER)");
        std::stringstream ss_targ_ul_w, ss_targ_ul_m;
        ss_targ_ul_w << std::fixed << std::setprecision(3) << report.underlayer.target_W_kN << " kN";
        log("   Target Weight (M50 / 10)  : " + ss_targ_ul_w.str());
        
        ss_targ_ul_m << std::fixed << std::setprecision(1) << report.underlayer.target_M50_kg << " kg";
        log("   Target Mass (M50 / 10)    : " + ss_targ_ul_m.str());
        log(std::string(40, '-'));

        if (!report.underlayer.design_valid && defaults.use_en13383) {
             log("   [WARNING] No suitable standard underlayer grading found.");
        } else {
             log(fmt_line("Adopted rock grading", report.underlayer.grading_name));

             std::stringstream ss_min; 
             ss_min << std::fixed << std::setprecision(2) << report.underlayer.w_min_kn << " kN (" 
                    << std::fixed << std::setprecision(0) << report.underlayer.w_min_kg << " kg)";
             log(fmt_line("Grading Min (Lower Limit)", ss_min.str()));

             std::stringstream ss_max; 
             ss_max << std::fixed << std::setprecision(2) << report.underlayer.w_max_kn << " kN (" 
                    << std::fixed << std::setprecision(0) << report.underlayer.w_max_kg << " kg)";
             log(fmt_line("Grading Max (Upper Limit)", ss_max.str()));

             std::stringstream ss_rep;
             ss_rep << std::fixed << std::setprecision(1) << report.underlayer.m_mean_kg << " kg";
             log(fmt_line("Representative M50", ss_rep.str()));

             std::stringstream ss_dn;
             ss_dn << std::fixed << std::setprecision(3) << report.underlayer.actual_dn << " m";
             log(fmt_line("Nominal Diameter (Dn_rock)", ss_dn.str()));

             std::stringstream ss_thick;
             ss_thick << std::fixed << std::setprecision(2) << report.underlayer.thickness << " m";
             log(fmt_line("Double Layer Thickness", ss_thick.str()));

             std::stringstream ss_pack; 
             ss_pack << std::fixed << std::setprecision(2) << report.underlayer.packing_density;
             log(fmt_line("Packing Density [rocks/100m2]", ss_pack.str()));
        }
        log(std::string(95, '='));
    }

    void save_report(std::string filepath="output.txt") {
        std::ofstream outfile(filepath);
        if (outfile.is_open()) {
            for (const auto& line : log_buffer) outfile << line << "\n";
            outfile.close();
            // GUI does not print to cout here, handled by messagebox or status
        } else {
            std::cerr << "Error writing file." << std::endl;
        }
    }
};

// ==============================================================================
// GUI IMPLEMENTATION
// ==============================================================================

static std::wstring fix_newlines_for_edit_control(const std::string &text) {
    std::wstring out;
    out.reserve(text.size() + 100);
    for (char c : text) {
        if (c == '\n') {
            out.push_back(L'\r');
        }
        out.push_back(static_cast<wchar_t>(c));
    }
    return out;
}

// --- FILE OUTPUT HELPER ---
void SaveReportToFile(const std::string& report_content) {
    // Overwrite "output.txt"
    std::ofstream file("output.txt");
    if (file.is_open()) {
        file << report_content; 
        file.close();
    }
}

// --- Win32 GUI Specifics ---

#define IDC_EDIT_HS 101
#define IDC_EDIT_TM 102
#define IDC_EDIT_HTOE 103
#define IDC_EDIT_SLOPE 104
#define IDC_EDIT_FORE 105
#define IDC_EDIT_RHO_R 106
#define IDC_EDIT_P 107
#define IDC_EDIT_DRATIO 108
#define IDC_EDIT_SD 109
#define IDC_EDIT_DUR 110
#define IDC_CHK_EN13383 111
#define IDC_COMBO_FORMULA 112
#define IDC_BUTTON_COMPUTE 113
#define IDC_OUTPUT 114

HWND hEditHs, hEditTm, hEditHToe, hEditSlope, hEditFore, hEditRhoR, hEditP, hEditDRatio, hEditSd, hEditDur;
HWND hChkEN13383, hComboFormula, hOutput;

RockSlopeCalculator calculator;
HFONT hMonoFont = NULL;
HFONT hUIFont = NULL;

void CreateControls(HWND hwnd) {
    hUIFont = CreateFontW(22, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                          DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                          DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");

    int x_label = 10;
    int x_input = 220;
    int w_label = 200;
    int w_input = 100;
    int y = 10;
    int step = 38;

    auto fmt = [](double val, int prec=2) {
        std::wstringstream ss;
        ss << std::fixed << std::setprecision(prec) << val;
        return ss.str();
    };

    auto CreateLabel = [&](const wchar_t* text, int x, int y) {
        HWND h = CreateWindowW(L"STATIC", text, WS_CHILD | WS_VISIBLE, x, y, w_label, 20, hwnd, NULL, NULL, NULL);
        SendMessageW(h, WM_SETFONT, (WPARAM)hUIFont, TRUE);
    };

    auto CreateInput = [&](const wchar_t* def, int id, int y) -> HWND {
        HWND h = CreateWindowW(L"EDIT", def, WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL, 
                              x_input, y, w_input, 25, hwnd, (HMENU)(INT_PTR)id, NULL, NULL);
        SendMessageW(h, WM_SETFONT, (WPARAM)hUIFont, TRUE);
        return h;
    };

    // Hs
    CreateLabel(L"Hs (Significant Wave) [m]:", x_label, y);
    hEditHs = CreateInput(fmt(calculator.defaults.Hs).c_str(), IDC_EDIT_HS, y);
    y += step;

    // Tm
    CreateLabel(L"Tm-1,0 (Spectral Period) [s]:", x_label, y);
    hEditTm = CreateInput(fmt(calculator.defaults.Tm10).c_str(), IDC_EDIT_TM, y);
    y += step;

    // h_toe
    CreateLabel(L"h_toe (Depth at Toe) [m]:", x_label, y);
    hEditHToe = CreateInput(fmt(calculator.defaults.h_toe).c_str(), IDC_EDIT_HTOE, y);
    y += step;

    // Slope
    CreateLabel(L"Structure Slope (cot alpha):", x_label, y);
    hEditSlope = CreateInput(fmt(calculator.defaults.slope_m).c_str(), IDC_EDIT_SLOPE, y);
    y += step;

    // Foreshore
    CreateLabel(L"Foreshore Slope (cot beta):", x_label, y);
    hEditFore = CreateInput(fmt(calculator.defaults.foreshore_m).c_str(), IDC_EDIT_FORE, y);
    y += step;

    // Rho R
    CreateLabel(L"Rock Density [kg/m3]:", x_label, y);
    hEditRhoR = CreateInput(fmt(calculator.defaults.rho_r, 0).c_str(), IDC_EDIT_RHO_R, y);
    y += step;

    // P
    CreateLabel(L"Permeability P (-):", x_label, y);
    hEditP = CreateInput(fmt(calculator.defaults.P).c_str(), IDC_EDIT_P, y);
    y += step;

    // D Ratio
    CreateLabel(L"Core/Armor D Ratio (-):", x_label, y);
    hEditDRatio = CreateInput(fmt(calculator.defaults.D_ratio).c_str(), IDC_EDIT_DRATIO, y);
    y += step;

    // Sd
    CreateLabel(L"Damage Level Sd (-):", x_label, y);
    hEditSd = CreateInput(fmt(calculator.defaults.Sd).c_str(), IDC_EDIT_SD, y);
    y += step;

    // Duration
    CreateLabel(L"Storm Duration [h]:", x_label, y);
    hEditDur = CreateInput(fmt(calculator.defaults.duration).c_str(), IDC_EDIT_DUR, y);
    y += step;

    // EN13383 Checkbox
    hChkEN13383 = CreateWindowW(L"BUTTON", L"Use EN13383 Standard Grading", WS_VISIBLE | WS_CHILD | BS_AUTOCHECKBOX,
                 x_label, y, 300, 20, hwnd, (HMENU)IDC_CHK_EN13383, NULL, NULL);
    SendMessageW(hChkEN13383, WM_SETFONT, (WPARAM)hUIFont, TRUE);
    CheckDlgButton(hwnd, IDC_CHK_EN13383, calculator.defaults.use_en13383 ? BST_CHECKED : BST_UNCHECKED);
    y += step + 10;

    // Formula Selection Label
    CreateLabel(L"Formula Preference:", x_label, y);
    y += 25;

    // Formula ComboBox
    hComboFormula = CreateWindowW(L"COMBOBOX", NULL, 
        WS_CHILD | WS_VISIBLE | CBS_DROPDOWNLIST | WS_VSCROLL, 
        x_label, y, 320, 300, hwnd, (HMENU)IDC_COMBO_FORMULA, NULL, NULL);
    SendMessageW(hComboFormula, WM_SETFONT, (WPARAM)hUIFont, TRUE);

    // --- COMBOBOX ORDER (Must Match FormulaResult Vector Order) ---
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Automatic (Hydraulic Zone Logic)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Hudson (1959)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Van der Meer (2021)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Van Gent Modified (2003)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Van Gent Simplified (2003)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Eldrup & Andersen (2019)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Etemad-Shahidi (2020)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Scaravaglione (Mod. VG 2025)");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"Scaravaglione (Mod. ES 2025)");
    SendMessageW(hComboFormula, CB_SETCURSEL, 0, 0); 
    
    y += 40;

    // Compute Button
    HWND hBtn = CreateWindowW(L"BUTTON", L"Calculate", WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
                 x_label, y, 200, 40, hwnd, (HMENU)IDC_BUTTON_COMPUTE, NULL, NULL);
    SendMessageW(hBtn, WM_SETFONT, (WPARAM)hUIFont, TRUE);

    // Output Box
    hOutput = CreateWindowW(L"EDIT", L"", WS_CHILD | WS_VISIBLE | WS_BORDER |
                           ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL | ES_READONLY | WS_HSCROLL | ES_AUTOHSCROLL,
                           350, 10, 1150, 740, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);

    hMonoFont = CreateFontW(20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                           DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                           DEFAULT_QUALITY, FIXED_PITCH | FF_DONTCARE, L"Courier New");
    SendMessageW(hOutput, WM_SETFONT, (WPARAM)hMonoFont, TRUE);
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_CREATE:
        CreateControls(hwnd);
        break;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_BUTTON_COMPUTE) {
            wchar_t buffer[64];
            Inputs inputs; 

            // Parse Inputs
            GetWindowTextW(hEditHs, buffer, 63); inputs.Hs = _wtof(buffer);
            GetWindowTextW(hEditTm, buffer, 63); inputs.Tm10 = _wtof(buffer);
            GetWindowTextW(hEditHToe, buffer, 63); inputs.h_toe = _wtof(buffer);
            GetWindowTextW(hEditSlope, buffer, 63); inputs.slope_m = _wtof(buffer);
            GetWindowTextW(hEditFore, buffer, 63); inputs.foreshore_m = _wtof(buffer);
            GetWindowTextW(hEditRhoR, buffer, 63); inputs.rho_r = _wtof(buffer);
            GetWindowTextW(hEditP, buffer, 63); inputs.P = _wtof(buffer);
            GetWindowTextW(hEditDRatio, buffer, 63); inputs.D_ratio = _wtof(buffer);
            GetWindowTextW(hEditSd, buffer, 63); inputs.Sd = _wtof(buffer);
            GetWindowTextW(hEditDur, buffer, 63); inputs.duration = _wtof(buffer);
            
            // Hardcoded constant defaults
            inputs.rho_w = 1025.0; 
            
            // Checkbox
            inputs.use_en13383 = (IsDlgButtonChecked(hwnd, IDC_CHK_EN13383) == BST_CHECKED);
            
            // ComboBox (Override Index)
            int selIndex = (int)SendMessageW(hComboFormula, CB_GETCURSEL, 0, 0);
            inputs.formula_override_index = selIndex; 

            if (inputs.Hs <= 0 || inputs.Tm10 <= 0 || inputs.h_toe <= 0) {
                 MessageBoxW(hwnd, L"Please enter positive values for basic hydraulic parameters.", L"Input Error", MB_ICONERROR | MB_OK);
                 break;
            }

            try {
                // Run Calculation
                FullReport r = calculator.solve(inputs);
                
                // GUI MANUAL SELECTION REPLACEMENT
                if (inputs.formula_override_index > 0) {
                    std::vector<FormulaResult> valid_results;
                    for (const auto& res : r.comparison) {
                        if (res.valid) valid_results.push_back(res);
                    }
                    
                    // Check bounds (Index 1-based from combobox)
                    if (inputs.formula_override_index <= valid_results.size()) {
                        r.recommended = valid_results[inputs.formula_override_index - 1];
                        std::string msg = "\n[MANUAL OVERRIDE] User switched selection to: ";
                        std::string uname = r.recommended.name;
                        std::transform(uname.begin(), uname.end(), uname.begin(), ::toupper);
                        msg += uname;
                        calculator.log(msg);
                    }
                }
                
                calculator.finalize_design(r);

                // Get Output text from the calculator class buffer
                std::string narrowReport = calculator.get_gui_output();
                
                // Save to file output.txt
                calculator.save_report("output.txt");

                // Display in GUI
                std::wstring guiText = fix_newlines_for_edit_control(narrowReport);
                SetWindowTextW(hOutput, guiText.c_str());
            } catch (const std::exception& e) {
                const char* what = e.what();
                std::wstring errMsg(what, what + strlen(what));
                MessageBoxW(hwnd, errMsg.c_str(), L"Calculation Error", MB_ICONERROR | MB_OK);
            }
        }
        break;

    case WM_DESTROY:
        if (hMonoFont) DeleteObject(hMonoFont);
        if (hUIFont) DeleteObject(hUIFont);
        PostQuitMessage(0);
        break;

    default:
        return DefWindowProcW(hwnd, msg, wParam, lParam);
    }
    return 0;
}

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow) {
    const wchar_t CLASS_NAME[] = L"RockSlopeCalcWindow";

    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(WNDCLASSEXW);
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.hIconSm = LoadIcon(NULL, IDI_APPLICATION);

    if (!RegisterClassExW(&wc)) {
        return 0;
    }

    HWND hwnd = CreateWindowExW(
        0, CLASS_NAME, L"Rock Slope Stability Calculator",
        WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME & ~WS_MAXIMIZEBOX,
        CW_USEDEFAULT, CW_USEDEFAULT, 1500, 750,
        NULL, NULL, hInstance, NULL);

    if (!hwnd) {
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    MSG msg = {};
    while (GetMessageW(&msg, NULL, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }

    return static_cast<int>(msg.wParam);
}