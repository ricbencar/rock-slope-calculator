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
//
//    1. U.S. Army Corps of Engineers. (1984). "Shore Protection Manual." Vol. I & II. 
//       Coastal Engineering Research Center, Vicksburg, MS. 
//       Link: https://usace.contentdm.oclc.org/digital/collection/p16021coll11/id/1934/ 
//
//    2. Van der Meer, J. W. (1987). "Stability of breakwater armour layers—Design formulae." 
//       Coastal Engineering, 11(3), 219-239.
//       Link: https://doi.org/10.1016/0378-3839(87)90013-5
//
//    3. Van der Meer, J. W. (1988). "Rock Slopes and Gravel Beaches Under Wave Attack." 
//       Doctoral Thesis, Delft University of Technology.
//       [Foundational work for modern stability formulas]
//       Link: https://repository.tudelft.nl/islandora/object/uuid:404b5nec-hm
//
//    4. Van Gent, M. R. A. (1995). "Wave interaction with permeable coastal structures."
//       Doctoral Thesis, Delft University of Technology.
//       Link: https://repository.tudelft.nl/islandora/object/uuid:7bbff8e4-215d-4bfc-a3af-51cdecb754bd
//
//    5. U.S. Army Corps of Engineers (USACE). (2002). "Coastal Engineering Manual." 
//       Engineer Manual 1110-2-1100, Washington, D.C.
//       Link: https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/
//
//    6. Van Gent, M. R. A., Smale, A. J., & Kuiper, C. (2003). "Stability of rock slopes 
//       with shallow foreshores." Proceedings of Coastal Structures 2003, Portland, OR, 100-112.
//       Link: https://doi.org/10.1061/40733(147)9
//
//    7. Van Gent, M. R. A. (2004). "On the stability of rock slopes." 
//       Environmentally Friendly Coastal Protection: Proceedings of the NATO Advanced 
//       Research Workshop, Varna, Bulgaria.
//       Link: https://doi.org/10.1007/1-4020-3301-X_12
//
//    8. CIRIA, CUR, CETMEF. (2007). "The Rock Manual. The Use of Rock in Hydraulic Engineering." 
//       (2nd edition). C683, CIRIA, London.
//       Link: https://www.ciria.org/ItemDetail?iProductCode=C683
//
//    9. CEN (2013). "EN 13383-1:2013 Armourstone - Part 1: Specification."
//       European Committee for Standardization.
//       Link: https://standards.iteh.ai/catalog/standards/cen/5f6b770f-38ba-4320-ad1a-83d0e29f7db2/en-13383-1-2013
//
//   10. Eldrup, M. R., & Lykke Andersen, T. (2019). "Extension of shallow water rock armour 
//       stability formulae to nonlinear waves." Coastal Engineering, 153, 103536.
//       Link: https://doi.org/10.1016/j.coastaleng.2019.103536
//
//   11. Etemad-Shahidi, A., Bali, M., & Van Gent, M. R. A. (2020). "On the stability of rock 
//       armored rubble mound structures." Coastal Engineering, 158, 103655.
//       Link: https://doi.org/10.1016/j.coastaleng.2020.103655
//
//   12. Van der Meer, J. W. (2021). "Rock armour slope stability under wave attack; the 
//       Van der Meer Formula revisited." Journal of Coastal and Hydraulic Structures, 1.
//       Link: https://doi.org/10.48438/jchs.2021.0008
//
//   13. Van der Meer, J. W., Lykke Andersen, T., & Roge Eldrup, M. (2024). "Rock Armour 
//       Slope Stability under Wave Attack in Shallow Water." Journal of Coastal and 
//       Hydraulic Structures, 4.
//       Link: https://doi.org/10.59490/jchs.2024.0035
//
//   14. Scaravaglione, G., Marino, S., Francone, A., Leone, E., Damiani, L., Tomasicchio, 
//       G. R., Van Gent, M. R. A., & Saponieri, A. (2025). "The influence of shallow water 
//       on rock armour stability." Coastal Engineering, 197, 104657.
//       Link: https://doi.org/10.1016/j.coastaleng.2024.104657
//
// 5. COMPILATION INSTRUCTIONS (MinGW on Windows):
//
//    g++ -O3 -std=c++17 rock_slope_calculator_cli.cpp \
//    -o rock_slope_calculator_cli.exe -static -static-libgcc -static-libstdc++
//
// 6. EXECUTION:
//
//    Interactive Mode:
//      rock_slope_calculator_cli.exe
//
//    Command Line Mode:
//      rock_slope_calculator_cli.exe [Hs] [Tm] [h_toe] [slope] [fore] [rho_r] [P] [D_ratio] [Sd] [Duration] [UseEN13383]
//
//    Example:
//      rock_slope_calculator_cli.exe 2.5 10.0 6.0 2.0 30.0 2650 0.4 0.3 2.0 6.0 true
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
    double NLL;      // Nominal Lower Limit (kg)
    double NUL;      // Nominal Upper Limit (kg)
    double M50;      // Representative M50 = 0.5 * (NLL + NUL)
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
    bool manual_selection; // Interactive manual selection flag
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

        // Initialize EN 13383 Database 
        // Logic: NLL and NUL taken from Category Name for LMA/HMA.
        // M50 calculated as 0.5 * (NLL + NUL).
        standard_gradings = {
            // Name            NLL (kg)   NUL (kg)   M50 (Calc)
            {"CP 45/125",      0.4,       1.2,       0.8},
            {"CP 63/180",      1.2,       3.8,       2.5},
            {"CP 90/250",      3.1,       9.3,       6.2},
            {"CP 45/180",      0.4,       1.2,       0.8}, // Kept small range conservative
            {"CP 90/180",      2.1,       2.8,       2.45},
            
            // Light Mass Armour (LMA) - NLL/NUL from name
            {"LMA 5-40",       5.0,       40.0,      22.5},
            {"LMA 10-60",      10.0,      60.0,      35.0},
            {"LMA 15-120",     15.0,      120.0,     67.5},
            {"LMA 40-200",     40.0,      200.0,     120.0},
            {"LMA 60-300",     60.0,      300.0,     180.0},
            {"LMA 15-300",     15.0,      300.0,     157.5},
            
            // Heavy Mass Armour (HMA) - NLL/NUL from name
            {"HMA 300-1000",   300.0,     1000.0,    650.0},
            {"HMA 1000-3000",  1000.0,    3000.0,    2000.0},
            {"HMA 3000-6000",  3000.0,    6000.0,    4500.0},
            {"HMA 6000-10000", 6000.0,    10000.0,   8000.0},
            {"HMA 10000-15000",10000.0,   15000.0,   12500.0}
        };

        // Sort by M50 for selection logic
        std::sort(standard_gradings.begin(), standard_gradings.end(), 
            [](const GradingDef& a, const GradingDef& b) {
                return a.M50 < b.M50;
            });
    }

    void log(std::string message) {
        std::cout << message << std::endl;
        log_buffer.push_back(message);
    }

    // --- Wave Mechanics ---

    double solve_wavelength(double T, double h) {
        if (h <= 0) return 0.0;
        double L0 = G * T * T / (2 * M_PI);
        double k0h = 2 * M_PI * h / L0;
        
        // Initial Guess (Carvalho 2006)
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
        
        // Celerity
        if (in.Tm10 > 0) h.C = h.L_toe / in.Tm10;
        else h.C = 0;
        
        // Group Celerity
        double k = (h.L_toe > 0) ? (2 * M_PI / h.L_toe) : 0;
        double kh = k * in.h_toe;
        double n;
        if (kh > 20) n = 0.5;
        else if (kh <= 0) n = 1.0;
        else n = 0.5 * (1 + (2 * kh) / std::sinh(2 * kh));
        h.Cg = n * h.C;

        // Steepness & Surf Similarity
        h.s_m10 = in.Hs / h.L0;
        h.s_local = (h.L_toe > 0) ? (in.Hs / h.L_toe) : 0;
        
        if (h.s_m10 > 0) 
            h.xi_m10 = std::tan(dp.alpha_rad) / std::sqrt(h.s_m10);
        else 
            h.xi_m10 = 0;

        // Relative Depth & Breaker Type
        h.rel_depth = (in.Hs > 0) ? (in.h_toe / in.Hs) : 999.0;

        if (h.xi_m10 < 0.5) h.breaker_type = "Spilling";
        else if (h.xi_m10 < 1.8) h.breaker_type = "Plunging";
        else if (h.xi_m10 < 3.0) h.breaker_type = "Surging";
        else h.breaker_type = "Collapsing/Surging";

        // Zone Classification (Output Logic based on Relative Depth)
        if (h.rel_depth > 3.0) h.zone_desc = "ZONE 1: Deep to Intermediate (h/Hm0 > 3.0)";
        else if (h.rel_depth > 1.5) h.zone_desc = "ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0)";
        else if (h.rel_depth > 0.5) h.zone_desc = "ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5)";
        else h.zone_desc = "ZONE 4: Extremely Shallow Water (h/Hm0 <= 0.5)";

        return h;
    }

    // --- Stability Formulas ---

    FormulaResult calc_hudson(const Inputs& p, const DerivedParams& dp, const Hydraulics& h) {
        // Hudson Kd based on Zone Classification depending on Relative Depth
        // Matching Python logic
        double Kd;
        if (h.rel_depth > 3.0) Kd = 4.0;
        else if (h.rel_depth > 1.5) Kd = 3.5;
        else if (h.rel_depth > 0.5) Kd = 3.0;
        else Kd = 2.0;

        double Dn = (1.27 * p.Hs) / (dp.Delta * std::pow(Kd * dp.cot_alpha, 1.0/3.0));
        double Ns = (1.27 * p.Hs) / (dp.Delta * Dn); // Re-calculate Ns for consistency
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
        // H2%/Hs ratio estimation based on relative depth (matching Python logic)
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

        // Helper to find results safely
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
        // FormulaResult eldrup = get_res("Eldrup");

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
            
            // Van der Meer Analysis
            log.push_back("   **A. Van der Meer (2021 Rewritten) [RECOMMENDED]**");
            log.push_back("      * **Advantages:** This formula is the modernized industry standard.");
            log.push_back("        Van der Meer (2021) rewrote the original formula to use the spectral period (Tm-1,0),");
            log.push_back("        eliminating the influence of spectral shape.");
            log.push_back("        Van der Meer et al. (2024) confirmed its validity for h/Hm0 > 1.5, preferring Hm0 over H1/3 for nonlinear waves.");
            log.push_back("      * **Physics:** It correctly assumes a Rayleigh distribution of wave heights, aligning with the actual deep-water statistics.");

            // Van Gent Analysis
            log.push_back("   **B. Van Gent Modified (2003)**");
            log.push_back("      * **Context:** This formula incorporates the ratio H2%/Hs.");
            log.push_back("        In deep water, with H2%/Hs = 1.4, this formula essentially converges closely with the Van der Meer predictions.");
            log.push_back("        However, its specific calibration was focused on the effects of shallow foreshores.");

            // Etemad-Shahidi Analysis
            log.push_back("   **C. Etemad-Shahidi (2020)**");
            log.push_back("      * **Comparison:** Etemad-Shahidi (2020) provides a robust formula validated for both deep and shallow water.");
            log.push_back("        It introduces a physical permeability parameter (D_core/D_armor) to replace the nominal P factor,");
            log.push_back("        reducing uncertainty. However, Van der Meer remains the primary standard for deep water.");

            // Divergence Check
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
            
            // Van der Meer Analysis
            log.push_back("   **A. Van der Meer (2021)**");
            log.push_back("      * **Advantages:** Van der Meer et al. (2024) extensively re-analyzed shallow water data and concluded");
            log.push_back("        that the rewritten Van der Meer formula (using Tm-1,0) is valid down to h/Hm0 = 1.5.");
            log.push_back("        It performs reasonably well, with slightly less reliability in the 1.0 < h/Hm0 < 1.5 range.");
            log.push_back("      * **Note:** For nonlinear waves in this zone, using Hm0 is preferred over H1/3 for nonlinear waves to avoid deviations.");

            // Van Gent Mod Analysis
            log.push_back("   **B. Van Gent Modified (2003)**");
            log.push_back("      * **Constraint:** This formula explicitly relies on the ratio H2%/Hs. Research by Van der Meer et al. (2024)");
            log.push_back("        highlights that predicting H2% accurately in this transition zone (where the ratio dips to ~1.2)");
            log.push_back("        is notoriously inaccurate without physical modeling. The formula is valid, but the input uncertainty is high.");

            // Van Gent Simp Analysis
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
        // Selection Rule: NLL < Target M50 < NUL
        // Tie-breaker: Choose class with smaller (NUL - NLL)
        
        GradingDef selected = {"", 0, 0, 0};
        bool found = false;
        double min_range_width = std::numeric_limits<double>::max();
        
        for (const auto& g : standard_gradings) {
            // Check containment: Target must be strictly inside nominal limits
            if (target_mass > g.NLL && target_mass < g.NUL) {
                
                double current_range = g.NUL - g.NLL;
                
                // Update if this is the first match OR if this range is tighter (smaller)
                if (current_range < min_range_width) {
                    min_range_width = current_range;
                    selected = g;
                    found = true;
                }
            }
        }

        if (found) {
            ld.grading_name = selected.name;
            ld.m_mean_kg = selected.M50;
            
            // Store Nominal Limits for report
            ld.w_min_kg = selected.NLL; 
            ld.w_max_kg = selected.NUL;
            
            // Recalculate weights in kN for internal consistency
            ld.w_min_kn = selected.NLL * G / 1000.0;
            ld.w_max_kn = selected.NUL * G / 1000.0;
            ld.w_mean_kn = ld.m_mean_kg * G / 1000.0;
            
            ld.actual_dn = std::pow(ld.w_mean_kn / gamma_r, 1.0/3.0);
            ld.design_valid = true;
        } else {
            ld.grading_name = "No Standard Fit (Target outside NLL-NUL)";
            ld.design_valid = false;
            // Note: design_valid = false triggers the custom grading fallback below
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

    void handle_manual_selection(FullReport& report) {
        if (!report.inputs.manual_selection) return;

        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "   MANUAL SELECTION MODE" << std::endl;
        std::cout << std::string(80, '=') << std::endl;
        
        std::vector<FormulaResult> valid_results;
        for (const auto& res : report.comparison) {
            if (res.valid) valid_results.push_back(res);
        }
        
        for (size_t i = 0; i < valid_results.size(); ++i) {
            std::cout << "   [" << (i + 1) << "] " << valid_results[i].name << std::endl;
        }
        
        bool valid_choice = false;
        while (!valid_choice) {
            std::cout << "\nEnter the number of your preferred formula (1-" << valid_results.size() << "): ";
            std::string input_str;
            std::getline(std::cin, input_str);
            try {
                size_t choice = std::stoul(input_str);
                if (choice >= 1 && choice <= valid_results.size()) {
                    report.recommended = valid_results[choice - 1];
                    std::string msg = "\n[MANUAL OVERRIDE] User switched selection to: ";
                    // Convert to upper case for log display
                    std::string uname = report.recommended.name;
                    std::transform(uname.begin(), uname.end(), uname.begin(), ::toupper);
                    msg += uname;
                    
                    log(msg);
                    valid_choice = true;
                } else {
                    std::cout << "Invalid selection. Try again." << std::endl;
                }
            } catch (...) {
                std::cout << "Please enter a valid number." << std::endl;
            }
        }
    }

    void finalize_design(FullReport& report) {
        // 1. Design Primary Armor Layer
        double target_dn = report.recommended.Dn50;
        double target_mass = defaults.rho_r * std::pow(target_dn, 3);
        
        report.armor_layer = design_layer(target_mass, target_dn, true);
        
        // --- Report Armor Layer Header ---
        log("\n" + std::string(95, '='));
        if (defaults.use_en13383) log("   4. ROCK ARMOUR LAYER DESIGN (EN 13383 Standard)");
        else log("   4. ROCK ARMOUR LAYER DESIGN (Custom Grading)");
        log(std::string(95, '='));
        log("PRIMARY ARMOR LAYER");
        
        // Helper for formatting lines
        auto fmt_line = [](std::string lab, std::string val) {
            std::stringstream ss;
            ss << "   " << std::left << std::setw(36) << lab << ": " << val;
            return ss.str();
        };

        // Print Theoretical Requirements (Always visible)
        std::stringstream ss_tw, ss_tm, ss_td;
        ss_tw << std::fixed << std::setprecision(2) << report.armor_layer.target_W_kN << " kN";
        log("   Theoretical Required W    : " + ss_tw.str());
        
        ss_tm << std::fixed << std::setprecision(0) << report.armor_layer.target_M50_kg << " kg";
        log("   Theoretical Required M50  : " + ss_tm.str());
        
        ss_td << std::fixed << std::setprecision(3) << report.armor_layer.target_Dn_m << " m";
        log("   Theoretical Required Dn50 : " + ss_td.str());
        
        log(std::string(40, '-'));

        // --- Armor Layer Details ---
        if (!report.armor_layer.design_valid && defaults.use_en13383) {
             log("   [WARNING] No standard EN13383 grading found for this mass.");
             // We do not return here to ensure Underlayer is still calculated/attempted
        } else {
             log(fmt_line("Adopted rock grading", report.armor_layer.grading_name));

             if (defaults.use_en13383) {
                 // === EN 13383 SPECIFIC FORMAT ===
                 double NLL = report.armor_layer.w_min_kg; // Stored in w_min slot
                 double NUL = report.armor_layer.w_max_kg; // Stored in w_max slot
                 
                 double ELL = 0.7 * NLL;
                 double EUL = 1.5 * NUL;
                 double rep_M50 = 0.5 * (NLL + NUL);

                 // Representative M50
                 std::stringstream ss_rep;
                 ss_rep << std::fixed << std::setprecision(1) << rep_M50 << " kg";
                 log(fmt_line("Representative M50", ss_rep.str()));

                 // Nominal Limits
                 std::stringstream ss_nll, ss_nul; 
                 ss_nll << std::fixed << std::setprecision(1) << NLL << " kg";
                 log(fmt_line("Nominal lower limit (NLL)", ss_nll.str()));

                 ss_nul << std::fixed << std::setprecision(1) << NUL << " kg";
                 log(fmt_line("Nominal upper limit (NUL)", ss_nul.str()));

                 // Extreme Limits
                 std::stringstream ss_ell, ss_eul; 
                 ss_ell << std::fixed << std::setprecision(1) << ELL << " kg";
                 log(fmt_line("Extreme lower limit (ELL)", ss_ell.str()));

                 ss_eul << std::fixed << std::setprecision(1) << EUL << " kg";
                 log(fmt_line("Extreme upper limit (EUL)", ss_eul.str()));

             } else {
                 // === CUSTOM GRADING FORMAT (Original) ===
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
             }

             // Common Geometry Outputs
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

        // 2. Design Underlayer
        // Rule: M50_underlayer = M50_armor / 10
        double ref_mass_armor = (defaults.use_en13383) 
            ? 0.5 * (report.armor_layer.w_min_kg + report.armor_layer.w_max_kg) // Recalculate Rep M50
            : report.armor_layer.m_mean_kg;
            
        double target_mass_ul = ref_mass_armor / 10.0;
        double target_dn_ul = std::pow(target_mass_ul / defaults.rho_r, 1.0/3.0);
        report.underlayer = design_layer(target_mass_ul, target_dn_ul, false);

        log("UNDERLAYER (FILTER LAYER)");
        std::stringstream ss_targ_ul_w, ss_targ_ul_m;
        ss_targ_ul_w << std::fixed << std::setprecision(3) << report.underlayer.target_W_kN << " kN";
        log("   Target Weight (M50 / 10)  : " + ss_targ_ul_w.str());
        
        ss_targ_ul_m << std::fixed << std::setprecision(1) << report.underlayer.target_M50_kg << " kg";
        log("   Target Mass (M50 / 10)    : " + ss_targ_ul_m.str());
        log(std::string(40, '-'));

        // --- Underlayer Details ---
        if (!report.underlayer.design_valid && defaults.use_en13383) {
             log("   [WARNING] No suitable standard underlayer grading found.");
        } else {
             log(fmt_line("Adopted rock grading", report.underlayer.grading_name));

             if (defaults.use_en13383) {
                 // === EN 13383 SPECIFIC FORMAT ===
                 double NLL = report.underlayer.w_min_kg;
                 double NUL = report.underlayer.w_max_kg;
                 
                 double ELL = 0.7 * NLL;
                 double EUL = 1.5 * NUL;
                 double rep_M50 = 0.5 * (NLL + NUL);

                 std::stringstream ss_rep;
                 ss_rep << std::fixed << std::setprecision(1) << rep_M50 << " kg";
                 log(fmt_line("Representative M50", ss_rep.str()));

                 std::stringstream ss_nll, ss_nul; 
                 ss_nll << std::fixed << std::setprecision(1) << NLL << " kg";
                 log(fmt_line("Nominal lower limit (NLL)", ss_nll.str()));

                 ss_nul << std::fixed << std::setprecision(1) << NUL << " kg";
                 log(fmt_line("Nominal upper limit (NUL)", ss_nul.str()));

                 std::stringstream ss_ell, ss_eul; 
                 ss_ell << std::fixed << std::setprecision(1) << ELL << " kg";
                 log(fmt_line("Extreme lower limit (ELL)", ss_ell.str()));

                 ss_eul << std::fixed << std::setprecision(1) << EUL << " kg";
                 log(fmt_line("Extreme upper limit (EUL)", ss_eul.str()));

             } else {
                 // === CUSTOM GRADING FORMAT (Original) ===
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
             }

             // Common Geometry Outputs
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
            std::cout << "\n[System] Report saved to " <<
#if defined(_WIN32)
            // Cheap way to get absolute path approx
            "C:\\...\\" + filepath 
#else
            filepath
#endif
            << std::endl; // Matches python abs path msg loosely
        } else {
            std::cerr << "Error writing file." << std::endl;
        }
    }
};

// ==============================================================================
// 5. MAIN EXECUTION
// ==============================================================================
int main(int argc, char* argv[]) {
    RockSlopeCalculator calc;
    Inputs in = calc.defaults;

    auto get_input_str = [](std::string prompt, std::string default_val) -> std::string {
        std::cout << prompt << " [default: " << default_val << "]: ";
        std::string s;
        std::getline(std::cin, s);
        if (!s.empty()) {
            // Trim whitespace
            s.erase(0, s.find_first_not_of(" \t\n\r"));
            s.erase(s.find_last_not_of(" \t\n\r") + 1);
        }
        if (s.empty()) return default_val;
        return s;
    };

    auto get_input_double = [&](std::string prompt, double default_val) -> double {
        std::stringstream ss; ss << default_val;
        std::string s = get_input_str(prompt, ss.str());
        try { return std::stod(s); } catch(...) { return default_val; }
    };

    auto get_input_bool = [&](std::string prompt, bool default_val) -> bool {
        std::string def_str = default_val ? "True" : "False";
        std::string s = get_input_str(prompt, def_str);
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        if (s == "true" || s == "1" || s == "t" || s == "yes" || s == "y") return true;
        if (s == "false" || s == "0" || s == "f" || s == "no" || s == "n") return false;
        return default_val;
    };

    if (argc >= 10) {
        // CLI Mode
        try {
            in.Hs = std::stod(argv[1]);
            in.Tm10 = std::stod(argv[2]);
            in.h_toe = std::stod(argv[3]);
            in.slope_m = std::stod(argv[4]);
            in.foreshore_m = std::stod(argv[5]);
            in.rho_r = std::stod(argv[6]);
            in.P = std::stod(argv[7]);
            in.D_ratio = std::stod(argv[8]);
            in.Sd = std::stod(argv[9]);
            if (argc >= 11) in.duration = std::stod(argv[10]);
            if (argc >= 12) {
                 std::string s = argv[11];
                 std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
                 in.use_en13383 = (s=="true"||s=="1");
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing args: " << e.what() << "\n";
            return 1;
        }
    } else {
        // Interactive Mode
        std::cout << "\n" << std::string(80, '=') << "\n";
        std::cout << "   ROCK SLOPE STABILITY CALCULATOR\n";
        std::cout << std::string(80, '=') << "\n";
        std::cout << "Press [Enter] to accept the default value shown in brackets.\n\n";

        in.Hs = get_input_double("Significant Wave Height Hm0 (at toe) [m]", in.Hs);
        in.Tm10 = get_input_double("Spectral Period Tm-1,0 (at toe) [s]", in.Tm10);
        in.h_toe = get_input_double("Water Depth at Toe h [m]", in.h_toe);
        in.slope_m = get_input_double("Structure Slope (m:1) [cot alpha]", in.slope_m);
        in.foreshore_m = get_input_double("Foreshore Slope (m_f:1) [cot beta]", in.foreshore_m);
        in.rho_r = get_input_double("Rock Density [kg/m3]", in.rho_r);
        in.P = get_input_double("Notional Permeability P (0.4=Permeable)", in.P);
        in.D_ratio = get_input_double("Core/Armor Diameter Ratio (0.3=Typical)", in.D_ratio);
        in.Sd = get_input_double("Design Damage Level S (2.0=Start)", in.Sd);
        in.duration = get_input_double("Storm Duration [hours]", in.duration);
        in.use_en13383 = get_input_bool("Use EN13383 Standard Grading? (True/False)", in.use_en13383);
        in.manual_selection = get_input_bool("Choose stability formula instead of automatic (True/False)", false);
    }

    try {
        FullReport r = calc.solve(in);
        calc.handle_manual_selection(r);
        calc.finalize_design(r);
        calc.save_report("output.txt");
    } catch (const std::exception& e) {
        std::cerr << "Calculation Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}