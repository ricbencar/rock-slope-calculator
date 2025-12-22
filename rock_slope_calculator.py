import math
import sys
import os

# ======================================================================================
# PROGRAM DESCRIPTION & METHODOLOGY
# ======================================================================================
#
# 1. PURPOSE:
#    This software calculates the required size and weight of rock armor units for 
#    rubble mound breakwaters and revetments. It performs a comprehensive analysis 
#    of hydraulic stability across different water depth zones (Deep, Shallow, 
#    Very Shallow, and Swash zones), selecting the most appropriate empirical formula 
#    based on the hydraulic regime.
#
# 2. METHODOLOGY & LOGIC:
#    The calculator implements a "Multi-Model Consensus" approach. It computes stability 
#    using several state-of-the-art empirical formulas derived from extensive hydraulic 
#    model testing. It then evaluates the hydraulic context (Relative Depth h/Hm0) 
#    to recommend the scientifically most accurate formula.
#
#    The logic follows a 4-Step Process:
#    a. Input Processing: Parsing wave data, geometry, and material properties.
#    b. Hydraulic Analysis: Computing local wavelength, celerity, wave steepness, 
#       surf similarity (Iribarren number), and determining the breaker type.
#    c. Stability Calculation: Running multiple empirical models (Hudson, Van der Meer, 
#       Van Gent, Etemad-Shahidi, Eldrup & Andersen, Scaravaglione).
#    d. Intelligent Selection: Automatically detecting the "Hydraulic Zone" 
#       (Deep vs. Shallow vs. Swash) to recommend the most valid result.
#
# 3. HYDRAULIC ZONES & FORMULA SELECTION STRATEGY:
#    - ZONE 1: Deep/Intermediate (h/Hm0 > 3.0) -> Recommends Van der Meer (2021).
#    - ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0) -> Recommends Van der Meer (2021) 
#      or Van Gent (2003), depending on spectral shape.
#    - ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5) -> Recommends Scaravaglione (2025) 
#      (Modified ES) to account for wave breaking and horizontal stability trends.
#    - ZONE 4: Extremely Shallow/Swash (h/Hm0 <= 0.5) -> Recommends Scaravaglione (2025) 
#      (Modified VG) to handle infragravity dominance and high turbulence.
#
# 4. BIBLIOGRAPHY & SCIENTIFIC REFERENCES:
#
#    1. U.S. Army Corps of Engineers. (1984). "Shore Protection Manual." Vol. I & II. 
#       Coastal Engineering Research Center, Vicksburg, MS. 
#       Link: https://usace.contentdm.oclc.org/digital/collection/p16021coll11/id/1934/ 
#
#    2. Van der Meer, J. W. (1987). "Stability of breakwater armour layers—Design formulae." 
#       Coastal Engineering, 11(3), 219-239.
#       Link: https://doi.org/10.1016/0378-3839(87)90013-5
#
#    3. Van der Meer, J. W. (1988). "Rock Slopes and Gravel Beaches Under Wave Attack." 
#       Doctoral Thesis, Delft University of Technology.
#       [Foundational work for modern stability formulas]
#       Link: https://repository.tudelft.nl/islandora/object/uuid:404b5nec-hm
#
#    4. Van Gent, M. R. A. (1995). "Wave interaction with permeable coastal structures."
#       Doctoral Thesis, Delft University of Technology.
#       Link: https://repository.tudelft.nl/islandora/object/uuid:7bbff8e4-215d-4bfc-a3af-51cdecb754bd
#
#    5. U.S. Army Corps of Engineers (USACE). (2002). "Coastal Engineering Manual." 
#       Engineer Manual 1110-2-1100, Washington, D.C.
#       Link: https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/
#
#    6. Van Gent, M. R. A., Smale, A. J., & Kuiper, C. (2003). "Stability of rock slopes 
#       with shallow foreshores." Proceedings of Coastal Structures 2003, Portland, OR, 100-112.
#       Link: https://doi.org/10.1061/40733(147)9
#
#    7. Van Gent, M. R. A. (2004). "On the stability of rock slopes." 
#       Environmentally Friendly Coastal Protection: Proceedings of the NATO Advanced 
#       Research Workshop, Varna, Bulgaria.
#       Link: https://doi.org/10.1007/1-4020-3301-X_12
#
#    8. CIRIA, CUR, CETMEF. (2007). "The Rock Manual. The Use of Rock in Hydraulic Engineering." 
#       (2nd edition). C683, CIRIA, London.
#       Link: https://www.ciria.org/ItemDetail?iProductCode=C683
#
#    9. CEN (2013). "EN 13383-1:2013 Armourstone - Part 1: Specification."
#       European Committee for Standardization.
#       Link: https://standards.iteh.ai/catalog/standards/cen/5f6b770f-38ba-4320-ad1a-83d0e29f7db2/en-13383-1-2013
#
#   10. Eldrup, M. R., & Lykke Andersen, T. (2019). "Extension of shallow water rock armour 
#       stability formulae to nonlinear waves." Coastal Engineering, 153, 103536.
#       Link: https://doi.org/10.1016/j.coastaleng.2019.103536
#
#   11. Etemad-Shahidi, A., Bali, M., & Van Gent, M. R. A. (2020). "On the stability of rock 
#       armored rubble mound structures." Coastal Engineering, 158, 103655.
#       Link: https://doi.org/10.1016/j.coastaleng.2020.103655
#
#   12. Van der Meer, J. W. (2021). "Rock armour slope stability under wave attack; the 
#       Van der Meer Formula revisited." Journal of Coastal and Hydraulic Structures, 1.
#       Link: https://doi.org/10.48438/jchs.2021.0008
#
#   13. Van der Meer, J. W., Lykke Andersen, T., & Roge Eldrup, M. (2024). "Rock Armour 
#       Slope Stability under Wave Attack in Shallow Water." Journal of Coastal and 
#       Hydraulic Structures, 4.
#       Link: https://doi.org/10.59490/jchs.2024.0035
#
#   14. Scaravaglione, G., Marino, S., Francone, A., Leone, E., Damiani, L., Tomasicchio, 
#       G. R., Van Gent, M. R. A., & Saponieri, A. (2025). "The influence of shallow water 
#       on rock armour stability." Coastal Engineering, 197, 104657.
#       Link: https://doi.org/10.1016/j.coastaleng.2024.104657
#
# ======================================================================================

# ==============================================================================
# 1. USER CONFIGURATION (DEFAULTS)
# ==============================================================================
DEFAULT_CONFIG = {
    # --- Hydraulic Parameters ---
    'Hs': 2.0,                  # Significant Wave Height Hm0 at toe (m)
    'Tm10': 12.0,               # Spectral Period Tm-1,0 at toe (s)
    'h_toe': 7.0,               # Water Depth at Toe h (m)
    'duration': 6.0,            # Storm Duration (hours)
    
    # --- Structure Geometry ---
    'slope_m': 2.0,             # Structure Slope cot(alpha) (e.g. 2.0 for 1:2.0)
    'foreshore_m': 30.0,        # Foreshore Slope cot(beta) (e.g. 30 for 1:30)
    
    # --- Material Properties ---
    'rho_r': 2650.0,            # Rock Density (kg/m3)
    'rho_w': 1025.0,            # Water Density (kg/m3)
    'P_notional': 0.4,          # Notional Permeability P (0.1=Imp, 0.4=Typical, 0.5=Perm, 0.6=Hom)
    'D_ratio': 0.3,             # Core/Armor Diameter Ratio (Dn50_core / Dn50_armor) - Used for Cp calc

    # --- Design Criteria ---
    'Sd': 2.0,                  # Design Damage Level S (2.0 = Start of Damage)
    'grading_EN13383': True,    # Use Standard EN13383 Grading (True) or Custom Power Law (False)
    
    # --- System Settings ---
    'output_file': "output.txt"
}

# ==============================================================================
# 2. DATABASES & STANDARDS
# ==============================================================================

class RockStandards:
    """
    Database of Standard Rock Gradings according to EN 13383-1:2013.
    Values represent Mass in kg.
    Used to automatically select standard grading classes based on calculated mass.
    
    Reference: CEN (2013). "EN 13383-1:2013 Armourstone - Part 1: Specification."
    """
    GRADINGS = [
        # --- Coarse Gradings (CP) ---
        {"name": "CP 45/125",       "min": 0.4,   "max": 1.2,    "M50": 0.8},
        {"name": "CP 63/180",       "min": 1.2,   "max": 3.8,    "M50": 2.5},
        {"name": "CP 90/250",       "min": 3.1,   "max": 9.3,    "M50": 6.2},
        {"name": "CP 45/180",       "min": 0.4,   "max": 1.2,    "M50": 0.8}, 
        {"name": "CP 90/180",       "min": 2.1,   "max": 2.8,    "M50": 2.45},

        # --- Light Mass Armourstone (LMA) ---
        {"name": "LMA 5-40",        "min": 10,    "max": 20,     "M50": 15},
        {"name": "LMA 10-60",       "min": 20,    "max": 35,     "M50": 27.5},
        {"name": "LMA 15-120",      "min": 35,    "max": 60,     "M50": 47.5},
        {"name": "LMA 40-200",      "min": 80,    "max": 120,    "M50": 100},
        {"name": "LMA 60-300",      "min": 120,   "max": 190,    "M50": 155},
        {"name": "LMA 15-300",      "min": 45,    "max": 135,    "M50": 90},

        # --- Heavy Mass Armourstone (HMA) ---
        {"name": "HMA 300-1000",    "min": 540,   "max": 690,    "M50": 615},
        {"name": "HMA 1000-3000",   "min": 1700,  "max": 2100,   "M50": 1900},
        {"name": "HMA 3000-6000",   "min": 4200,  "max": 4800,   "M50": 4500},
        {"name": "HMA 6000-10000",  "min": 7500,  "max": 8500,   "M50": 8000},
        {"name": "HMA 10000-15000", "min": 12000, "max": 13000,  "M50": 12500}
    ]

    @staticmethod
    def get_grading(required_M50_kg):
        """
        Selects the lightest standard grading where Class_M50 >= Required_M50.
        This ensures the selected rock class is sufficient for stability.
        
        Args:
            required_M50_kg (float): Theoretical required mass in kg.
            
        Returns:
            dict: The selected grading dictionary or None if no suitable class exists.
        """
        sorted_grades = sorted(RockStandards.GRADINGS, key=lambda x: x.get('M50', float('inf')))
        for grade in sorted_grades:
            if grade.get('M50', 0) >= required_M50_kg:
                return grade
        return None 

    @staticmethod
    def get_closest_grading(target_M50_kg):
        """
        Selects the grading whose M50 is closest to the target M50.
        Typically used for sizing armourstone layers.
        
        Args:
            target_M50_kg (float): Target mass in kg.
            
        Returns:
            dict: The closest grading dictionary.
        """
        sorted_grades = sorted(RockStandards.GRADINGS, key=lambda x: x.get('M50', float('inf')))
        closest_grade = None
        min_diff = float('inf')
        for grade in sorted_grades:
            val = grade.get('M50')
            if val is None: continue
            diff = abs(val - target_M50_kg)
            if diff < min_diff:
                min_diff = diff
                closest_grade = grade
        return closest_grade

class CoastalConstants:
    """
    Physical constants used throughout the calculations.
    """
    G = 9.80665  # Standard gravity (m/s^2)

# ==============================================================================
# 3. CORE LOGIC CLASSES
# ==============================================================================

class DesignParameters:
    """
    Container for all design input parameters.
    Calculates derived geometry and material properties (e.g., density ratio Delta, permeability Cp).
    
    References:
    - Etemad-Shahidi et al. (2020) for the physical permeability coefficient (Cp).
    """
    def __init__(self, Hs, Tm10, h_toe, slope_m, foreshore_m, P, D_ratio, Sd, N, rho_r, rho_w):
        # Hydraulic Inputs
        self.Hs = Hs
        self.Tm10 = Tm10
        self.h_toe = h_toe
        self.N = N
        self.Sd = Sd
        
        # Geometry Inputs
        self.slope_m = slope_m
        self.foreshore_m = foreshore_m
        
        # Material Inputs
        self.rho_r = rho_r
        self.rho_w = rho_w
        self.P = P  
        self.D_ratio = D_ratio 
        
        # Derived Geometry (Slope Angle)
        self.cot_alpha = slope_m
        self.alpha_rad = math.atan(1.0 / slope_m)
        self.alpha_deg = math.degrees(self.alpha_rad)
        
        # Derived Material Properties (Relative Density)
        self.Delta = (rho_r / rho_w) - 1.0
        
        # Derived Physical Permeability Coefficient Cp
        # Cp relates the core-to-armor size ratio to hydraulic stability.
        self.Cp = (1 + (self.D_ratio**0.3))**0.6

class WaveMechanics:
    """
    Handles hydraulic wave transformation and parameter calculation.
    
    Key Functions:
    - solve_wavelength: Iterative solution for linear dispersion relation.
    - analyze: Comprehensive calculation of wave kinematics and breaking characteristics.
    
    References:
    - Coastal Engineering Manual (2002).
    - Goda (2010) for wave parameter definitions.
    """
    @staticmethod
    def solve_wavelength(T, h):
        """
        Calculates Wavelength L using Newton-Raphson iteration for the dispersion relation.
        
        L = (g * T^2 / 2*pi) * tanh(2*pi*h / L)
        
        Args:
            T (float): Wave period (s).
            h (float): Water depth (m).
            
        Returns:
            float: Wavelength L (m).
        """
        pi = math.pi
        L0 = CoastalConstants.G * (T**2) / (2 * pi) # Deep water wavelength
        if h <= 0: return 0.0

        k0h = 2 * pi * h / L0
        # Initial Guess approximation (Carvalho 2006)
        term = (math.pow(6.0/5.0, k0h)) * math.sqrt(k0h)
        L = L0 * math.tanh(term)
        
        # Newton-Raphson Iteration
        dL = 1.0
        delta = 0.00000001
        iter_count = 0
        
        while abs(dL / L) >= delta and iter_count < 100:
            val = 2 * pi * h / L
            f1 = L - L0 * math.tanh(val)
            val_delta = 2 * pi * h / (L + delta)
            f2 = (L + delta) - L0 * math.tanh(val_delta)
            denom = f2 - f1
            if denom == 0: break
            dL = delta * f1 / denom
            L = L - dL
            iter_count += 1
        return L

    @staticmethod
    def analyze(params):
        """
        Performs full hydraulic analysis at the structure toe.
        Computes L, C, Cg, Deep Water Steepness, Local Steepness, and Surf Similarity.
        Determines the breaker type (Surging, Plunging, etc.).
        """
        g = CoastalConstants.G
        
        # Basic Wave Parameters
        L0 = (g * params.Tm10**2) / (2 * math.pi)
        L_toe = WaveMechanics.solve_wavelength(params.Tm10, params.h_toe)
        
        # Celerity
        if params.Tm10 > 0:
            C = L_toe / params.Tm10
        else:
            C = 0
            
        # Group Celerity (Cg)
        k = 2 * math.pi / L_toe if L_toe > 0 else 0
        kh = k * params.h_toe
        if kh > 20: n = 0.5
        elif kh <= 0: n = 1.0
        else: n = 0.5 * (1 + (2 * kh) / math.sinh(2 * kh))
        Cg = n * C

        # Wave Steepness
        s_m10 = params.Hs / L0 # Deep water steepness
        s_local = params.Hs / L_toe if L_toe > 0 else 0 # Local steepness at toe

        # Surf Similarity Parameter (Iribarren Number) based on deep water steepness
        if s_m10 > 0:
            xi_m10 = math.tan(params.alpha_rad) / math.sqrt(s_m10)
        else:
            xi_m10 = 0
            
        # Relative Depth
        if params.Hs > 0:
            rel_depth = params.h_toe / params.Hs
        else:
            rel_depth = 999.0
        
        # Breaker Type Classification
        if xi_m10 < 0.5: breaker_type = "Spilling"
        elif xi_m10 < 1.8: breaker_type = "Plunging"
        elif xi_m10 < 3.0: breaker_type = "Surging"
        else: breaker_type = "Collapsing/Surging"

        return {
            'L0': L0, 'L_toe': L_toe, 'C': C, 'Cg': Cg,
            's_m10': s_m10, 's_local': s_local,
            'xi_m10': xi_m10, 'rel_depth': rel_depth,
            'breaker_type': breaker_type
        }

class StabilityCalculator:
    """
    Implements various empirical stability formulas from the scientific literature.
    Each method calculates the required nominal diameter (Dn50) and Stability Number (Ns).
    """
    @staticmethod
    def get_mass(Dn50, rho_r):
        """Calculates mass from nominal diameter and density."""
        return rho_r * (Dn50**3)

    @staticmethod
    def calc_hudson_kd(Hs, Delta, Dn50, cot_alpha):
        """Back-calculates the Hudson stability coefficient Kd."""
        if Dn50 <= 0: return 0.0
        Ns = (1.27 * Hs) / (Delta * Dn50)
        return (Ns**3) / cot_alpha

    @staticmethod
    def hudson(p, mech):
        """
        Hudson (1959) Formula.
        Source: Shore Protection Manual (1984) / Coastal Engineering Manual (2002).
        Legacy formula used primarily as a baseline comparison.
        """
        # Hudson Kd based on Zone Classification depending on Relative Depth
        rel_depth = mech['rel_depth']
        if rel_depth > 3.0:
            Kd = 4.0 
        elif 1.5 < rel_depth <= 3.0:
            Kd = 3.5 
        elif 0.5 < rel_depth <= 1.5:
            Kd = 3.0
        else:
            Kd = 2.0

        Dn50 = (1.27 * p.Hs) / (p.Delta * (Kd * p.cot_alpha)**(1.0/3.0))
        Ns = (1.27 * p.Hs) / (p.Delta * Dn50)
        return Dn50, Ns, "Hudson (1959) - Legacy Ref"

    @staticmethod
    def vdm_2021(p, mech):
        """
        Van der Meer (2021) Formula.
        Source: Van der Meer, J. W. (2021). "Rock armour slope stability under wave attack..."
        Rewritten formula utilizing Spectral Period (Tm-1,0) for deep/intermediate water.
        Differentiates between Plunging and Surging waves.
        """
        c_pl = 6.49; c_su = 0.97
        
        # Calculate critical surf similarity parameter for transition
        term_crit = (c_pl / c_su) * (p.P**0.31) * math.sqrt(math.tan(p.alpha_rad))
        xi_cr = term_crit ** (1.0 / (p.P + 0.5))
        
        damage_term = (p.Sd / math.sqrt(p.N))**0.2
        
        if mech['xi_m10'] < xi_cr:
            # Plunging Waves
            Ns = c_pl * (p.P**0.18) * damage_term * (mech['xi_m10']**-0.5)
            note = f"Van der Meer 2021 (Plunging: xi < {xi_cr:.2f})"
        else:
            # Surging Waves
            Ns = c_su * (p.P**-0.13) * damage_term * math.sqrt(p.cot_alpha) * (mech['xi_m10']**p.P)
            note = f"Van der Meer 2021 (Surging: xi > {xi_cr:.2f})"
            
        Dn50 = p.Hs / (p.Delta * Ns)
        return Dn50, Ns, note

    @staticmethod
    def van_gent_2003_modified(p, mech):
        """
        Van Gent Modified (2003) Formula.
        Source: Van Gent et al. (2003). "Stability of rock slopes with shallow foreshores."
        Incorporates the influence of shallow water wave height distribution (H2%/Hs).
        """
        # Estimate H2%/Hs ratio based on relative depth (Van der Meer et al., 2024 analysis of Van Gent data)
        # Deep water: ratio approx 1.4
        # Shallow transition (h/Hm0 approx 1.5): ratio dips to approx 1.2
        # Very shallow: ratio recovers
        if mech['rel_depth'] >= 3.0:
            h2_ratio = 1.4
        elif mech['rel_depth'] < 1.5:
            # Rebounding in very shallow water, linear approx back to 1.4 at h=0
            h2_ratio = 1.2 + 0.2 * (1.5 - mech['rel_depth']) / 1.5
            h2_ratio = min(h2_ratio, 1.4)
        else:
            # Linear interpolation between 3.0 (1.4) and 1.5 (1.2)
            h2_ratio = 1.2 + 0.2 * (mech['rel_depth'] - 1.5) / 1.5

        c_pl = 8.4
        c_su = 1.3
        
        # Critical surf similarity
        # formula: xi_c = ( (c_pl/c_su) * P^0.31 * sqrt(tan_alpha) ) ^ (1/(P+0.5))
        # Note: Standard VdM uses 6.49/0.97. Van Gent mod uses 8.4/1.3.
        term_crit = (c_pl / c_su) * (p.P**0.31) * math.sqrt(math.tan(p.alpha_rad))
        xi_cr = term_crit ** (1.0 / (p.P + 0.5))

        damage_term = (p.Sd / math.sqrt(p.N))**0.2
        ratio_term = (h2_ratio)**(-1)

        if mech['xi_m10'] < xi_cr:
            # Plunging
            Ns = 8.4 * (p.P**0.18) * damage_term * (mech['xi_m10']**-0.5) * ratio_term
            note = f"Van Gent Modified (2003) Plunging (H2%/Hs={h2_ratio:.2f})"
        else:
            # Surging
            # Exponent for xi is P (same as VdM 1988/2021 structure)
            Ns = 1.3 * (p.P**-0.13) * damage_term * math.sqrt(p.cot_alpha) * (mech['xi_m10']**p.P) * ratio_term
            note = f"Van Gent Modified (2003) Surging (H2%/Hs={h2_ratio:.2f})"

        Dn50 = p.Hs / (p.Delta * Ns)
        return Dn50, Ns, note

    @staticmethod
    def van_gent_2003_simplified(p, mech):
        """
        Van Gent et al. (2003) Simplified Formula.
        Source: Van Gent et al. (2003). "Stability of rock slopes with shallow foreshores."
        Calibrated for shallow foreshores where wave height distribution is truncated.
        """
        c_VG = 1.75
        damage_term = (p.Sd / math.sqrt(p.N))**0.2
        perm_term = 1.0 + p.D_ratio
        
        Ns = c_VG * math.sqrt(p.cot_alpha) * perm_term * damage_term
        Dn50 = p.Hs / (p.Delta * Ns)
        return Dn50, Ns, "Van Gent Simplified (2003)"

    @staticmethod
    def eldrup_andersen_2019(p, mech):
        """
        Eldrup & Andersen (2019) Formula.
        Source: Eldrup & Andersen (2019). "Extension of shallow water to nonlinear waves."
        Address stability in shallow water focusing on highly nonlinear waves.
        """
        c_EA1 = 4.5; c_EA2 = 3.1
        damage_term = (p.Sd / math.sqrt(p.N))**0.2
        
        # Plunging condition
        Ns_pl = c_EA1 * damage_term * (1.6**p.P) * (mech['xi_m10']**(0.4*p.P - 0.67))
        
        # Surging condition
        cot_term = min(p.cot_alpha, 2.0)**0.23
        Ns_su = c_EA2 * damage_term * (p.P**0.17) * cot_term
        
        if mech['xi_m10'] < 2.5: 
            Ns = Ns_pl; note = "Eldrup (Plunging)"
        else:
            Ns = Ns_su; note = "Eldrup (Surging)"
            
        Dn50 = p.Hs / (p.Delta * Ns)
        return Dn50, Ns, note

    @staticmethod
    def etemad_shahidi_2020(p, mech):
        """
        Etemad-Shahidi et al. (2020) Formula.
        Source: Etemad-Shahidi et al. (2020). "On the stability of rock armored rubble mound structures."
        Unified formula for deep/shallow water using physical permeability parameter (Cp).
        Includes foreshore slope correction for depth-limited conditions.
        """
        c_ES1 = 4.5; c_ES2 = 3.9
        m = 1.0 / p.foreshore_m
        term_Nw = p.N**(-1.0/10.0)
        term_Sd = p.Sd**(1.0/6.0)
        
        # Foreshore correction
        if mech['rel_depth'] < 3.0:
            foreshore_factor = max(0.1, (1 - 3*m))
            note_suffix = " (Depth-Limited)"
        else:
            foreshore_factor = 1.0
            note_suffix = " (Deep)"

        if mech['xi_m10'] < 1.8:
            # Plunging
            Ns = c_ES1 * p.Cp * term_Nw * term_Sd * (mech['xi_m10']**(-7.0/12.0)) * foreshore_factor
            note = "Etemad-Shahidi (Plunging)" + note_suffix
        else:
            # Surging
            Ns = c_ES2 * p.Cp * term_Nw * term_Sd * (mech['xi_m10']**(-1.0/3.0)) * foreshore_factor
            note = "Etemad-Shahidi (Surging)" + note_suffix
            
        Dn50 = p.Hs / (p.Delta * Ns)
        return Dn50, Ns, note

    @staticmethod
    def scaravaglione_mod_vg(p, mech):
        """
        Scaravaglione (Modified VG 2025) Formula.
        Source: Scaravaglione et al. (2025). "The influence of shallow water on rock armour stability."
        Modified Van Gent for extremely shallow water / swash zone.
        Recalibrated coefficient c_VG = 3.3.
        """
        c_VG_new = 3.3
        perm_term = 1.0 + p.D_ratio
        s_safe = max(mech['s_m10'], 0.005)
        steep_term = s_safe**0.1
        damage_term = (p.Sd / math.sqrt(p.N))**0.2 
        
        Ns = c_VG_new * math.sqrt(p.cot_alpha) * perm_term * steep_term * damage_term
        
        # Apply safety cap on Kd to prevent unrealistic stability in extreme shallows
        Kd_calc = (Ns**3) / p.cot_alpha
        note = "Scaravaglione (Mod. VG)"
        if Kd_calc > 5.0:
            Ns = (5.0 * p.cot_alpha)**(1.0/3.0)
            note += " [Capped Kd=5.0]"
            
        Dn50 = p.Hs / (p.Delta * Ns)
        return Dn50, Ns, note

    @staticmethod
    def scaravaglione_mod_es(p, mech):
        """
        Scaravaglione (Modified ES 2025) Formula.
        Source: Scaravaglione et al. (2025). "The influence of shallow water on rock armour stability."
        Modified Etemad-Shahidi for very shallow water to address over-prediction.
        Decouples wave steepness term.
        """
        if mech['xi_m10'] < 1.8:
            return None, None, "N/A (Req. xi > 1.8)"
            
        c_ES_new = 3.55
        term_N = p.N**(-0.1)
        term_Sd = p.Sd**(1.0/6.0)
        term_cot = p.cot_alpha**(1.0/3.0)
        s_safe = max(mech['s_m10'], 0.005)
        term_s = s_safe**(1.0/20.0) 
        
        Ns = c_ES_new * p.Cp * term_N * term_cot * term_Sd * term_s
        
        # Apply safety cap on Kd
        Kd_calc = (Ns**3) / p.cot_alpha
        note = "Scaravaglione (Mod. ES)"
        if Kd_calc > 5.0:
            Ns = (5.0 * p.cot_alpha)**(1.0/3.0)
            note += " [Capped Kd=5.0]"
            
        Dn50 = p.Hs / (p.Delta * Ns)
        return Dn50, Ns, note

# ==============================================================================
# 4. INTELLIGENT SELECTOR & REPORTING
# ==============================================================================

class IntelligentDesignEngine:
    """
    The core decision-making engine. It compares results from all stability calculators,
    analyzes the hydraulic context (water depth, wave steepness), and recommends the
    most scientifically appropriate formula based on the detected zone.
    """
    def __init__(self, output_file):
        self.output_file = output_file
        self.log_buffer = []

    def log(self, message=""):
        """Appends output to buffer and prints to console."""
        print(message)
        self.log_buffer.append(message)

    def save_logs(self):
        """Writes the buffered log to the output file."""
        try:
            with open(self.output_file, "w", encoding="utf-8") as f:
                f.write("\n".join(self.log_buffer))
            print(f"\n[System] Report saved to {os.path.abspath(self.output_file)}")
        except Exception as e:
            print(f"Error saving file: {e}")

    def log_inputs(self, p):
        """Logs the initial user input parameters."""
        self.log("ROCK SLOPE STABILITY CALCULATOR")
        self.log("\n" + "="*95)
        self.log("   1. DESIGN INPUT PARAMETERS")
        self.log("="*95)
        self.log(f"{'PARAMETER':<35} | {'VALUE':<10} | {'UNIT'}")
        self.log("-" * 95)
        self.log(f"{'Significant Wave Height (Hs)':<35} | {p.Hs:<10.2f} | m")
        self.log(f"{'Spectral Period (Tm-1,0)':<35} | {p.Tm10:<10.2f} | s")
        self.log(f"{'Water Depth (h_toe)':<35} | {p.h_toe:<10.2f} | m")
        self.log(f"{'Structure Slope':<35} | 1:{p.slope_m:<8.2f} | (V:H)")
        self.log(f"{'Foreshore Slope':<35} | 1:{p.foreshore_m:<8.2f} | (V:H)")
        self.log(f"{'Permeability (P_notional)':<35} | {p.P:<10.2f} | (-)")
        self.log(f"{'Physical Permeability (Cp)':<35} | {p.Cp:<10.2f} | (-)")
        self.log(f"{'Damage Level (S)':<35} | {p.Sd:<10.1f} | (-)")
        self.log("-" * 95)

    def log_hydraulics(self, mech):
        """Logs the calculated hydraulic parameters and identifies the hydraulic zone."""
        self.log("\n" + "="*95)
        self.log("   2. CALCULATED HYDRAULIC PARAMETERS")
        self.log("="*95)
        self.log(f"{'PARAMETER':<40} | {'VALUE':<10} | {'UNIT'}")
        self.log("-" * 95)
        self.log(f"{'Deep Water Wavelength (L0)':<40} | {mech['L0']:<10.2f} | m")
        self.log(f"{'Wavelength at Toe (L_toe)':<40} | {mech['L_toe']:<10.2f} | m")
        self.log(f"{'Wave Celerity at Toe (C)':<40} | {mech['C']:<10.2f} | m/s")
        self.log(f"{'Group Celerity at Toe (Cg)':<40} | {mech['Cg']:<10.2f} | m/s")
        self.log(f"{'Deep Water Steepness (s_m-1,0)':<40} | {mech['s_m10']:<10.5f} | (-)")
        self.log(f"{'Local Wave Steepness (s_local)':<40} | {mech['s_local']:<10.5f} | (-)")
        self.log(f"{'Surf Similarity (xi_m-1,0)':<40} | {mech['xi_m10']:<10.2f} | (-)")
        self.log(f"{'Breaker Type (Visual/Physical)':<40} | {mech['breaker_type']:<10} | (-)")
        self.log(f"{'Relative Depth (h/Hm0)':<40} | {mech['rel_depth']:<10.2f} | (-)")
        self.log(f"{'Relative Depth (h/L0)':<40} | {(mech['rel_depth']*mech['s_m10']):<10.3f} | (-)")
        
        # Zone Classification Output logic based on Relative Depth
        rel_depth = mech['rel_depth']
        if rel_depth > 3.0:
            zone = "ZONE 1: Deep to Intermediate (h/Hm0 > 3.0)"
        elif 1.5 < rel_depth <= 3.0:
            zone = "ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0)"
        elif 0.5 < rel_depth <= 1.5:
            zone = "ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5)"
        else:
            zone = "ZONE 4: Extremely Shallow Water (h/Hm0 <= 0.5)"
            
        self.log(f"{'Hydraulic Zone':<40} | {zone:<40}")
        self.log("-" * 95)

    def evaluate_and_recommend(self, params, mech, results):
        """
        Analyzes hydraulic context to select the best formula and provides detailed justification.
        References specific scientific conclusions from the bibliography.
        """
        self.log("\n" + "="*95)
        self.log("   3. FORMULA SELECTION & JUSTIFICATION")
        self.log("="*95)

        rel_depth = mech['rel_depth'] # h / Hm0
        justification = []
        rec_res = None
        
        # --- Helper to find results safely ---
        def get_res(name_sub):
            for r in results:
                if name_sub in r['name']: return r
            return None

        res_vdm = get_res("Van der Meer (2021)")
        res_vg_mod = get_res("Van Gent Modified (2003)")
        res_vg_simp = get_res("Van Gent Simplified")
        res_es = get_res("Etemad-Shahidi")
        res_mod_vg = get_res("Mod. VG")
        res_mod_es = get_res("Mod. ES")
        res_eldrup = get_res("Eldrup")

        # --- EXTENSIVE JUSTIFICATION LOGIC ---
        
        # ZONE 1: Deep to Intermediate (h/Hm0 > 3.0)
        if rel_depth > 3.0:
            regime = "Deep / Intermediate Water"
            rec_res = res_vdm
            justification.append(f"### 1. HYDRAULIC CONTEXT: {regime}")
            justification.append(f"   The structure is located in deep water relative to the wave height (h/Hm0 = {rel_depth:.2f} > 3.0).")
            justification.append("   In this regime, the wave height distribution strictly follows the **Rayleigh distribution**.")
            justification.append("   Key characteristics:")
            justification.append("     * The ratio H2%/Hs is constant at approximately 1.4.")
            justification.append("     * Wave breaking is limited to whitecapping or direct interaction with the armor layer.")
            justification.append("     * The spectral shape is standard (JONSWAP/Pierson-Moskowitz), and energy transfer to low-frequencies is minimal.")
            justification.append("")
            justification.append("### 2. FORMULA COMPARISON & ANALYSIS")
            
            # Van der Meer Analysis
            justification.append("   **A. Van der Meer (2021 Rewritten) [RECOMMENDED]**")
            justification.append("      * **Advantages:** This formula is the modernized industry standard.")
            justification.append("        Van der Meer (2021) rewrote the original formula to use the spectral period (Tm-1,0),")
            justification.append("        eliminating the influence of spectral shape.")
            justification.append("        Van der Meer et al. (2024) confirmed its validity for h/Hm0 > 1.5, preferring Hm0 over H1/3 for nonlinear waves.")
            justification.append("      * **Physics:** It correctly assumes a Rayleigh distribution of wave heights, aligning with the actual deep-water statistics.")
            
            # Van Gent Analysis
            justification.append("   **B. Van Gent Modified (2003)**")
            justification.append("      * **Context:** This formula incorporates the ratio H2%/Hs.")
            justification.append("        In deep water, with H2%/Hs = 1.4, this formula essentially converges closely with the Van der Meer predictions.")
            justification.append("        However, its specific calibration was focused on the effects of shallow foreshores.")
            
            # Etemad-Shahidi Analysis
            justification.append("   **C. Etemad-Shahidi (2020)**")
            justification.append("      * **Comparison:** Etemad-Shahidi (2020) provides a robust formula validated for both deep and shallow water.")
            justification.append("        It introduces a physical permeability parameter (D_core/D_armor) to replace the nominal P factor,")
            justification.append("        reducing uncertainty. However, Van der Meer remains the primary standard for deep water.")
            
            # Check for divergence
            if res_es and res_vdm and res_vdm['Dn'] and res_es['Dn']:
                diff_pct = abs(res_vdm['Dn'] - res_es['Dn']) / res_vdm['Dn']
                if diff_pct > 0.10:
                    justification.append("      * **Note on Divergence:** The result deviates from Van der Meer here. This typically occurs in the 'transition zone'")
                    justification.append("        of the surf similarity parameter (xi approx 2.0 - 4.5). Etemad-Shahidi transitions to Surging physics")
                    justification.append("        earlier (xi > 1.8), predicting higher stability, whereas Van der Meer maintains Plunging physics")
                    justification.append("        (lower stability) until a higher critical threshold. Van der Meer is more conservative here.")
                else:
                    justification.append("        It typically converges with Van der Meer here.")
            
            justification.append("")
            justification.append("### 3. FINAL JUSTIFICATION")
            justification.append("   **Use [Van der Meer (2021 Rewritten)]**.")
            justification.append("   It provides the most theoretically consistent result for non-depth-limited waves.")
            
            if res_es and res_vdm and res_vdm['Dn'] and res_es['Dn']:
                diff = abs(res_vdm['Dn'] - res_es['Dn'])
                justification.append(f"   *Verification:* Etemad-Shahidi yields Dn50 = {res_es['Dn']:.3f}m (Difference: {diff:.3f}m).")

        # ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0)
        elif 1.5 < rel_depth <= 3.0:
            regime = "Shallow Water (Transition Zone)"
            rec_res = res_vdm
            justification.append(f"### 1. HYDRAULIC CONTEXT: {regime}")
            justification.append(f"   The structure is in the transition zone (1.5 < h/Hm0 = {rel_depth:.2f} <= 3.0).")
            justification.append("   Key characteristics:")
            justification.append("     * **Spectral Truncation:** The largest waves in the spectrum break on the foreshore.")
            justification.append("     * **Distribution Shift:** The wave height distribution deviates from Rayleigh; H2%/Hm0 drops below 1.4.")
            justification.append("     * **Shoaling:** Significant shoaling modifies the wave shape before impact, creating peaked crests and flat troughs.")
            justification.append("")
            justification.append("### 2. FORMULA COMPARISON & ANALYSIS")
            
            # Van der Meer Analysis
            justification.append("   **A. Van der Meer (2021)**")
            justification.append("      * **Advantages:** Van der Meer et al. (2024) extensively re-analyzed shallow water data and concluded")
            justification.append("        that the rewritten Van der Meer formula (using Tm-1,0) is valid down to h/Hm0 = 1.5.")
            justification.append("        It performs reasonably well, with slightly less reliability in the 1.0 < h/Hm0 < 1.5 range.")
            justification.append("      * **Note:** For nonlinear waves in this zone, using Hm0 is preferred over H1/3 to avoid deviations.")
            
            # Van Gent Modified Analysis
            justification.append("   **B. Van Gent Modified (2003)**")
            justification.append("      * **Constraint:** This formula explicitly relies on the ratio H2%/Hs. Research by Van der Meer et al. (2024)")
            justification.append("        highlights that predicting H2% accurately in this transition zone (where the ratio dips to ~1.2)")
            justification.append("        is notoriously inaccurate without physical modeling. The formula is valid, but the input uncertainty is high.")
            
            # Van Gent Simplified Analysis
            justification.append("   **C. Van Gent et al. (2003) Simplified**")
            justification.append("      * **Context:** This formula was specifically derived for shallow foreshores.")
            justification.append("        However, Van der Meer et al. (2024) found that the simplified formula often does not match the data")
            justification.append("        in the surging domain as well as the rewritten Van der Meer formula.")
            
            justification.append("")
            justification.append("### 3. FINAL JUSTIFICATION")
            justification.append("   **Use [Van der Meer (2021 Rewritten)]**.")
            justification.append("   Recent research (2024) confirms its validity in this depth range (h/Hm0 > 1.5), favoring it over simplified methods")
            justification.append("   due to the uncertainties in predicting H2% required for the Van Gent Modified formula.")

        # ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5)
        elif 0.5 < rel_depth <= 1.5:
            regime = "Very Shallow Water (Surf Zone)"
            # Primary: Mod ES (Scaravaglione), Fallback: Eldrup or VG
            if res_mod_es and res_mod_es['Dn'] is not None:
                rec_res = res_mod_es
                justification.append(f"### 1. HYDRAULIC CONTEXT: {regime}")
                justification.append(f"   The structure is in the surf zone (0.5 < h/Hm0 = {rel_depth:.2f} <= 1.5).")
                justification.append("   Key characteristics:")
                justification.append("     * **Severe Breaking:** Waves are constantly breaking.")
                justification.append("     * **Saturation:** Wave height is depth-limited (H ~ 0.5h). Increasing offshore energy does not increase load.")
                justification.append("     * **Infragravity Dominance:** Scaravaglione et al. (2025) and VdM (2024) note that infragravity waves")
                justification.append("       begin to dominate the spectrum, causing Tm-1,0 to increase massively (up to 4x).")
                justification.append("     * **Formula Deviation:** Standard formulas fail here because the stability curves flatten out (Horizontal Trend).")
                justification.append("")
                justification.append("### 2. FORMULA COMPARISON & ANALYSIS")
                
                # Standard Formulas Analysis
                justification.append("   **A. Standard Formulas (VdM, VG, Standard ES)**")
                justification.append("      * **Failure Mode:** Scaravaglione et al. (2025) demonstrated that these formulas fail to converge here.")
                justification.append("        Using the inflated spectral period results in over-predicted stability (for surging) or under-predicted (for plunging).")
                
                # Modified ES Analysis
                justification.append("   **B. Scaravaglione (Modified ES 2025) [RECOMMENDED]**")
                justification.append("      * **Advantages:** This formula explicitly **decouples** the wave steepness term (s^0.05) from the structure slope term.")
                justification.append("      * **Physics:** It is calibrated specifically for surging/bore-like waves in the surf zone using new coefficients (c_ES,new=3.55).")
                justification.append("      * **Safety Cap Applied:** The formula has a weak dependence on steepness. For very long swells, it may predict")
                justification.append("        theoretically high stability (Kd > 10). This system has capped Kd at 5.0 to ensure physical realism.")
                
                justification.append("")
                justification.append("### 3. FINAL JUSTIFICATION")
                justification.append("   **Use [Scaravaglione (Modified ES 2025)]**.")
                justification.append("   This represents the state-of-the-art for broken waves in very shallow water, correcting the overestimation of damage.")
            else:
                rec_res = res_vg_simp
                justification.append(f"### 1. HYDRAULIC CONTEXT: {regime} (Plunging)")
                justification.append("   The structure is in the surf zone, but conditions are **Plunging** (xi < 1.8).")
                justification.append("   The Modified ES formula is only calibrated for surging/bore conditions.")
                justification.append("")
                justification.append("### 3. FINAL JUSTIFICATION")
                justification.append("   **Use [Van Gent Simplified (2003)]**.")
                justification.append("   It acts as a robust fallback. Caution is advised as damage may be underpredicted for impermeable structures.")

        # ZONE 4: Extremely Shallow Water (h/Hm0 <= 0.5)
        else:
            regime = "Extremely Shallow Water (Swash Zone)"
            rec_res = res_mod_vg
            justification.append(f"### 1. HYDRAULIC CONTEXT: {regime}")
            justification.append(f"   The structure is located effectively in the **Swash Zone** (h/Hm0 = {rel_depth:.2f} <= 0.5).")
            justification.append("   Key characteristics:")
            justification.append("     * **Aeration:** High air entrainment reduces the effective fluid density and buoyancy of the rocks.")
            justification.append("     * **Impact:** Wave impact is characterized by a high-velocity turbulent bore.")
            justification.append("     * **Hydrostatics:** The hydrostatic cushioning effect is negligible.")
            justification.append("")
            justification.append("### 2. FORMULA COMPARISON & ANALYSIS")
            
            # Standard VG Analysis
            justification.append("   **A. Standard Van Gent (2003)**")
            justification.append("      * **Disadvantages:** Using the standard coefficient (1.75) here is **Unsafe**.")
            justification.append("        Scaravaglione et al. (2025) showed that stability is significantly lower than predicted by intermediate-depth formulas")
            justification.append("        due to the lack of buoyancy and intense turbulence.")
            
            # Modified VG Analysis
            justification.append("   **B. Scaravaglione (Modified VG 2025) [RECOMMENDED]**")
            justification.append("      * **Advantages:** This formula uses a recalibrated coefficient (C_VG = 3.3 instead of 1.75).")
            justification.append("      * **Physics:** It explicitly accounts for the increased instability in the swash zone, correcting the underestimation")
            justification.append("        of damage by the original VG formula in this specific regime.")
            
            justification.append("")
            justification.append("### 3. FINAL JUSTIFICATION")
            justification.append("   **Use [Scaravaglione (Modified VG 2025)]**.")
            justification.append("   It provides the necessary safety margin for swash zone instability where standard formulas fail.")

        # --- COMPARISON TABLE ---
        self.log(f"COMPARISON OF RESULTS:")
        self.log(f"{'Method':<30} | {'Ns (-)':<8} | {'Dn50 (m)':<10} | {'M50 (kg)':<10} | {'Kd_eq (-)':<8} | {'NOTES'}")
        self.log("-" * 115) 
        
        for res in results:
            if res['Dn'] is None: continue
            mass = StabilityCalculator.get_mass(res['Dn'], params.rho_r)
            # Removed the marker as requested
            self.log(f"{res['name']:<30} | {res['Ns']:<8.4f} | {res['Dn']:<10.3f} | {mass:<10.0f} | {res['Kd']:<8.2f} | {res['note']}")

        self.log("\nJUSTIFICATION & ANALYSIS:")
        for line in justification:
            self.log(line)
        self.log("-" * 95)
        
        return rec_res

    def present_layer_design(self, params, target_dn, target_mass, grading_EN13383=True):
        """
        Generates the design output for Primary Armor and Underlayers.
        Calculates weights, diameters, and packing densities.
        Uses EN 13383 standard gradings if enabled, or a custom power law distribution.
        """
        
        # --- 1. PRIMARY ARMOR LAYER ---
        target_weight_kn = target_mass * CoastalConstants.G / 1000.0
        
        self.log("\n" + "="*95)
        if grading_EN13383:
            self.log("   4. ROCK ARMOUR LAYER DESIGN (EN 13383 Standard)")
        else:
            self.log("   4. ROCK ARMOUR LAYER DESIGN (Custom Grading)")
        self.log("="*95)
        
        self.log("PRIMARY ARMOR LAYER")
        self.log(f"   Theoretical Required W    : {target_weight_kn:.2f} kN")
        self.log(f"   Theoretical Required M50  : {target_mass:.0f} kg")
        self.log(f"   Theoretical Required Dn50 : {target_dn:.3f} m")
        self.log("-" * 40)

        g = CoastalConstants.G
        gamma_r = params.rho_r * g / 1000.0
        
        # --- LOGIC BRANCH: EN13383 vs Custom ---
        if grading_EN13383:
            # Select appropriate standard grading from database
            grading_armor = RockStandards.get_grading(target_mass)
            if not grading_armor:
                self.log("   [WARNING] No standard EN13383 grading found for this mass.")
                # Return empty values if no grading is found to prevent crash
                w_min_kn = 0; w_max_kn = 0; m_mean_kg = 0; w_mean_kn = 0; actual_dn = 0
                return
            else:
                w_min_kn = grading_armor.get('min', 0) * g / 1000.0
                w_max_kn = grading_armor.get('max', 0) * g / 1000.0
                m_mean_kg = grading_armor.get('M50', 0)
                w_mean_kn = m_mean_kg * g / 1000.0
                actual_dn = (w_mean_kn / gamma_r)**(1.0/3.0)
                grading_name = grading_armor['name']
                
        else:
            # --- CUSTOM POWER LAW CALCULATION ---
            # Used when standard grading is disabled.
            # Using 'x' as Theoretical Required W (in kN)
            x_val = target_weight_kn
            
            # Grading Min Params (Formula Coefficients)
            a_min = 1.056832014477894E+00
            b_min = 1.482769823574055E+00
            c_min = -2.476127406338004E-01
            
            # Formula: (Theoretical Required W) * a/(1+(x/b)**c)
            w_min_kn = target_weight_kn * a_min / (1 + (x_val / b_min)**c_min)
            
            # Grading Max Params (Formula Coefficients)
            a_max = 1.713085676568561E+00
            b_max = 2.460481255856126E+05
            c_max = 1.327263214034671E-01
            
            w_max_kn = target_weight_kn * a_max / (1 + (x_val / b_max)**c_max)
            
            # Design values based on target
            w_mean_kn = target_weight_kn
            m_mean_kg = target_mass
            actual_dn = target_dn
            grading_name = "Custom Grading"

        # Calculate Layer Geometry
        layer_thickness = 2 * 1.0 * actual_dn # n*kt*Dn50 (n=2 layers)
        packing_density_per_m2 = 2 * 1.0 * (1 - 0.30) / (actual_dn**2) # Approx porosity 30%
        packing_density_100m2 = packing_density_per_m2 * 100
        
        # Calculate weights in kg for display
        w_min_kg = (w_min_kn * 1000) / g
        w_max_kg = (w_max_kn * 1000) / g

        self.log(f"   Adopted rock grading                : {grading_name}")
        self.log(f"   Grading Min (Lower Limit)           : {w_min_kn:.2f} kN ({w_min_kg:.0f} kg)")
        self.log(f"   Grading Max (Upper Limit)           : {w_max_kn:.2f} kN ({w_max_kg:.0f} kg)")
        self.log(f"   Representative M50                  : {m_mean_kg:.0f} kg")
        self.log(f"   Nominal Diameter (Dn_rock)          : {actual_dn:.3f} m")
        self.log(f"   Double Layer Thickness              : {layer_thickness:.2f} m")
        self.log(f"   Packing Density [rocks/100m2]       : {packing_density_100m2:.2f}")
        self.log("-" * 95)

        # --- 2. UNDERLAYER ---
        # Standard design practice: Filter layer mass is typically M50/10 to M50/15
        target_mass_ul = m_mean_kg / 10.0
        target_weight_kn_ul = target_mass_ul * CoastalConstants.G / 1000.0
        
        self.log("UNDERLAYER (FILTER LAYER)")
        self.log(f"   Target Weight (M50 / 10)  : {target_weight_kn_ul:.3f} kN")
        self.log(f"   Target Mass (M50 / 10)    : {target_mass_ul:.1f} kg")
        self.log("-" * 40)
        
        if grading_EN13383:
            grading_ul = RockStandards.get_closest_grading(target_mass_ul)
            if grading_ul:
                w_min_kn_ul = grading_ul.get('min', 0) * g / 1000.0
                w_max_kn_ul = grading_ul.get('max', 0) * g / 1000.0
                m_mean_kg_ul = grading_ul.get('M50', 0)
                w_mean_kn_ul = m_mean_kg_ul * g / 1000.0
                actual_dn_ul = (w_mean_kn_ul / gamma_r)**(1.0/3.0)
                grading_name_ul = grading_ul['name']
            else:
                 self.log("   [WARNING] No suitable standard underlayer grading found.")
                 return
        else:
            # Custom Approach for Underlayer (Same Logic as Main Armor)
            x_val_ul = target_weight_kn_ul
            
            # Grading Min Params
            a_min = 1.056832014477894E+00
            b_min = 1.482769823574055E+00
            c_min = -2.476127406338004E-01
            
            w_min_kn_ul = target_weight_kn_ul * a_min / (1 + (x_val_ul / b_min)**c_min)
            
            # Grading Max Params
            a_max = 1.713085676568561E+00
            b_max = 2.460481255856126E+05
            c_max = 1.327263214034671E-01
            
            w_max_kn_ul = target_weight_kn_ul * a_max / (1 + (x_val_ul / b_max)**c_max)
            
            w_mean_kn_ul = target_weight_kn_ul
            m_mean_kg_ul = target_mass_ul
            actual_dn_ul = (w_mean_kn_ul / gamma_r)**(1.0/3.0)
            grading_name_ul = "Custom Grading Underlayer"

        layer_thickness_ul = 2 * 1.0 * actual_dn_ul
        packing_density_per_m2_ul = 2 * 1.0 * (1 - 0.30) / (actual_dn_ul**2)
        packing_density_100m2_ul = packing_density_per_m2_ul * 100
        
        # Calculate weights in kg for display
        w_min_kg_ul = (w_min_kn_ul * 1000) / g
        w_max_kg_ul = (w_max_kn_ul * 1000) / g
        
        self.log(f"   Adopted rock grading                : {grading_name_ul}")
        self.log(f"   Grading Min (Lower Limit)           : {w_min_kn_ul:.2f} kN ({w_min_kg_ul:.0f} kg)")
        self.log(f"   Grading Max (Upper Limit)           : {w_max_kn_ul:.2f} kN ({w_max_kg_ul:.0f} kg)")
        self.log(f"   Representative M50                  : {m_mean_kg_ul:.1f} kg")
        self.log(f"   Nominal Diameter (Dn_rock)          : {actual_dn_ul:.3f} m")
        self.log(f"   Double Layer Thickness              : {layer_thickness_ul:.2f} m")
        self.log(f"   Packing Density [rocks/100m2]       : {packing_density_100m2_ul:.2f}")
            
        self.log("="*95 + "\n")

# ==============================================================================
# 5. MAIN EXECUTION
# ==============================================================================

def get_input(prompt, default_value, data_type=float):
    """
    Utility function to get user input with a default value.
    Supports boolean input conversion (yes/no, true/false).
    """
    user_val = input(f"{prompt} [default: {default_value}]: ").strip()
    if not user_val: return default_value
    try:
        # Handle Boolean explicitly for the new flag
        if data_type == bool:
            return user_val.lower() in ['true', '1', 't', 'yes', 'y']
        return data_type(user_val)
    except: return default_value

def main():
    """
    Main program entry point.
    1. Collects user input for hydraulic and geometric parameters.
    2. Instantiates calculation objects.
    3. Runs hydraulic analysis.
    4. Runs all stability formulas.
    5. Evaluates and selects the best result.
    6. Presents the final design.
    """
    print("\n" + "="*80)
    print("   ROCK SLOPE STABILITY CALCULATOR")
    print("="*80)
    print("Press [Enter] to accept the default value shown in brackets.\n")
    
    # --- Step 1: Input Collection ---
    Hs = get_input("Significant Wave Height Hm0 (at toe) [m]", DEFAULT_CONFIG['Hs'])
    Tm10 = get_input("Spectral Period Tm-1,0 (at toe) [s]", DEFAULT_CONFIG['Tm10'])
    h_toe = get_input("Water Depth at Toe h [m]", DEFAULT_CONFIG['h_toe'])
    slope_m = get_input("Structure Slope (m:1) [cot alpha]", DEFAULT_CONFIG['slope_m']) 
    foreshore_m = get_input("Foreshore Slope (m_f:1) [cot beta]", DEFAULT_CONFIG['foreshore_m'])
    rho_r = get_input("Rock Density [kg/m3]", DEFAULT_CONFIG['rho_r'])
    P = get_input("Notional Permeability P (0.4=Permeable)", DEFAULT_CONFIG['P_notional'])
    D_ratio = get_input("Core/Armor Diameter Ratio (0.3=Typical)", DEFAULT_CONFIG['D_ratio'])
    Sd = get_input("Design Damage Level S (2.0=Start)", DEFAULT_CONFIG['Sd'])
    duration = get_input("Storm Duration [hours]", DEFAULT_CONFIG['duration'])
    
    # Grading Selection Flag
    grading_EN13383 = get_input("Use EN13383 Standard Grading? (True/False)", DEFAULT_CONFIG['grading_EN13383'], bool)

    # Manual Selection Flag
    manual_selection = get_input("Choose stability formula instead of automatic (True/False)", False, bool)
    
    # Calculate Number of Waves
    N_waves = (duration * 3600.0) / Tm10
    
    # --- Step 2: Instantiation ---
    p = DesignParameters(Hs, Tm10, h_toe, slope_m, foreshore_m, P, D_ratio, Sd, N_waves, rho_r, DEFAULT_CONFIG['rho_w'])
    
    calc = StabilityCalculator()
    engine = IntelligentDesignEngine(DEFAULT_CONFIG['output_file'])
    
    engine.log_inputs(p)
    
    # --- Step 3: Hydraulic Analysis ---
    mech = WaveMechanics.analyze(p)
    engine.log_hydraulics(mech)

    # --- Step 4: Stability Calculation (Run all models) ---
    results = []
    
    dn, ns, note = calc.hudson(p, mech)
    kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
    results.append({"name": "Hudson (1959)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})
    
    dn, ns, note = calc.vdm_2021(p, mech)
    kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
    results.append({"name": "Van der Meer (2021)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})

    dn, ns, note = calc.van_gent_2003_modified(p, mech)
    kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
    results.append({"name": "Van Gent Modified (2003)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})

    dn, ns, note = calc.van_gent_2003_simplified(p, mech)
    kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
    results.append({"name": "Van Gent Simplified (2003)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})
    
    dn, ns, note = calc.eldrup_andersen_2019(p, mech)
    kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
    results.append({"name": "Eldrup & Andersen (2019)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})

    dn, ns, note = calc.etemad_shahidi_2020(p, mech)
    kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
    results.append({"name": "Etemad-Shahidi (2020)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})

    dn, ns, note = calc.scaravaglione_mod_vg(p, mech)
    kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
    results.append({"name": "Scaravaglione (Mod. VG 2025)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})

    dn, ns, note = calc.scaravaglione_mod_es(p, mech)
    if dn is not None:
        kd = calc.calc_hudson_kd(p.Hs, p.Delta, dn, p.cot_alpha)
        results.append({"name": "Scaravaglione (Mod. ES 2025)", "Dn": dn, "Ns": ns, "Kd": kd, "note": note})

    # --- Step 5: Evaluation & Recommendation / Manual Override ---
    best_res = None
    
    if manual_selection:
        print("\n" + "="*80)
        print("   MANUAL SELECTION MODE")
        print("="*80)
        
        # Filter only valid results (dn is not None) just in case
        available_formulas = [r for r in results if r.get('Dn') is not None]
        
        for i, res in enumerate(available_formulas):
            print(f"   [{i+1}] {res['name']}")
        
        valid_choice = False
        while not valid_choice:
            try:
                choice_idx = int(input(f"\nEnter the number of your preferred formula (1-{len(available_formulas)}): "))
                if 1 <= choice_idx <= len(available_formulas):
                    best_res = available_formulas[choice_idx - 1]
                    engine.log(f"\n[MANUAL OVERRIDE] User switched selection to: {best_res['name'].upper()}")
                    valid_choice = True
                else:
                    print("Invalid selection. Try again.")
            except ValueError:
                print("Please enter a valid number.")
    
    else:
        # Automatic Evaluation
        best_res = engine.evaluate_and_recommend(p, mech, results)

    # --- Step 6: Final Design Presentation ---
    if best_res:
        m50 = calc.get_mass(best_res['Dn'], p.rho_r)
        engine.present_layer_design(p, best_res['Dn'], m50, grading_EN13383)
        
    engine.save_logs()

if __name__ == "__main__":
    main()