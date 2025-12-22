! ======================================================================================
! PROGRAM DESCRIPTION & METHODOLOGY
! ======================================================================================
!
! 1. PURPOSE:
!    This software calculates the required size and weight of rock armor units for 
!    rubble mound breakwaters and revetments. It performs a comprehensive analysis 
!    of hydraulic stability across different water depth zones (Deep, Shallow, 
!    Very Shallow, and Swash zones), selecting the most appropriate empirical formula 
!    based on the hydraulic regime.
!
! 2. METHODOLOGY & LOGIC:
!    The calculator implements a "Multi-Model Consensus" approach. It computes stability 
!    using several state-of-the-art empirical formulas derived from extensive hydraulic 
!    model testing. It then evaluates the hydraulic context (Relative Depth h/Hm0) 
!    to recommend the scientifically most accurate formula.
!
!    The logic follows a 4-Step Process:
!    a. Input Processing: Parsing wave data, geometry, and material properties.
!    b. Hydraulic Analysis: Computing local wavelength, celerity, wave steepness, 
!       surf similarity (Iribarren number), and determining the breaker type.
!    c. Stability Calculation: Running multiple empirical models (Hudson, Van der Meer, 
!       Van Gent, Etemad-Shahidi, Eldrup & Andersen, Scaravaglione).
!    d. Intelligent Selection: Automatically detecting the "Hydraulic Zone" 
!       (Deep vs. Shallow vs. Swash) to recommend the most valid result.
!
! 3. HYDRAULIC ZONES & FORMULA SELECTION STRATEGY:
!    - ZONE 1: Deep/Intermediate (h/Hm0 > 3.0) -> Recommends Van der Meer (2021).
!    - ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0) -> Recommends Van der Meer (2021) 
!      or Van Gent (2003), depending on spectral shape.
!    - ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5) -> Recommends Scaravaglione (2025) 
!      (Modified ES) to account for wave breaking and horizontal stability trends.
!    - ZONE 4: Extremely Shallow/Swash (h/Hm0 <= 0.5) -> Recommends Scaravaglione (2025) 
!      (Modified VG) to handle infragravity dominance and high turbulence.
!
! 4. BIBLIOGRAPHY & SCIENTIFIC REFERENCES:
!
!    1. U.S. Army Corps of Engineers. (1984). "Shore Protection Manual." Vol. I & II. 
!       Coastal Engineering Research Center, Vicksburg, MS. 
!       Link: https://usace.contentdm.oclc.org/digital/collection/p16021coll11/id/1934/ 
!
!    2. Van der Meer, J. W. (1987). "Stability of breakwater armour layers—Design formulae." 
!       Coastal Engineering, 11(3), 219-239.
!       Link: https://doi.org/10.1016/0378-3839(87)90013-5
!
!    3. Van der Meer, J. W. (1988). "Rock Slopes and Gravel Beaches Under Wave Attack." 
!       Doctoral Thesis, Delft University of Technology.
!       [Foundational work for modern stability formulas]
!       Link: https://repository.tudelft.nl/islandora/object/uuid:404b5nec-hm
!
!    4. Van Gent, M. R. A. (1995). "Wave interaction with permeable coastal structures."
!       Doctoral Thesis, Delft University of Technology.
!       Link: https://repository.tudelft.nl/islandora/object/uuid:7bbff8e4-215d-4bfc-a3af-51cdecb754bd
!
!    5. U.S. Army Corps of Engineers (USACE). (2002). "Coastal Engineering Manual." 
!       Engineer Manual 1110-2-1100, Washington, D.C.
!       Link: https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/
!
!    6. Van Gent, M. R. A., Smale, A. J., & Kuiper, C. (2003). "Stability of rock slopes 
!       with shallow foreshores." Proceedings of Coastal Structures 2003, Portland, OR, 100-112.
!       Link: https://doi.org/10.1061/40733(147)9
!
!    7. Van Gent, M. R. A. (2004). "On the stability of rock slopes." 
!       Environmentally Friendly Coastal Protection: Proceedings of the NATO Advanced 
!       Research Workshop, Varna, Bulgaria.
!       Link: https://doi.org/10.1007/1-4020-3301-X_12
!
!    8. CIRIA, CUR, CETMEF. (2007). "The Rock Manual. The Use of Rock in Hydraulic Engineering." 
!       (2nd edition). C683, CIRIA, London.
!       Link: https://www.ciria.org/ItemDetail?iProductCode=C683
!
!    9. CEN (2013). "EN 13383-1:2013 Armourstone - Part 1: Specification."
!       European Committee for Standardization.
!       Link: https://standards.iteh.ai/catalog/standards/cen/5f6b770f-38ba-4320-ad1a-83d0e29f7db2/en-13383-1-2013
!
!   10. Eldrup, M. R., & Lykke Andersen, T. (2019). "Extension of shallow water rock armour 
!       stability formulae to nonlinear waves." Coastal Engineering, 153, 103536.
!       Link: https://doi.org/10.1016/j.coastaleng.2019.103536
!
!   11. Etemad-Shahidi, A., Bali, M., & Van Gent, M. R. A. (2020). "On the stability of rock 
!       armored rubble mound structures." Coastal Engineering, 158, 103655.
!       Link: https://doi.org/10.1016/j.coastaleng.2020.103655
!
!   12. Van der Meer, J. W. (2021). "Rock armour slope stability under wave attack; the 
!       Van der Meer Formula revisited." Journal of Coastal and Hydraulic Structures, 1.
!       Link: https://doi.org/10.48438/jchs.2021.0008
!
!   13. Van der Meer, J. W., Lykke Andersen, T., & Roge Eldrup, M. (2024). "Rock Armour 
!       Slope Stability under Wave Attack in Shallow Water." Journal of Coastal and 
!       Hydraulic Structures, 4.
!       Link: https://doi.org/10.59490/jchs.2024.0035
!
!   14. Scaravaglione, G., Marino, S., Francone, A., Leone, E., Damiani, L., Tomasicchio, 
!       G. R., Van Gent, M. R. A., & Saponieri, A. (2025). "The influence of shallow water 
!       on rock armour stability." Coastal Engineering, 197, 104657.
!       Link: https://doi.org/10.1016/j.coastaleng.2024.104657
!
! 5. COMPILATION INSTRUCTIONS (MinGW on Windows):
!
!    gfortran -O3 -std=f2008 -o rock_slope_calculator_cli.exe rock_slope_calculator_cli.f90 -static -static-libgfortran -static-libgcc
!
! 6. EXECUTION:
!
!    Interactive Mode:
!      rock_slope_calculator_cli.exe
!
!    Command Line Mode:
!      rock_slope_calculator_cli.exe [Hs] [Tm] [h_toe] [slope] [fore] [rho_r] [P] [D_ratio] [Sd] [Duration] [UseEN13383]
!
!    Example:
!      rock_slope_calculator_cli.exe 2.5 10.0 6.0 2.0 30.0 2650 0.4 0.3 2.0 6.0 true
!
! ======================================================================================

PROGRAM RockSlopeCalculator
    IMPLICIT NONE

    ! ----------------------------------------------------------------------
    ! CONSTANTS
    ! ----------------------------------------------------------------------
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp
    REAL(dp), PARAMETER :: g = 9.80665_dp

    ! ----------------------------------------------------------------------
    ! DATA STRUCTURES
    ! ----------------------------------------------------------------------

    TYPE :: GradingDef
        CHARACTER(LEN=64) :: name
        REAL(dp) :: min_mass ! kg
        REAL(dp) :: max_mass ! kg
        REAL(dp) :: M50      ! kg
    END TYPE GradingDef

    TYPE :: Inputs
        REAL(dp) :: Hs          ! Significant Wave Height (m)
        REAL(dp) :: Tm10        ! Spectral Period Tm-1,0 (s)
        REAL(dp) :: h_toe       ! Water Depth at Toe (m)
        REAL(dp) :: slope_m     ! Structure Slope cot(alpha)
        REAL(dp) :: foreshore_m ! Foreshore Slope cot(beta)
        REAL(dp) :: rho_r       ! Rock Density (kg/m3)
        REAL(dp) :: rho_w       ! Water Density (kg/m3) - Default 1025
        REAL(dp) :: P           ! Notional Permeability
        REAL(dp) :: D_ratio     ! Core/Armor Diameter Ratio
        REAL(dp) :: Sd          ! Design Damage Level
        REAL(dp) :: duration    ! Storm Duration (hours)
        LOGICAL  :: use_en13383 ! Use Standard Grading
        LOGICAL  :: manual_selection ! Interactive manual selection flag
    END TYPE Inputs

    TYPE :: DerivedParams
        REAL(dp) :: cot_alpha
        REAL(dp) :: alpha_rad
        REAL(dp) :: alpha_deg
        REAL(dp) :: Delta       ! Relative buoyant density
        REAL(dp) :: Cp          ! Physical permeability coefficient
        REAL(dp) :: N_waves     ! Number of waves
    END TYPE DerivedParams

    TYPE :: Hydraulics
        REAL(dp) :: L0
        REAL(dp) :: L_toe
        REAL(dp) :: C
        REAL(dp) :: Cg
        REAL(dp) :: s_m10       ! Deep water steepness
        REAL(dp) :: s_local     ! Local steepness
        REAL(dp) :: xi_m10      ! Surf similarity parameter
        REAL(dp) :: rel_depth   ! h / Hm0
        CHARACTER(LEN=32) :: breaker_type
        CHARACTER(LEN=100) :: zone_desc
    END TYPE Hydraulics

    TYPE :: FormulaResult
        CHARACTER(LEN=64) :: name
        REAL(dp) :: Dn50
        REAL(dp) :: Ns
        REAL(dp) :: Kd
        CHARACTER(LEN=100) :: note
        LOGICAL :: valid        ! If formula calculation was successful (Dn > 0)
    END TYPE FormulaResult

    TYPE :: LayerDesign
        CHARACTER(LEN=32) :: layer_name ! "Primary Armor" or "Underlayer"
        CHARACTER(LEN=64) :: grading_name
        REAL(dp) :: target_W_kN
        REAL(dp) :: target_M50_kg
        REAL(dp) :: target_Dn_m
        
        REAL(dp) :: w_min_kn
        REAL(dp) :: w_max_kn
        REAL(dp) :: w_min_kg
        REAL(dp) :: w_max_kg
        
        REAL(dp) :: m_mean_kg
        REAL(dp) :: w_mean_kn
        REAL(dp) :: actual_dn
        REAL(dp) :: thickness
        REAL(dp) :: packing_density ! rocks/100m2
        
        LOGICAL :: design_valid ! If a valid grading/design was found
    END TYPE LayerDesign

    TYPE :: FullReport
        TYPE(Inputs) :: inputs
        TYPE(DerivedParams) :: derived
        TYPE(Hydraulics) :: hydro
        TYPE(FormulaResult), ALLOCATABLE :: comparison(:)
        TYPE(FormulaResult) :: recommended
        CHARACTER(LEN=300), ALLOCATABLE :: justification(:)
        TYPE(LayerDesign) :: armor_layer
        TYPE(LayerDesign) :: underlayer
    END TYPE FullReport

    ! ----------------------------------------------------------------------
    ! VARIABLES
    ! ----------------------------------------------------------------------
    TYPE(GradingDef) :: standard_gradings(16)
    TYPE(Inputs) :: user_inputs
    TYPE(FullReport) :: final_report
    INTEGER :: num_args
    CHARACTER(LEN=100) :: arg_val
    CHARACTER(LEN=100) :: buffer
    INTEGER :: i, iostat

    ! ----------------------------------------------------------------------
    ! INITIALIZATION
    ! ----------------------------------------------------------------------
    
    ! Initialize EN 13383 Database (Mass in kg) - Sorted by M50
    ! Coarse / Light Gradings
    standard_gradings(1) = GradingDef("CP 45/125", 0.4_dp, 1.2_dp, 0.8_dp)
    standard_gradings(2) = GradingDef("CP 45/180", 0.4_dp, 1.2_dp, 0.8_dp) ! M50=0.8 duplicates, logic handles sort
    standard_gradings(3) = GradingDef("CP 90/180", 2.1_dp, 2.8_dp, 2.45_dp)
    standard_gradings(4) = GradingDef("CP 63/180", 1.2_dp, 3.8_dp, 2.5_dp)
    standard_gradings(5) = GradingDef("CP 90/250", 3.1_dp, 9.3_dp, 6.2_dp)
    
    ! Light Mass Armourstone (LMA)
    standard_gradings(6) = GradingDef("LMA 5-40", 10.0_dp, 20.0_dp, 15.0_dp)
    standard_gradings(7) = GradingDef("LMA 10-60", 20.0_dp, 35.0_dp, 27.5_dp)
    standard_gradings(8) = GradingDef("LMA 15-120", 35.0_dp, 60.0_dp, 47.5_dp)
    standard_gradings(9) = GradingDef("LMA 15-300", 45.0_dp, 135.0_dp, 90.0_dp)
    standard_gradings(10) = GradingDef("LMA 40-200", 80.0_dp, 120.0_dp, 100.0_dp)
    standard_gradings(11) = GradingDef("LMA 60-300", 120.0_dp, 190.0_dp, 155.0_dp)
    
    ! Heavy Mass Armourstone (HMA)
    standard_gradings(12) = GradingDef("HMA 300-1000", 540.0_dp, 690.0_dp, 615.0_dp)
    standard_gradings(13) = GradingDef("HMA 1000-3000", 1700.0_dp, 2100.0_dp, 1900.0_dp)
    standard_gradings(14) = GradingDef("HMA 3000-6000", 4200.0_dp, 4800.0_dp, 4500.0_dp)
    standard_gradings(15) = GradingDef("HMA 6000-10000", 7500.0_dp, 8500.0_dp, 8000.0_dp)
    standard_gradings(16) = GradingDef("HMA 10000-15000", 12000.0_dp, 13000.0_dp, 12500.0_dp)

    ! Bubble sort standard_gradings by M50 just to be safe (though defined mostly sorted)
    CALL sort_gradings(standard_gradings)

    ! Default Inputs
    user_inputs%Hs = 2.0_dp
    user_inputs%Tm10 = 12.0_dp
    user_inputs%h_toe = 7.0_dp
    user_inputs%slope_m = 2.0_dp
    user_inputs%foreshore_m = 30.0_dp
    user_inputs%rho_r = 2650.0_dp
    user_inputs%rho_w = 1025.0_dp
    user_inputs%P = 0.4_dp
    user_inputs%D_ratio = 0.3_dp
    user_inputs%Sd = 2.0_dp
    user_inputs%duration = 6.0_dp
    user_inputs%use_en13383 = .TRUE.
    user_inputs%manual_selection = .FALSE.

    ! ==============================================================================
    ! MAIN EXECUTION BLOCK
    ! ==============================================================================
    
    num_args = COMMAND_ARGUMENT_COUNT()

    IF (num_args >= 9) THEN
        ! CLI Mode
        CALL GET_COMMAND_ARGUMENT(1, arg_val); READ(arg_val, *) user_inputs%Hs
        CALL GET_COMMAND_ARGUMENT(2, arg_val); READ(arg_val, *) user_inputs%Tm10
        CALL GET_COMMAND_ARGUMENT(3, arg_val); READ(arg_val, *) user_inputs%h_toe
        CALL GET_COMMAND_ARGUMENT(4, arg_val); READ(arg_val, *) user_inputs%slope_m
        CALL GET_COMMAND_ARGUMENT(5, arg_val); READ(arg_val, *) user_inputs%foreshore_m
        CALL GET_COMMAND_ARGUMENT(6, arg_val); READ(arg_val, *) user_inputs%rho_r
        CALL GET_COMMAND_ARGUMENT(7, arg_val); READ(arg_val, *) user_inputs%P
        CALL GET_COMMAND_ARGUMENT(8, arg_val); READ(arg_val, *) user_inputs%D_ratio
        CALL GET_COMMAND_ARGUMENT(9, arg_val); READ(arg_val, *) user_inputs%Sd
        
        IF (num_args >= 10) THEN
            CALL GET_COMMAND_ARGUMENT(10, arg_val); READ(arg_val, *) user_inputs%duration
        END IF
        
        IF (num_args >= 11) THEN
            CALL GET_COMMAND_ARGUMENT(11, arg_val)
            CALL to_lower(arg_val)
            IF (TRIM(arg_val) == "true" .OR. TRIM(arg_val) == "1") THEN
                user_inputs%use_en13383 = .TRUE.
            ELSE
                user_inputs%use_en13383 = .FALSE.
            END IF
        END IF
    ELSE
        ! Interactive Mode
        PRINT *, ""
        PRINT *, "================================================================================"
        PRINT *, "   ROCK SLOPE STABILITY CALCULATOR"
        PRINT *, "================================================================================"
        PRINT *, "Press [Enter] to accept the default value shown in brackets."
        PRINT *, ""

        CALL get_param("Significant Wave Height Hm0 (at toe) [m]", user_inputs%Hs)
        CALL get_param("Spectral Period Tm-1,0 (at toe) [s]", user_inputs%Tm10)
        CALL get_param("Water Depth at Toe h [m]", user_inputs%h_toe)
        CALL get_param("Structure Slope (m:1) [cot alpha]", user_inputs%slope_m)
        CALL get_param("Foreshore Slope (m_f:1) [cot beta]", user_inputs%foreshore_m)
        CALL get_param("Rock Density [kg/m3]", user_inputs%rho_r)
        CALL get_param("Notional Permeability P (0.4=Permeable)", user_inputs%P)
        CALL get_param("Core/Armor Diameter Ratio (0.3=Typical)", user_inputs%D_ratio)
        CALL get_param("Design Damage Level S (2.0=Start)", user_inputs%Sd)
        CALL get_param("Storm Duration [hours]", user_inputs%duration)
        
        CALL get_bool_param("Use EN13383 Standard Grading? (True/False)", user_inputs%use_en13383)
        CALL get_bool_param("Choose stability formula instead of automatic (True/False)", user_inputs%manual_selection)
    END IF

    CALL solve(user_inputs, final_report)
    CALL generate_report_file(final_report, "output.txt")

CONTAINS

    ! --- Helpers ---

    SUBROUTINE sort_gradings(gradings)
        TYPE(GradingDef), INTENT(INOUT) :: gradings(:)
        TYPE(GradingDef) :: temp
        INTEGER :: i, j
        DO i = 1, SIZE(gradings) - 1
            DO j = i + 1, SIZE(gradings)
                IF (gradings(i)%M50 > gradings(j)%M50) THEN
                    temp = gradings(i)
                    gradings(i) = gradings(j)
                    gradings(j) = temp
                END IF
            END DO
        END DO
    END SUBROUTINE sort_gradings

    SUBROUTINE to_lower(str)
        CHARACTER(LEN=*), INTENT(INOUT) :: str
        INTEGER :: i
        DO i = 1, LEN(str)
            IF (IACHAR(str(i:i)) >= IACHAR('A') .AND. IACHAR(str(i:i)) <= IACHAR('Z')) THEN
                str(i:i) = ACHAR(IACHAR(str(i:i)) + 32)
            END IF
        END DO
    END SUBROUTINE to_lower

    SUBROUTINE get_param(prompt, val)
        CHARACTER(LEN=*), INTENT(IN) :: prompt
        REAL(dp), INTENT(INOUT) :: val
        CHARACTER(LEN=100) :: buffer_local
        INTEGER :: iostatus_val
        
        WRITE(*, '(A, " [default: ", F0.2, "]: ")', ADVANCE='NO') prompt, val
        READ(*, '(A)') buffer_local
        
        IF (LEN_TRIM(buffer_local) > 0) THEN
            READ(buffer_local, *, IOSTAT=iostatus_val) val
            IF (iostatus_val /= 0) THEN
                ! If read fails, keep default (simplistic handling)
            END IF
        END IF
    END SUBROUTINE get_param

    SUBROUTINE get_bool_param(prompt, val)
        CHARACTER(LEN=*), INTENT(IN) :: prompt
        LOGICAL, INTENT(INOUT) :: val
        CHARACTER(LEN=100) :: buffer_local
        CHARACTER(LEN=5) :: def_str
        
        IF (val) THEN
            def_str = "True"
        ELSE
            def_str = "False"
        END IF

        WRITE(*, '(A, " [default: ", A, "]: ")', ADVANCE='NO') prompt, TRIM(def_str)
        READ(*, '(A)') buffer_local
        
        IF (LEN_TRIM(buffer_local) > 0) THEN
            CALL to_lower(buffer_local)
            IF (TRIM(buffer_local) == "true" .OR. TRIM(buffer_local) == "1" .OR. &
                TRIM(buffer_local) == "t" .OR. TRIM(buffer_local) == "yes" .OR. TRIM(buffer_local) == "y") THEN
                val = .TRUE.
            ELSEIF (TRIM(buffer_local) == "false" .OR. TRIM(buffer_local) == "0" .OR. &
                    TRIM(buffer_local) == "f" .OR. TRIM(buffer_local) == "no" .OR. TRIM(buffer_local) == "n") THEN
                val = .FALSE.
            END IF
        END IF
    END SUBROUTINE get_bool_param

    ! --- Wave Mechanics ---

    FUNCTION solve_wavelength(T, h) RESULT(L)
        REAL(dp), INTENT(IN) :: T, h
        REAL(dp) :: L
        REAL(dp) :: L0, k0h, term, dL, val, f1, val_delta, f2, denom
        REAL(dp), PARAMETER :: tol = 1.0e-8_dp
        INTEGER :: iter

        IF (h <= 0.0_dp) THEN
            L = 0.0_dp
            RETURN
        END IF

        L0 = (g * T**2) / (2.0_dp * PI)
        k0h = 2.0_dp * PI * h / L0
        
        ! Initial Guess (Carvalho 2006)
        term = (6.0_dp/5.0_dp)**k0h * SQRT(k0h)
        L = L0 * TANH(term)
        
        ! Newton-Raphson
        dL = 1.0_dp
        iter = 0
        DO WHILE (ABS(dL/L) >= tol .AND. iter < 100)
            val = 2.0_dp * PI * h / L
            f1 = L - L0 * TANH(val)
            
            val_delta = 2.0_dp * PI * h / (L + tol)
            f2 = (L + tol) - L0 * TANH(val_delta)
            
            denom = f2 - f1
            IF (denom == 0.0_dp) EXIT
            
            dL = tol * f1 / denom
            L = L - dL
            iter = iter + 1
        END DO
    END FUNCTION solve_wavelength

    FUNCTION analyze_hydraulics(in, dp_params) RESULT(h)
        TYPE(Inputs), INTENT(IN) :: in
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics) :: h
        
        REAL(dp) :: k, kh, n

        h%L0 = (g * in%Tm10**2) / (2.0_dp * PI)
        h%L_toe = solve_wavelength(in%Tm10, in%h_toe)
        
        IF (in%Tm10 > 0.0_dp) THEN
            h%C = h%L_toe / in%Tm10
        ELSE
            h%C = 0.0_dp
        END IF
        
        ! Group Celerity
        IF (h%L_toe > 0.0_dp) THEN
            k = (2.0_dp * PI) / h%L_toe
        ELSE
            k = 0.0_dp
        END IF
        
        kh = k * in%h_toe
        IF (kh > 20.0_dp) THEN
            n = 0.5_dp
        ELSEIF (kh <= 0.0_dp) THEN
            n = 1.0_dp
        ELSE
            n = 0.5_dp * (1.0_dp + (2.0_dp * kh) / SINH(2.0_dp * kh))
        END IF
        h%Cg = n * h%C

        ! Steepness
        h%s_m10 = in%Hs / h%L0
        IF (h%L_toe > 0.0_dp) THEN
            h%s_local = in%Hs / h%L_toe
        ELSE
            h%s_local = 0.0_dp
        END IF
        
        IF (h%s_m10 > 0.0_dp) THEN
            h%xi_m10 = TAN(dp_params%alpha_rad) / SQRT(h%s_m10)
        ELSE
            h%xi_m10 = 0.0_dp
        END IF

        ! Relative Depth
        IF (in%Hs > 0.0_dp) THEN
            h%rel_depth = in%h_toe / in%Hs
        ELSE
            h%rel_depth = 999.0_dp
        END IF

        IF (h%xi_m10 < 0.5_dp) THEN
            h%breaker_type = "Spilling"
        ELSEIF (h%xi_m10 < 1.8_dp) THEN
            h%breaker_type = "Plunging"
        ELSEIF (h%xi_m10 < 3.0_dp) THEN
            h%breaker_type = "Surging"
        ELSE
            h%breaker_type = "Collapsing/Surging"
        END IF

        IF (h%rel_depth > 3.0_dp) THEN
            h%zone_desc = "ZONE 1: Deep to Intermediate (h/Hm0 > 3.0)"
        ELSEIF (h%rel_depth > 1.5_dp) THEN
            h%zone_desc = "ZONE 2: Shallow Water (1.5 < h/Hm0 <= 3.0)"
        ELSEIF (h%rel_depth > 0.5_dp) THEN
            h%zone_desc = "ZONE 3: Very Shallow Water (0.5 < h/Hm0 <= 1.5)"
        ELSE
            h%zone_desc = "ZONE 4: Extremely Shallow Water (h/Hm0 <= 0.5)"
        END IF
    END FUNCTION analyze_hydraulics

    ! --- Formula Functions ---

    FUNCTION calc_hudson(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: Kd_local
        
        ! Hudson Kd based on Zone Classification depending on Relative Depth
        IF (h%rel_depth > 3.0_dp) THEN
            Kd_local = 4.0_dp
        ELSEIF (h%rel_depth > 1.5_dp) THEN
            Kd_local = 3.5_dp
        ELSEIF (h%rel_depth > 0.5_dp) THEN
            Kd_local = 3.0_dp
        ELSE
            Kd_local = 2.0_dp
        END IF

        res%name = "Hudson (1959)"
        res%Dn50 = (1.27_dp * p%Hs) / (dp_params%Delta * (Kd_local * dp_params%cot_alpha)**(1.0_dp/3.0_dp))
        res%Ns = (1.27_dp * p%Hs) / (dp_params%Delta * res%Dn50)
        res%Kd = Kd_local
        res%note = "Hudson (1959) - Legacy Ref"
        res%valid = .TRUE.
    END FUNCTION calc_hudson

    FUNCTION calc_vdm_2021(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: c_pl, c_su, term_crit, xi_cr, damage_term
        CHARACTER(LEN=20) :: xi_str
        
        c_pl = 6.49_dp
        c_su = 0.97_dp
        term_crit = (c_pl / c_su) * (p%P**0.31_dp) * SQRT(TAN(dp_params%alpha_rad))
        xi_cr = term_crit**(1.0_dp / (p%P + 0.5_dp))
        damage_term = (p%Sd / SQRT(dp_params%N_waves))**0.2_dp
        
        res%name = "Van der Meer (2021)"
        
        WRITE(xi_str, '(F0.2)') xi_cr

        IF (h%xi_m10 < xi_cr) THEN
            res%Ns = c_pl * (p%P**0.18_dp) * damage_term * (h%xi_m10**(-0.5_dp))
            res%note = "Van der Meer 2021 (Plunging: xi < " // TRIM(xi_str) // ")"
        ELSE
            res%Ns = c_su * (p%P**(-0.13_dp)) * damage_term * SQRT(dp_params%cot_alpha) * (h%xi_m10**p%P)
            res%note = "Van der Meer 2021 (Surging: xi > " // TRIM(xi_str) // ")"
        END IF
        
        res%Dn50 = p%Hs / (dp_params%Delta * res%Ns)
        res%Kd = (res%Ns**3) / dp_params%cot_alpha
        res%valid = .TRUE.
    END FUNCTION calc_vdm_2021

    FUNCTION calc_van_gent_mod(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: h2_ratio, c_pl, c_su, term_crit, xi_cr, damage_term, ratio_term
        CHARACTER(LEN=20) :: h2_str
        
        ! H2%/Hs ratio estimation
        IF (h%rel_depth >= 3.0_dp) THEN
            h2_ratio = 1.4_dp
        ELSEIF (h%rel_depth < 1.5_dp) THEN
            h2_ratio = 1.2_dp + 0.2_dp * (1.5_dp - h%rel_depth) / 1.5_dp
            h2_ratio = MIN(h2_ratio, 1.4_dp)
        ELSE
            h2_ratio = 1.2_dp + 0.2_dp * (h%rel_depth - 1.5_dp) / 1.5_dp
        END IF

        c_pl = 8.4_dp
        c_su = 1.3_dp
        
        term_crit = (c_pl / c_su) * (p%P**0.31_dp) * SQRT(TAN(dp_params%alpha_rad))
        xi_cr = term_crit**(1.0_dp / (p%P + 0.5_dp))
        damage_term = (p%Sd / SQRT(dp_params%N_waves))**0.2_dp
        ratio_term = h2_ratio**(-1.0_dp)

        res%name = "Van Gent Modified (2003)"
        
        WRITE(h2_str, '(F0.2)') h2_ratio

        IF (h%xi_m10 < xi_cr) THEN
            res%Ns = 8.4_dp * (p%P**0.18_dp) * damage_term * (h%xi_m10**(-0.5_dp)) * ratio_term
            res%note = "Van Gent Modified (2003) Plunging (H2%/Hs=" // TRIM(h2_str) // ")"
        ELSE
            res%Ns = 1.3_dp * (p%P**(-0.13_dp)) * damage_term * SQRT(dp_params%cot_alpha) * (h%xi_m10**p%P) * ratio_term
            res%note = "Van Gent Modified (2003) Surging (H2%/Hs=" // TRIM(h2_str) // ")"
        END IF

        res%Dn50 = p%Hs / (dp_params%Delta * res%Ns)
        res%Kd = (res%Ns**3) / dp_params%cot_alpha
        res%valid = .TRUE.
    END FUNCTION calc_van_gent_mod

    FUNCTION calc_van_gent_simp(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: c_VG, damage_term, perm_term
        
        c_VG = 1.75_dp
        damage_term = (p%Sd / SQRT(dp_params%N_waves))**0.2_dp
        perm_term = 1.0_dp + p%D_ratio
        
        res%name = "Van Gent Simplified (2003)"
        res%Ns = c_VG * SQRT(dp_params%cot_alpha) * perm_term * damage_term
        res%Dn50 = p%Hs / (dp_params%Delta * res%Ns)
        res%Kd = (res%Ns**3) / dp_params%cot_alpha
        res%note = "Van Gent Simplified (2003)"
        res%valid = .TRUE.
    END FUNCTION calc_van_gent_simp

    FUNCTION calc_eldrup(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: c_EA1, c_EA2, damage_term, Ns_pl, cot_term, Ns_su
        
        c_EA1 = 4.5_dp
        c_EA2 = 3.1_dp
        damage_term = (p%Sd / SQRT(dp_params%N_waves))**0.2_dp
        
        Ns_pl = c_EA1 * damage_term * (1.6_dp**p%P) * (h%xi_m10**(0.4_dp*p%P - 0.67_dp))
        cot_term = MIN(dp_params%cot_alpha, 2.0_dp)**0.23_dp
        Ns_su = c_EA2 * damage_term * (p%P**0.17_dp) * cot_term
        
        res%name = "Eldrup & Andersen (2019)"
        
        IF (h%xi_m10 < 2.5_dp) THEN
            res%Ns = Ns_pl
            res%note = "Eldrup (Plunging)"
        ELSE
            res%Ns = Ns_su
            res%note = "Eldrup (Surging)"
        END IF
        
        res%Dn50 = p%Hs / (dp_params%Delta * res%Ns)
        res%Kd = (res%Ns**3) / dp_params%cot_alpha
        res%valid = .TRUE.
    END FUNCTION calc_eldrup

    FUNCTION calc_es_2020(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: c_ES1, c_ES2, m, term_Nw, term_Sd, fs_factor
        CHARACTER(LEN=32) :: suffix
        
        c_ES1 = 4.5_dp
        c_ES2 = 3.9_dp
        m = 1.0_dp / p%foreshore_m
        term_Nw = dp_params%N_waves**(-0.1_dp)
        term_Sd = p%Sd**(1.0_dp/6.0_dp)
        
        IF (h%rel_depth < 3.0_dp) THEN
            fs_factor = MAX(0.1_dp, 1.0_dp - 3.0_dp*m)
            suffix = " (Depth-Limited)"
        ELSE
            fs_factor = 1.0_dp
            suffix = " (Deep)"
        END IF
        
        res%name = "Etemad-Shahidi (2020)"
        
        IF (h%xi_m10 < 1.8_dp) THEN
            res%Ns = c_ES1 * dp_params%Cp * term_Nw * term_Sd * (h%xi_m10**(-7.0_dp/12.0_dp)) * fs_factor
            res%note = "Etemad-Shahidi (Plunging)" // TRIM(suffix)
        ELSE
            res%Ns = c_ES2 * dp_params%Cp * term_Nw * term_Sd * (h%xi_m10**(-1.0_dp/3.0_dp)) * fs_factor
            res%note = "Etemad-Shahidi (Surging)" // TRIM(suffix)
        END IF
        
        res%Dn50 = p%Hs / (dp_params%Delta * res%Ns)
        res%Kd = (res%Ns**3) / dp_params%cot_alpha
        res%valid = .TRUE.
    END FUNCTION calc_es_2020

    FUNCTION calc_mod_vg(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: c_VG_new, perm_term, s_safe, steep_term, damage_term, Kd_calc
        
        c_VG_new = 3.3_dp
        perm_term = 1.0_dp + p%D_ratio
        s_safe = MAX(h%s_m10, 0.005_dp)
        steep_term = s_safe**0.1_dp
        damage_term = (p%Sd / SQRT(dp_params%N_waves))**0.2_dp
        
        res%name = "Scaravaglione (Mod. VG 2025)"
        res%Ns = c_VG_new * SQRT(dp_params%cot_alpha) * perm_term * steep_term * damage_term
        
        Kd_calc = (res%Ns**3) / dp_params%cot_alpha
        res%note = "Scaravaglione (Mod. VG)"
        
        IF (Kd_calc > 5.0_dp) THEN
            res%Ns = (5.0_dp * dp_params%cot_alpha)**(1.0_dp/3.0_dp)
            res%note = TRIM(res%note) // " [Capped Kd=5.0]"
        END IF
        
        res%Dn50 = p%Hs / (dp_params%Delta * res%Ns)
        res%Kd = (res%Ns**3) / dp_params%cot_alpha
        res%valid = .TRUE.
    END FUNCTION calc_mod_vg

    FUNCTION calc_mod_es(p, dp_params, h) RESULT(res)
        TYPE(Inputs), INTENT(IN) :: p
        TYPE(DerivedParams), INTENT(IN) :: dp_params
        TYPE(Hydraulics), INTENT(IN) :: h
        TYPE(FormulaResult) :: res
        
        REAL(dp) :: c_ES_new, term_N, term_Sd, term_cot, s_safe, term_s, Kd_calc
        
        IF (h%xi_m10 < 1.8_dp) THEN
            res%name = "Scaravaglione (Mod. ES 2025)"
            res%Dn50 = 0.0_dp
            res%Ns = 0.0_dp
            res%Kd = 0.0_dp
            res%note = "N/A (Req. xi > 1.8)"
            res%valid = .FALSE.
            RETURN
        END IF
        
        c_ES_new = 3.55_dp
        term_N = dp_params%N_waves**(-0.1_dp)
        term_Sd = p%Sd**(1.0_dp/6.0_dp)
        term_cot = dp_params%cot_alpha**(1.0_dp/3.0_dp)
        s_safe = MAX(h%s_m10, 0.005_dp)
        term_s = s_safe**(1.0_dp/20.0_dp)
        
        res%name = "Scaravaglione (Mod. ES 2025)"
        res%Ns = c_ES_new * dp_params%Cp * term_N * term_cot * term_Sd * term_s
        
        Kd_calc = (res%Ns**3) / dp_params%cot_alpha
        res%note = "Scaravaglione (Mod. ES)"
        
        IF (Kd_calc > 5.0_dp) THEN
            res%Ns = (5.0_dp * dp_params%cot_alpha)**(1.0_dp/3.0_dp)
            res%note = TRIM(res%note) // " [Capped Kd=5.0]"
        END IF
        
        res%Dn50 = p%Hs / (dp_params%Delta * res%Ns)
        res%Kd = (res%Ns**3) / dp_params%cot_alpha
        res%valid = .TRUE.
    END FUNCTION calc_mod_es

    ! --- Design Logic ---

    FUNCTION design_layer(target_mass, target_dn, is_armor, in) RESULT(ld)
        REAL(dp), INTENT(IN) :: target_mass, target_dn
        LOGICAL, INTENT(IN) :: is_armor
        TYPE(Inputs), INTENT(IN) :: in
        TYPE(LayerDesign) :: ld
        
        REAL(dp) :: gamma_r, w_min, w_max, diff, min_diff
        REAL(dp) :: final_w_mean, final_w_min, final_w_max, final_M50
        CHARACTER(LEN=64) :: selected_name
        LOGICAL :: found
        INTEGER :: i
        REAL(dp) :: x_val, a_min, b_min, c_min, a_max, b_max, c_max, porosity

        IF (is_armor) THEN
            ld%layer_name = "Primary Armor"
        ELSE
            ld%layer_name = "Underlayer"
        END IF
        
        ld%target_M50_kg = target_mass
        ld%target_Dn_m = target_dn
        ld%target_W_kN = target_mass * g / 1000.0_dp
        
        gamma_r = in%rho_r * g / 1000.0_dp
        
        ld%design_valid = .FALSE.
        IF (in%use_en13383) THEN
            selected_name = ""
            found = .FALSE.
            IF (is_armor) THEN
                ! For Armor: Select lightest where M50 >= Target
                DO i = 1, SIZE(standard_gradings)
                    IF (standard_gradings(i)%M50 >= target_mass) THEN
                        selected_name = standard_gradings(i)%name
                        final_M50 = standard_gradings(i)%M50
                        final_w_min = standard_gradings(i)%min_mass * g / 1000.0_dp
                        final_w_max = standard_gradings(i)%max_mass * g / 1000.0_dp
                        found = .TRUE.
                        EXIT
                    END IF
                END DO
            ELSE
                ! For Underlayer: Select closest M50
                min_diff = HUGE(1.0_dp)
                DO i = 1, SIZE(standard_gradings)
                    diff = ABS(standard_gradings(i)%M50 - target_mass)
                    IF (diff < min_diff) THEN
                        min_diff = diff
                        selected_name = standard_gradings(i)%name
                        final_M50 = standard_gradings(i)%M50
                        final_w_min = standard_gradings(i)%min_mass * g / 1000.0_dp
                        final_w_max = standard_gradings(i)%max_mass * g / 1000.0_dp
                        found = .TRUE.
                    END IF
                END DO
            END IF
            
            IF (found) THEN
                ld%grading_name = selected_name
                ld%m_mean_kg = final_M50
                ld%w_min_kn = final_w_min
                ld%w_max_kn = final_w_max
                ld%w_min_kg = final_w_min * 1000.0_dp / g
                ld%w_max_kg = final_w_max * 1000.0_dp / g
                ld%w_mean_kn = ld%m_mean_kg * g / 1000.0_dp
                ld%actual_dn = (ld%w_mean_kn / gamma_r)**(1.0_dp/3.0_dp)
                ld%design_valid = .TRUE.
            ELSE
                ld%grading_name = "No Standard Fit"
                ld%design_valid = .FALSE.
            END IF
        END IF
        
        IF (.NOT. in%use_en13383 .OR. .NOT. ld%design_valid) THEN
            IF (is_armor) THEN
                ld%grading_name = "Custom Grading"
            ELSE
                ld%grading_name = "Custom Grading Underlayer"
            END IF
            
            x_val = ld%target_W_kN
            
            ! Grading Min Params
            a_min = 1.056832014477894E+00_dp
            b_min = 1.482769823574055E+00_dp
            c_min = -2.476127406338004E-01_dp
            
            ld%w_min_kn = x_val * a_min / (1.0_dp + (x_val / b_min)**c_min)
            
            ! Grading Max Params
            a_max = 1.713085676568561E+00_dp
            b_max = 2.460481255856126E+05_dp
            c_max = 1.327263214034671E-01_dp
            
            ld%w_max_kn = x_val * a_max / (1.0_dp + (x_val / b_max)**c_max)
            
            ld%w_mean_kn = x_val
            ld%m_mean_kg = target_mass
            ld%actual_dn = target_dn
            ld%w_min_kg = (ld%w_min_kn * 1000.0_dp) / g
            ld%w_max_kg = (ld%w_max_kn * 1000.0_dp) / g
            ld%design_valid = .TRUE.
        END IF
        
        ! Geometry
        ld%thickness = 2.0_dp * 1.0_dp * ld%actual_dn
        porosity = 0.30_dp
        ld%packing_density = 100.0_dp * 2.0_dp * 1.0_dp * (1.0_dp - porosity) / (ld%actual_dn**2)
    END FUNCTION design_layer

    ! --- Main Logic ---

    SUBROUTINE add_justification(report, line)
        TYPE(FullReport), INTENT(INOUT) :: report
        CHARACTER(LEN=*), INTENT(IN) :: line
        CHARACTER(LEN=300), ALLOCATABLE :: temp(:)
        INTEGER :: n
        
        IF (.NOT. ALLOCATED(report%justification)) THEN
            ALLOCATE(report%justification(1))
            n = 0
        ELSE
            n = SIZE(report%justification)
            ALLOCATE(temp(n))
            temp = report%justification
            DEALLOCATE(report%justification)
            ALLOCATE(report%justification(n + 1))
            report%justification(1:n) = temp
            DEALLOCATE(temp)
        END IF
        
        report%justification(n + 1) = line
    END SUBROUTINE add_justification

    FUNCTION find_result(report, key) RESULT(res)
        TYPE(FullReport), INTENT(IN) :: report
        CHARACTER(LEN=*), INTENT(IN) :: key
        TYPE(FormulaResult) :: res
        INTEGER :: i
        
        res%valid = .FALSE.
        res%name = ""
        
        DO i = 1, SIZE(report%comparison)
            IF (INDEX(report%comparison(i)%name, key) > 0 .AND. report%comparison(i)%valid) THEN
                res = report%comparison(i)
                RETURN
            END IF
        END DO
    END FUNCTION find_result

    SUBROUTINE solve(in, report)
        TYPE(Inputs), INTENT(IN) :: in
        TYPE(FullReport), INTENT(OUT) :: report
        
        TYPE(DerivedParams) :: dp_local
        TYPE(Hydraulics) :: h
        TYPE(FormulaResult) :: vdm, vg_mod, vg_simp, es, mod_vg, mod_es
        CHARACTER(LEN=100) :: str_buf, str_buf2
        
        report%inputs = in
        
        ! Derived Params
        dp_local%cot_alpha = in%slope_m
        dp_local%alpha_rad = ATAN(1.0_dp / in%slope_m)
        dp_local%alpha_deg = dp_local%alpha_rad * 180.0_dp / PI
        dp_local%Delta = (in%rho_r / in%rho_w) - 1.0_dp
        dp_local%Cp = (1.0_dp + in%D_ratio**0.3_dp)**0.6_dp
        dp_local%N_waves = (in%duration * 3600.0_dp) / in%Tm10
        report%derived = dp_local
        
        ! Hydraulics
        report%hydro = analyze_hydraulics(in, dp_local)
        
        ! Allocate comparison array (8 formulas now)
        ALLOCATE(report%comparison(8))
        
        report%comparison(1) = calc_hudson(in, dp_local, report%hydro)
        report%comparison(2) = calc_vdm_2021(in, dp_local, report%hydro)
        report%comparison(3) = calc_van_gent_mod(in, dp_local, report%hydro)
        report%comparison(4) = calc_van_gent_simp(in, dp_local, report%hydro)
        report%comparison(5) = calc_eldrup(in, dp_local, report%hydro)
        report%comparison(6) = calc_es_2020(in, dp_local, report%hydro)
        report%comparison(7) = calc_mod_vg(in, dp_local, report%hydro)
        report%comparison(8) = calc_mod_es(in, dp_local, report%hydro)
        
        ! Recommendation Logic
        vdm = find_result(report, "Van der Meer (2021)")
        vg_mod = find_result(report, "Van Gent Modified")
        vg_simp = find_result(report, "Van Gent Simplified")
        es = find_result(report, "Etemad-Shahidi")
        mod_vg = find_result(report, "Mod. VG")
        mod_es = find_result(report, "Mod. ES")
        
        ! Initialize justification list
        IF (ALLOCATED(report%justification)) DEALLOCATE(report%justification)
        
        IF (report%hydro%rel_depth > 3.0_dp) THEN
            ! ZONE 1: Deep
            report%recommended = vdm
            CALL add_justification(report, "### 1. HYDRAULIC CONTEXT: Deep / Intermediate Water")
            WRITE(str_buf, '(F0.2)') report%hydro%rel_depth
            CALL add_justification(report, "   The structure is located in deep water relative to the wave height (h/Hm0 = " // &
                TRIM(str_buf) // " > 3.0).")
            CALL add_justification(report, "   In this regime, the wave height distribution strictly follows the **Rayleigh distribution**.")
            CALL add_justification(report, "   Key characteristics:")
            CALL add_justification(report, "     * The ratio H2%/Hs is constant at approximately 1.4.")
            CALL add_justification(report, "     * Wave breaking is limited to whitecapping or direct interaction with the armor layer.")
            CALL add_justification(report, "     * The spectral shape is standard (JONSWAP/Pierson-Moskowitz), and energy transfer to low-frequencies is minimal.")
            CALL add_justification(report, "")
            CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
            
            ! Van der Meer
            CALL add_justification(report, "   **A. Van der Meer (2021 Rewritten) [RECOMMENDED]**")
            CALL add_justification(report, "      * **Advantages:** This formula is the modernized industry standard.")
            CALL add_justification(report, "        Van der Meer (2021) rewrote the original formula to use the spectral period (Tm-1,0),")
            CALL add_justification(report, "        eliminating the influence of spectral shape.")
            CALL add_justification(report, "        Van der Meer et al. (2024) confirmed its validity for h/Hm0 > 1.5, preferring Hm0 over H1/3 for nonlinear waves.")
            CALL add_justification(report, "      * **Physics:** It correctly assumes a Rayleigh distribution of wave heights, aligning with the actual deep-water statistics.")

            ! Van Gent Mod
            CALL add_justification(report, "   **B. Van Gent Modified (2003)**")
            CALL add_justification(report, "      * **Context:** This formula incorporates the ratio H2%/Hs.")
            CALL add_justification(report, "        In deep water, with H2%/Hs = 1.4, this formula essentially converges closely with the Van der Meer predictions.")
            CALL add_justification(report, "        However, its specific calibration was focused on the effects of shallow foreshores.")

            ! Etemad-Shahidi
            CALL add_justification(report, "   **C. Etemad-Shahidi (2020)**")
            CALL add_justification(report, "      * **Comparison:** Etemad-Shahidi (2020) provides a robust formula validated for both deep and shallow water.")
            CALL add_justification(report, "        It introduces a physical permeability parameter (D_core/D_armor) to replace the nominal P factor,")
            CALL add_justification(report, "        reducing uncertainty. However, Van der Meer remains the primary standard for deep water.")

            IF (es%valid .AND. vdm%valid .AND. vdm%Dn50 > 0.0_dp) THEN
                IF (ABS(vdm%Dn50 - es%Dn50) / vdm%Dn50 > 0.10_dp) THEN
                    CALL add_justification(report, "      * **Note on Divergence:** The result deviates from Van der Meer here. This typically occurs in the 'transition zone'")
                    CALL add_justification(report, "        of the surf similarity parameter (xi approx 2.0 - 4.5). Etemad-Shahidi transitions to Surging physics")
                    CALL add_justification(report, "        earlier (xi > 1.8), predicting higher stability, whereas Van der Meer maintains Plunging physics")
                    CALL add_justification(report, "        (lower stability) until a higher critical threshold. Van der Meer is more conservative here.")
                ELSE
                    CALL add_justification(report, "        It typically converges with Van der Meer here.")
                END IF
            END IF

            CALL add_justification(report, "")
            CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
            CALL add_justification(report, "   **Use [Van der Meer (2021 Rewritten)]**.")
            CALL add_justification(report, "   It provides the most theoretically consistent result for non-depth-limited waves.")
            
            IF (es%valid .AND. vdm%valid) THEN
                 WRITE(str_buf, '(F0.3)') es%Dn50
                 WRITE(str_buf2, '(F0.3)') ABS(vdm%Dn50 - es%Dn50)
                 CALL add_justification(report, "   *Verification:* Etemad-Shahidi yields Dn50 = " // TRIM(str_buf) // &
                    "m (Difference: " // TRIM(str_buf2) // "m).")
            END IF

        ELSEIF (report%hydro%rel_depth > 1.5_dp) THEN
            ! ZONE 2: Shallow
            report%recommended = vdm
            WRITE(str_buf, '(F0.2)') report%hydro%rel_depth
            CALL add_justification(report, "### 1. HYDRAULIC CONTEXT: Shallow Water (Transition Zone)")
            CALL add_justification(report, "   The structure is in the transition zone (1.5 < h/Hm0 = " // TRIM(str_buf) // &
                " <= 3.0).")
            CALL add_justification(report, "   Key characteristics:")
            CALL add_justification(report, "     * **Spectral Truncation:** The largest waves in the spectrum break on the foreshore.")
            CALL add_justification(report, "     * **Distribution Shift:** The wave height distribution deviates from Rayleigh; H2%/Hm0 drops below 1.4.")
            CALL add_justification(report, "     * **Shoaling:** Significant shoaling modifies the wave shape before impact, creating peaked crests and flat troughs.")
            CALL add_justification(report, "")
            CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
            
            ! VdM
            CALL add_justification(report, "   **A. Van der Meer (2021)**")
            CALL add_justification(report, "      * **Advantages:** Van der Meer et al. (2024) extensively re-analyzed shallow water data and concluded")
            CALL add_justification(report, "        that the rewritten Van der Meer formula (using Tm-1,0) is valid down to h/Hm0 = 1.5.")
            CALL add_justification(report, "        It performs reasonably well, with slightly less reliability in the 1.0 < h/Hm0 < 1.5 range.")
            CALL add_justification(report, "      * **Note:** For nonlinear waves in this zone, using Hm0 is preferred over H1/3 for nonlinear waves to avoid deviations.")

            ! VG Mod
            CALL add_justification(report, "   **B. Van Gent Modified (2003)**")
            CALL add_justification(report, "      * **Constraint:** This formula explicitly relies on the ratio H2%/Hs. Research by Van der Meer et al. (2024)")
            CALL add_justification(report, "        highlights that predicting H2% accurately in this transition zone (where the ratio dips to ~1.2)")
            CALL add_justification(report, "        is notoriously inaccurate without physical modeling. The formula is valid, but the input uncertainty is high.")

            ! VG Simp
            CALL add_justification(report, "   **C. Van Gent et al. (2003) Simplified**")
            CALL add_justification(report, "      * **Context:** This formula was specifically derived for shallow foreshores.")
            CALL add_justification(report, "        However, Van der Meer et al. (2024) found that the simplified formula often does not match the data")
            CALL add_justification(report, "        in the surging domain as well as the rewritten Van der Meer formula.")
            
            CALL add_justification(report, "")
            CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
            CALL add_justification(report, "   **Use [Van der Meer (2021 Rewritten)]**.")
            CALL add_justification(report, "   Recent research (2024) confirms its validity in this depth range (h/Hm0 > 1.5), favoring it over simplified methods")
            CALL add_justification(report, "   due to the uncertainties in predicting H2% required for the Van Gent Modified formula.")

        ELSEIF (report%hydro%rel_depth > 0.5_dp) THEN
            ! ZONE 3: Very Shallow
            WRITE(str_buf, '(F0.2)') report%hydro%rel_depth
            
            IF (mod_es%valid) THEN
                report%recommended = mod_es
                CALL add_justification(report, "### 1. HYDRAULIC CONTEXT: Very Shallow Water (Surf Zone)")
                CALL add_justification(report, "   The structure is in the surf zone (0.5 < h/Hm0 = " // TRIM(str_buf) // &
                    " <= 1.5).")
                CALL add_justification(report, "   Key characteristics:")
                CALL add_justification(report, "     * **Severe Breaking:** Waves are constantly breaking.")
                CALL add_justification(report, "     * **Saturation:** Wave height is depth-limited (H ~ 0.5h). Increasing offshore energy does not increase load.")
                CALL add_justification(report, "     * **Infragravity Dominance:** Scaravaglione et al. (2025) and VdM (2024) note that infragravity waves")
                CALL add_justification(report, "       begin to dominate the spectrum, causing Tm-1,0 to increase massively (up to 4x).")
                CALL add_justification(report, "     * **Formula Deviation:** Standard formulas fail here because the stability curves flatten out (Horizontal Trend).")
                CALL add_justification(report, "")
                CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
                
                CALL add_justification(report, "   **A. Standard Formulas (VdM, VG, Standard ES)**")
                CALL add_justification(report, "      * **Failure Mode:** Scaravaglione et al. (2025) demonstrated that these formulas fail to converge here.")
                CALL add_justification(report, "        Using the inflated spectral period results in over-predicted stability (for surging) or under-predicted (for plunging).")
                
                CALL add_justification(report, "   **B. Scaravaglione (Modified ES 2025) [RECOMMENDED]**")
                CALL add_justification(report, "      * **Advantages:** This formula explicitly **decouples** the wave steepness term (s^0.05) from the structure slope term.")
                CALL add_justification(report, "      * **Physics:** It is calibrated specifically for surging/bore-like waves in the surf zone using new coefficients (c_ES,new=3.55).")
                CALL add_justification(report, "      * **Safety Cap Applied:** The formula has a weak dependence on steepness. For very long swells, it may predict")
                CALL add_justification(report, "        theoretically high stability (Kd > 10). This system has capped Kd at 5.0 to ensure physical realism.")
                
                CALL add_justification(report, "")
                CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
                CALL add_justification(report, "   **Use [Scaravaglione (Modified ES 2025)]**.")
                CALL add_justification(report, "   This represents the state-of-the-art for broken waves in very shallow water, correcting the overestimation of damage.")
            ELSE
                report%recommended = vg_simp
                CALL add_justification(report, "### 1. HYDRAULIC CONTEXT: Very Shallow Water (Surf Zone) (Plunging)")
                CALL add_justification(report, "   The structure is in the surf zone, but conditions are **Plunging** (xi < 1.8).")
                CALL add_justification(report, "   The Modified ES formula is only calibrated for surging/bore conditions.")
                CALL add_justification(report, "")
                CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
                CALL add_justification(report, "   **Use [Van Gent Simplified (2003)]**.")
                CALL add_justification(report, "   It acts as a robust fallback. Caution is advised as damage may be underpredicted for impermeable structures.")
            END IF
        ELSE
            ! ZONE 4: Swash
            report%recommended = mod_vg
            WRITE(str_buf, '(F0.2)') report%hydro%rel_depth
            CALL add_justification(report, "### 1. HYDRAULIC CONTEXT: Extremely Shallow Water (Swash Zone)")
            CALL add_justification(report, "   The structure is located effectively in the **Swash Zone** (h/Hm0 = " // &
                TRIM(str_buf) // " <= 0.5).")
            CALL add_justification(report, "   Key characteristics:")
            CALL add_justification(report, "     * **Aeration:** High air entrainment reduces the effective fluid density and buoyancy of the rocks.")
            CALL add_justification(report, "     * **Impact:** Wave impact is characterized by a high-velocity turbulent bore.")
            CALL add_justification(report, "     * **Hydrostatics:** The hydrostatic cushioning effect is negligible.")
            CALL add_justification(report, "")
            CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
            
            CALL add_justification(report, "   **A. Standard Van Gent (2003)**")
            CALL add_justification(report, "      * **Disadvantages:** Using the standard coefficient (1.75) here is **Unsafe**.")
            CALL add_justification(report, "        Scaravaglione et al. (2025) showed that stability is significantly lower than predicted by intermediate-depth formulas")
            CALL add_justification(report, "        due to the lack of buoyancy and intense turbulence.")
            
            CALL add_justification(report, "   **B. Scaravaglione (Modified VG 2025) [RECOMMENDED]**")
            CALL add_justification(report, "      * **Advantages:** This formula uses a recalibrated coefficient (C_VG = 3.3 instead of 1.75).")
            CALL add_justification(report, "      * **Physics:** It explicitly accounts for the increased instability in the swash zone, correcting the underestimation")
            CALL add_justification(report, "        of damage by the original VG formula in this specific regime.")
            
            CALL add_justification(report, "")
            CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
            CALL add_justification(report, "   **Use [Scaravaglione (Modified VG 2025)]**.")
            CALL add_justification(report, "   It provides the necessary safety margin for swash zone instability where standard formulas fail.")
        END IF

        ! Interactive Manual Selection Logic
        IF (in%manual_selection) THEN
            PRINT *, ""
            PRINT *, "================================================================================"
            PRINT *, "   MANUAL SELECTION MODE"
            PRINT *, "================================================================================"
            
            DO i = 1, SIZE(report%comparison)
                IF (report%comparison(i)%valid) THEN
                    WRITE(*, '(A, I0, A, A)') "   [", i, "] ", TRIM(report%comparison(i)%name)
                END IF
            END DO
            
            PRINT *, ""
            WRITE(*, '(A)', ADVANCE='NO') "Enter the number of your preferred formula: "
            READ(*, '(A)') str_buf
            
            READ(str_buf, *, IOSTAT=iostat) i
            IF (iostat == 0 .AND. i >= 1 .AND. i <= SIZE(report%comparison)) THEN
                IF (report%comparison(i)%valid) THEN
                    report%recommended = report%comparison(i)
                    str_buf = ""
                    str_buf = CHAR(10) // "[MANUAL OVERRIDE] User switched selection to: " // TRIM(report%recommended%name)
                    CALL to_lower(str_buf) ! Note: to_upper implementation omitted, just printing as is or lowercase is acceptable fallback or adding to_upper helper.
                    ! Just using name as is for CLI consistency
                    str_buf = CHAR(10) // "[MANUAL OVERRIDE] User switched selection to: " // TRIM(report%recommended%name)
                    PRINT '(A)', TRIM(str_buf)
                    CALL add_justification(report, TRIM(str_buf))
                ELSE
                    PRINT *, "Invalid selection (formula not valid). Using automatic."
                END IF
            ELSE
                PRINT *, "Invalid selection. Using automatic."
            END IF
        END IF

        ! Calculate Layers
        report%armor_layer = design_layer(in%rho_r * (report%recommended%Dn50**3), report%recommended%Dn50, .TRUE., in)
        report%underlayer = design_layer(report%armor_layer%m_mean_kg / 10.0_dp, &
            (report%armor_layer%m_mean_kg / 10.0_dp / in%rho_r)**(1.0_dp/3.0_dp), .FALSE., in)

    END SUBROUTINE solve

    SUBROUTINE generate_report_file(report, filepath)
        TYPE(FullReport), INTENT(IN) :: report
        CHARACTER(LEN=*), INTENT(IN) :: filepath
        
        INTEGER :: u, iostatus_val, i
        CHARACTER(LEN=30) :: name_field
        CHARACTER(LEN=120) :: separator
        CHARACTER(LEN=120) :: separator_long
        TYPE(FormulaResult) :: res
        REAL(dp) :: mass
        
        separator = "-----------------------------------------------------------------------------------------------"
        separator_long = "--------------------------------------------------------------------------------" // &
                         "-----------------------------------"
        
        OPEN(NEWUNIT=u, FILE=filepath, STATUS='REPLACE', ACTION='WRITE', IOSTAT=iostatus_val)
        IF (iostatus_val /= 0) THEN
            PRINT *, "Error saving file."
            RETURN
        END IF
        
        WRITE(u, '(A)') "ROCK SLOPE STABILITY CALCULATOR"
        WRITE(u, '(A)') ""
        WRITE(u, '(A)') "==============================================================================================="
        WRITE(u, '(A)') "   1. DESIGN INPUT PARAMETERS"
        WRITE(u, '(A)') "==============================================================================================="
        WRITE(u, '(A35, " | ", A10, " | ", A)') "PARAMETER", "VALUE", "UNIT"
        WRITE(u, '(A)') TRIM(separator)
        WRITE(u, '(A35, " | ", F10.2, " | ", A)') "Significant Wave Height (Hs)", report%inputs%Hs, "m"
        WRITE(u, '(A35, " | ", F10.2, " | ", A)') "Spectral Period (Tm-1,0)", report%inputs%Tm10, "s"
        WRITE(u, '(A35, " | ", F10.2, " | ", A)') "Water Depth (h_toe)", report%inputs%h_toe, "m"
        WRITE(u, '(A35, " | 1:", F8.2, " | ", A)')  "Structure Slope", report%inputs%slope_m, "(V:H)"
        WRITE(u, '(A35, " | 1:", F8.2, " | ", A)')  "Foreshore Slope", report%inputs%foreshore_m, "(V:H)"
        WRITE(u, '(A35, " | ", F10.2, " | ", A)') "Permeability (P_notional)", report%inputs%P, "(-)"
        WRITE(u, '(A35, " | ", F10.2, " | ", A)') "Physical Permeability (Cp)", report%derived%Cp, "(-)"
        WRITE(u, '(A35, " | ", F10.1, " | ", A)') "Damage Level (S)", report%inputs%Sd, "(-)"
        WRITE(u, '(A)') TRIM(separator)
        WRITE(u, '(A)') ""
        
        WRITE(u, '(A)') "==============================================================================================="
        WRITE(u, '(A)') "   2. CALCULATED HYDRAULIC PARAMETERS"
        WRITE(u, '(A)') "==============================================================================================="
        WRITE(u, '(A40, " | ", A10, " | ", A)') "PARAMETER", "VALUE", "UNIT"
        WRITE(u, '(A)') TRIM(separator)
        WRITE(u, '(A40, " | ", F10.2, " | ", A)') "Deep Water Wavelength (L0)", report%hydro%L0, "m"
        WRITE(u, '(A40, " | ", F10.2, " | ", A)') "Wavelength at Toe (L_toe)", report%hydro%L_toe, "m"
        WRITE(u, '(A40, " | ", F10.2, " | ", A)') "Wave Celerity at Toe (C)", report%hydro%C, "m/s"
        WRITE(u, '(A40, " | ", F10.2, " | ", A)') "Group Celerity at Toe (Cg)", report%hydro%Cg, "m/s"
        WRITE(u, '(A40, " | ", F10.5, " | ", A)') "Deep Water Steepness (s_m-1,0)", report%hydro%s_m10, "(-)"
        WRITE(u, '(A40, " | ", F10.5, " | ", A)') "Local Wave Steepness (s_local)", report%hydro%s_local, "(-)"
        WRITE(u, '(A40, " | ", F10.2, " | ", A)') "Surf Similarity (xi_m-1,0)", report%hydro%xi_m10, "(-)"
        WRITE(u, '(A40, " | ", A10, " | ", A)')   "Breaker Type (Visual/Physical)", report%hydro%breaker_type, "(-)"
        WRITE(u, '(A40, " | ", F10.2, " | ", A)') "Relative Depth (h/Hm0)", report%hydro%rel_depth, "(-)"
        WRITE(u, '(A40, " | ", F10.3, " | ", A)') "Relative Depth (h/L0)", &
            report%hydro%rel_depth * report%hydro%s_m10, "(-)"
        WRITE(u, '(A40, " | ", A40)')             "Hydraulic Zone", report%hydro%zone_desc
        WRITE(u, '(A)') TRIM(separator)
        WRITE(u, '(A)') ""
        
        WRITE(u, '(A)') "==============================================================================================="
        WRITE(u, '(A)') "   3. FORMULA SELECTION & JUSTIFICATION"
        WRITE(u, '(A)') "==============================================================================================="
        WRITE(u, '(A)') "COMPARISON OF RESULTS:"
        WRITE(u, '(A30, " | ", A8, " | ", A10, " | ", A10, " | ", A8, " | ", A)') &
            "Method", "Ns (-)", "Dn50 (m)", "M50 (kg)", "Kd_eq (-)", "NOTES"
        WRITE(u, '(A)') TRIM(separator_long)
        
        DO i = 1, SIZE(report%comparison)
            res = report%comparison(i)
            IF (res%valid) THEN
                mass = report%inputs%rho_r * (res%Dn50**3)
                ! Use fixed width char buffer for left alignment
                name_field = res%name 
                WRITE(u, '(A30, " | ", F8.4, " | ", F10.3, " | ", F10.0, " | ", F8.2, " | ", A)') &
                    name_field, res%Ns, res%Dn50, mass, res%Kd, TRIM(res%note)
            END IF
        END DO
        WRITE(u, '(A)') ""
        WRITE(u, '(A)') "JUSTIFICATION & ANALYSIS:"
        IF (ALLOCATED(report%justification)) THEN
            DO i = 1, SIZE(report%justification)
                WRITE(u, '(A)') TRIM(report%justification(i))
            END DO
        END IF
        WRITE(u, '(A)') TRIM(separator)
        WRITE(u, '(A)') ""
        
        WRITE(u, '(A)') "==============================================================================================="
        IF (report%inputs%use_en13383) THEN
            WRITE(u, '(A)') "   4. ROCK ARMOUR LAYER DESIGN (EN 13383 Standard)"
        ELSE
            WRITE(u, '(A)') "   4. ROCK ARMOUR LAYER DESIGN (Custom Grading)"
        END IF
        WRITE(u, '(A)') "==============================================================================================="
        WRITE(u, '(A)') "PRIMARY ARMOR LAYER"
        WRITE(u, '(A, F0.2, A)') "   Theoretical Required W    : ", report%armor_layer%target_W_kN, " kN"
        WRITE(u, '(A, F0.0, A)') "   Theoretical Required M50  : ", report%armor_layer%target_M50_kg, " kg"
        WRITE(u, '(A, F0.3, A)') "   Theoretical Required Dn50 : ", report%armor_layer%target_Dn_m, " m"
        WRITE(u, '(A)') "----------------------------------------"
        
        IF (.NOT. report%armor_layer%design_valid .AND. report%inputs%use_en13383) THEN
             WRITE(u, '(A)') "   [WARNING] No standard EN13383 grading found for this mass."
        ELSE
             WRITE(u, '(A, A)')       "   Adopted rock grading                : ", TRIM(report%armor_layer%grading_name)
             WRITE(u, '(A, F0.2, A, F0.0, A)') "   Grading Min (Lower Limit)           : ", report%armor_layer%w_min_kn, &
                " kN (", report%armor_layer%w_min_kg, " kg)"
             WRITE(u, '(A, F0.2, A, F0.0, A)') "   Grading Max (Upper Limit)           : ", report%armor_layer%w_max_kn, &
                " kN (", report%armor_layer%w_max_kg, " kg)"
             WRITE(u, '(A, F0.0, A)') "   Representative M50                  : ", report%armor_layer%m_mean_kg, " kg"
             WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn_rock)          : ", report%armor_layer%actual_dn, " m"
             WRITE(u, '(A, F0.2, A)') "   Double Layer Thickness              : ", report%armor_layer%thickness, " m"
             WRITE(u, '(A, F0.2)')    "   Packing Density [rocks/100m2]       : ", report%armor_layer%packing_density
        END IF
        
        WRITE(u, '(A)') TRIM(separator)
        WRITE(u, '(A)') "UNDERLAYER (FILTER LAYER)"
        WRITE(u, '(A, F0.3, A)') "   Target Weight (M50 / 10)  : ", report%underlayer%target_W_kN, " kN"
        WRITE(u, '(A, F0.1, A)') "   Target Mass (M50 / 10)    : ", report%underlayer%target_M50_kg, " kg"
        WRITE(u, '(A)') "----------------------------------------"
        
        IF (.NOT. report%underlayer%design_valid .AND. report%inputs%use_en13383) THEN
             WRITE(u, '(A)') "   [WARNING] No suitable standard underlayer grading found."
        ELSE
             WRITE(u, '(A, A)')       "   Adopted rock grading                : ", TRIM(report%underlayer%grading_name)
             WRITE(u, '(A, F0.2, A, F0.0, A)') "   Grading Min (Lower Limit)           : ", report%underlayer%w_min_kn, &
                " kN (", report%underlayer%w_min_kg, " kg)"
             WRITE(u, '(A, F0.2, A, F0.0, A)') "   Grading Max (Upper Limit)           : ", report%underlayer%w_max_kn, &
                " kN (", report%underlayer%w_max_kg, " kg)"
             WRITE(u, '(A, F0.1, A)') "   Representative M50                  : ", report%underlayer%m_mean_kg, " kg"
             WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn_rock)          : ", report%underlayer%actual_dn, " m"
             WRITE(u, '(A, F0.2, A)') "   Double Layer Thickness              : ", report%underlayer%thickness, " m"
             WRITE(u, '(A, F0.2)')    "   Packing Density [rocks/100m2]       : ", report%underlayer%packing_density
        END IF
        WRITE(u, '(A)') "==============================================================================================="

        CLOSE(u)
        
        PRINT *, ""
        PRINT *, "[System] Report saved to ", TRIM(filepath)
    END SUBROUTINE generate_report_file

END PROGRAM RockSlopeCalculator