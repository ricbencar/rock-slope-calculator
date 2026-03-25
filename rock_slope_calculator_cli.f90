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
!      rock_slope_calculator_cli.exe [Hs] [Tm] [h_toe] [slope] [fore] [rho_r] [P] [D_ratio] [Sd] [Duration] [UseEN13383] [CustomFamily]
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
        REAL(dp) :: nll_kg   ! Nominal Lower Limit
        REAL(dp) :: nul_kg   ! Nominal Upper Limit
        REAL(dp) :: M50      ! Representative M50 = 0.5 * (NLL + NUL)
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
        CHARACTER(LEN=8) :: custom_family ! AUTO, HMA, LMA, or CP when using custom grading
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
        CHARACTER(LEN=64) :: breaker_type
        CHARACTER(LEN=128) :: zone_desc
    END TYPE Hydraulics

    TYPE :: FormulaResult
        CHARACTER(LEN=64) :: name
        REAL(dp) :: Dn50
        REAL(dp) :: Ns
        REAL(dp) :: Kd
        CHARACTER(LEN=160) :: note
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
        
        LOGICAL :: used_custom_interpolation
        CHARACTER(LEN=8) :: custom_family
        REAL(dp) :: custom_ratio_nul_nll
        CHARACTER(LEN=160) :: custom_ratio_note
        
        LOGICAL :: design_valid ! If a valid grading/design was found
    END TYPE LayerDesign

    TYPE :: FullReport
        TYPE(Inputs) :: inputs
        TYPE(DerivedParams) :: derived
        TYPE(Hydraulics) :: hydro
        TYPE(FormulaResult), ALLOCATABLE :: comparison(:)
        TYPE(FormulaResult) :: recommended
        CHARACTER(LEN=300), ALLOCATABLE :: justification(:)
        CHARACTER(LEN=200) :: manual_override_message
        TYPE(LayerDesign) :: armor_layer
        TYPE(LayerDesign) :: underlayer
    END TYPE FullReport

    ! ----------------------------------------------------------------------
    ! VARIABLES
    ! ----------------------------------------------------------------------
    TYPE(GradingDef) :: standard_gradings(17)
    TYPE(Inputs) :: user_inputs
    TYPE(FullReport) :: final_report
    INTEGER :: num_args
    CHARACTER(LEN=100) :: arg_val
    CHARACTER(LEN=100) :: buffer
    INTEGER :: i, iostat

    ! ----------------------------------------------------------------------
    ! INITIALIZATION
    ! ----------------------------------------------------------------------
    
    ! Initialize EN 13383 Database 
    ! Logic: NLL and NUL taken from the revised CP grading limits where applicable,
    ! and from the category name for LMA/HMA.
    ! M50 calculated as 0.5 * (NLL + NUL).
    
    standard_gradings(1) = GradingDef("CP32/90", 0.868_dp, 19.319_dp, 10.0935_dp)
    standard_gradings(2) = GradingDef("CP45/125", 2.415_dp, 51.758_dp, 27.0865_dp)
    standard_gradings(3) = GradingDef("CP63/180", 6.626_dp, 154.548_dp, 80.587_dp)
    standard_gradings(4) = GradingDef("CP90/250", 19.319_dp, 414.063_dp, 216.691_dp)
    standard_gradings(5) = GradingDef("CP45/180", 2.415_dp, 154.548_dp, 78.4815_dp)
    standard_gradings(6) = GradingDef("CP90/180", 19.319_dp, 154.548_dp, 86.9335_dp)
    
    ! Light Mass Armour (LMA) - NLL/NUL from name
    standard_gradings(7) = GradingDef("LMA 5-40", 5.0_dp, 40.0_dp, 22.5_dp)
    standard_gradings(8) = GradingDef("LMA 10-60", 10.0_dp, 60.0_dp, 35.0_dp)
    standard_gradings(9) = GradingDef("LMA 15-120", 15.0_dp, 120.0_dp, 67.5_dp)
    standard_gradings(10) = GradingDef("LMA 40-200", 40.0_dp, 200.0_dp, 120.0_dp)
    standard_gradings(11) = GradingDef("LMA 60-300", 60.0_dp, 300.0_dp, 180.0_dp)
    standard_gradings(12) = GradingDef("LMA 15-300", 15.0_dp, 300.0_dp, 157.5_dp)
    
    ! Heavy Mass Armour (HMA) - NLL/NUL from name
    standard_gradings(13) = GradingDef("HMA 300-1000", 300.0_dp, 1000.0_dp, 650.0_dp)
    standard_gradings(14) = GradingDef("HMA 1000-3000", 1000.0_dp, 3000.0_dp, 2000.0_dp)
    standard_gradings(15) = GradingDef("HMA 3000-6000", 3000.0_dp, 6000.0_dp, 4500.0_dp)
    standard_gradings(16) = GradingDef("HMA 6000-10000", 6000.0_dp, 10000.0_dp, 8000.0_dp)
    standard_gradings(17) = GradingDef("HMA 10000-15000", 10000.0_dp, 15000.0_dp, 12500.0_dp)

    ! Sort by M50 for selection logic
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
    user_inputs%custom_family = 'AUTO'

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

        IF (num_args >= 12) THEN
            CALL GET_COMMAND_ARGUMENT(12, arg_val)
            CALL to_upper(arg_val)
            arg_val = ADJUSTL(arg_val)
            IF (TRIM(arg_val) == 'HMA' .OR. TRIM(arg_val) == 'LMA' .OR. TRIM(arg_val) == 'CP' .OR. &
                TRIM(arg_val) == 'AUTO') THEN
                user_inputs%custom_family = TRIM(arg_val)
            ELSE
                user_inputs%custom_family = 'AUTO'
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
        IF (.NOT. user_inputs%use_en13383) THEN
            CALL get_string_param("Custom grading family [AUTO/HMA/LMA/CP]", user_inputs%custom_family)
            CALL to_upper(user_inputs%custom_family)
            IF (TRIM(user_inputs%custom_family) /= 'HMA' .AND. TRIM(user_inputs%custom_family) /= 'LMA' .AND. &
                TRIM(user_inputs%custom_family) /= 'CP' .AND. TRIM(user_inputs%custom_family) /= 'AUTO') THEN
                user_inputs%custom_family = 'AUTO'
            END IF
        END IF
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

    SUBROUTINE to_upper(str)
        CHARACTER(LEN=*), INTENT(INOUT) :: str
        INTEGER :: i
        DO i = 1, LEN(str)
            IF (IACHAR(str(i:i)) >= IACHAR('a') .AND. IACHAR(str(i:i)) <= IACHAR('z')) THEN
                str(i:i) = ACHAR(IACHAR(str(i:i)) - 32)
            END IF
        END DO
    END SUBROUTINE to_upper

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

    SUBROUTINE get_string_param(prompt, val)
        CHARACTER(LEN=*), INTENT(IN) :: prompt
        CHARACTER(LEN=*), INTENT(INOUT) :: val
        CHARACTER(LEN=100) :: buffer_local

        WRITE(*, '(A, " [default: ", A, "]: ")', ADVANCE='NO') prompt, TRIM(val)
        READ(*, '(A)') buffer_local

        IF (LEN_TRIM(buffer_local) > 0) THEN
            val = ADJUSTL(buffer_local)
        END IF
    END SUBROUTINE get_string_param

    FUNCTION pad_right(str, width) RESULT(out)
        CHARACTER(LEN=*), INTENT(IN) :: str
        INTEGER, INTENT(IN) :: width
        CHARACTER(LEN=width) :: out
        INTEGER :: n
        out = " "
        n = MIN(LEN_TRIM(str), width)
        IF (n > 0) out(1:n) = str(1:n)
    END FUNCTION pad_right

    FUNCTION real_str(val, prec) RESULT(out)
        REAL(dp), INTENT(IN) :: val
        INTEGER, INTENT(IN) :: prec
        CHARACTER(LEN=:), ALLOCATABLE :: out
        CHARACTER(LEN=64) :: tmp
        SELECT CASE (prec)
        CASE (0)
            WRITE(tmp, '(F12.0)') val
        CASE (1)
            WRITE(tmp, '(F12.1)') val
        CASE (2)
            WRITE(tmp, '(F12.2)') val
        CASE (3)
            WRITE(tmp, '(F12.3)') val
        CASE (4)
            WRITE(tmp, '(F12.4)') val
        CASE (5)
            WRITE(tmp, '(F12.5)') val
        CASE DEFAULT
            WRITE(tmp, '(ES16.8)') val
        END SELECT
        out = TRIM(ADJUSTL(tmp))
    END FUNCTION real_str

    FUNCTION int_str(i) RESULT(out)
        INTEGER, INTENT(IN) :: i
        CHARACTER(LEN=:), ALLOCATABLE :: out
        CHARACTER(LEN=32) :: tmp
        WRITE(tmp, '(I0)') i
        out = TRIM(tmp)
    END FUNCTION int_str

    FUNCTION fmt_row35(label, val, unit) RESULT(out)
        CHARACTER(LEN=*), INTENT(IN) :: label, val, unit
        CHARACTER(LEN=:), ALLOCATABLE :: out
        out = pad_right(label, 35) // ' | ' // pad_right(val, 10) // ' | ' // unit
    END FUNCTION fmt_row35

    FUNCTION fmt_row40(label, val, unit) RESULT(out)
        CHARACTER(LEN=*), INTENT(IN) :: label, val, unit
        CHARACTER(LEN=:), ALLOCATABLE :: out
        out = pad_right(label, 40) // ' | ' // pad_right(val, 10) // ' | ' // unit
    END FUNCTION fmt_row40

    FUNCTION fmt_formula_row(name, ns, dn, mass, kd, note) RESULT(out)
        CHARACTER(LEN=*), INTENT(IN) :: name, note
        REAL(dp), INTENT(IN) :: ns, dn, kd
        INTEGER, INTENT(IN) :: mass
        CHARACTER(LEN=:), ALLOCATABLE :: out
        out = pad_right(name, 30) // ' | ' // pad_right(real_str(ns, 4), 8) // ' | ' // &
              pad_right(real_str(dn, 3), 10) // ' | ' // pad_right(int_str(mass), 10) // ' | ' // &
              pad_right(real_str(kd, 2), 8) // ' | ' // TRIM(note)
    END FUNCTION fmt_formula_row

    FUNCTION fmt_detail_line(label, val) RESULT(out)
        CHARACTER(LEN=*), INTENT(IN) :: label, val
        CHARACTER(LEN=:), ALLOCATABLE :: out
        out = '   ' // pad_right(label, 36) // ': ' // TRIM(val)
    END FUNCTION fmt_detail_line

    SUBROUTINE emit_line(unit_no, line)
        INTEGER, INTENT(IN) :: unit_no
        CHARACTER(LEN=*), INTENT(IN) :: line
        WRITE(unit_no, '(A)') TRIM(line)
        WRITE(*, '(A)') TRIM(line)
    END SUBROUTINE emit_line

    FUNCTION grading_family(name) RESULT(out)
        CHARACTER(LEN=*), INTENT(IN) :: name
        CHARACTER(LEN=8) :: out
        CHARACTER(LEN=64) :: local
        local = ADJUSTL(name)
        out = 'UNKNOWN'
        IF (INDEX(local, 'HMA ') == 1) THEN
            out = 'HMA'
        ELSEIF (INDEX(local, 'LMA ') == 1) THEN
            out = 'LMA'
        ELSEIF (INDEX(local, 'CP') == 1) THEN
            out = 'CP'
        END IF
    END FUNCTION grading_family

    FUNCTION select_custom_family(target_mass) RESULT(out)
        REAL(dp), INTENT(IN) :: target_mass
        CHARACTER(LEN=8) :: out
        REAL(dp) :: safe_mass, best_score, score
        CHARACTER(LEN=8) :: fam
        INTEGER :: i

        safe_mass = MAX(target_mass, 1.0E-9_dp)
        best_score = HUGE(1.0_dp)
        out = 'LMA'

        DO i = 1, SIZE(standard_gradings)
            fam = grading_family(standard_gradings(i)%name)
            IF (TRIM(fam) /= 'UNKNOWN') THEN
                score = ABS(LOG(safe_mass) - LOG(MAX(standard_gradings(i)%M50, 1.0E-9_dp)))
                IF (TRIM(fam) == 'CP') score = score + 0.08_dp
                IF (score < best_score) THEN
                    best_score = score
                    out = fam
                END IF
            END IF
        END DO
    END FUNCTION select_custom_family

    FUNCTION interpolate_family_ratio(target_mass, family, note) RESULT(ratio)
        REAL(dp), INTENT(IN) :: target_mass
        CHARACTER(LEN=*), INTENT(IN) :: family
        CHARACTER(LEN=*), INTENT(OUT) :: note
        REAL(dp) :: ratio

        INTEGER, PARAMETER :: max_family = 32
        REAL(dp) :: x(max_family), y(max_family), xt, t, tmpx, tmpy
        CHARACTER(LEN=64) :: names(max_family), tmpname
        CHARACTER(LEN=8) :: fam
        INTEGER :: i, j, n

        note = ''
        n = 0
        DO i = 1, SIZE(standard_gradings)
            fam = grading_family(standard_gradings(i)%name)
            IF (TRIM(fam) == TRIM(family)) THEN
                n = n + 1
                x(n) = LOG(MAX(standard_gradings(i)%M50, 1.0E-9_dp))
                y(n) = LOG(MAX(standard_gradings(i)%nul_kg / standard_gradings(i)%nll_kg, 1.0_dp + 1.0E-9_dp))
                names(n) = standard_gradings(i)%name
            END IF
        END DO

        IF (n <= 0) THEN
            ratio = 3.0_dp
            note = 'fallback ratio R=NUL/NLL = 3.0 (no family data)'
            RETURN
        END IF

        DO i = 1, n - 1
            DO j = i + 1, n
                IF (x(i) > x(j)) THEN
                    tmpx = x(i); x(i) = x(j); x(j) = tmpx
                    tmpy = y(i); y(i) = y(j); y(j) = tmpy
                    tmpname = names(i); names(i) = names(j); names(j) = tmpname
                END IF
            END DO
        END DO

        xt = LOG(MAX(target_mass, 1.0E-9_dp))

        IF (n == 1) THEN
            ratio = EXP(y(1))
            note = TRIM(family) // ' family single-point ratio used'
            RETURN
        END IF

        IF (xt <= x(1)) THEN
            ratio = EXP(y(1))
            note = TRIM(family) // ' family ratio clamped to lower-end class ' // TRIM(names(1))
            RETURN
        END IF

        IF (xt >= x(n)) THEN
            ratio = EXP(y(n))
            note = TRIM(family) // ' family ratio clamped to upper-end class ' // TRIM(names(n))
            RETURN
        END IF

        DO i = 1, n - 1
            IF (xt >= x(i) .AND. xt <= x(i + 1)) THEN
                t = (xt - x(i)) / (x(i + 1) - x(i))
                ratio = EXP(y(i) + t * (y(i + 1) - y(i)))
                note = TRIM(family) // ' family ratio interpolated between ' // TRIM(names(i)) // ' and ' // TRIM(names(i + 1))
                RETURN
            END IF
        END DO

        ratio = EXP(y(n))
        note = TRIM(family) // ' family ratio fallback to upper-end class ' // TRIM(names(n))
    END FUNCTION interpolate_family_ratio

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

        REAL(dp) :: gamma_r, diff, min_diff, porosity
        REAL(dp) :: final_w_min, final_w_max, final_M50
        REAL(dp) :: ratio_nul_nll, nll_kg, nul_kg
        CHARACTER(LEN=64) :: selected_name
        CHARACTER(LEN=8) :: family
        CHARACTER(LEN=160) :: ratio_note
        LOGICAL :: found
        INTEGER :: i

        IF (is_armor) THEN
            ld%layer_name = 'Primary Armor'
        ELSE
            ld%layer_name = 'Underlayer'
        END IF

        ld%target_M50_kg = target_mass
        ld%target_Dn_m = target_dn
        ld%target_W_kN = target_mass * g / 1000.0_dp
        ld%used_custom_interpolation = .FALSE.
        ld%custom_family = ''
        ld%custom_ratio_nul_nll = 0.0_dp
        ld%custom_ratio_note = ''

        gamma_r = in%rho_r * g / 1000.0_dp

        ld%design_valid = .FALSE.
        IF (in%use_en13383) THEN
            selected_name = ''
            found = .FALSE.
            min_diff = HUGE(1.0_dp)

            DO i = 1, SIZE(standard_gradings)
                IF (target_mass > standard_gradings(i)%nll_kg .AND. &
                    target_mass < standard_gradings(i)%nul_kg) THEN

                    diff = standard_gradings(i)%nul_kg - standard_gradings(i)%nll_kg
                    IF (diff < min_diff) THEN
                        min_diff = diff
                        selected_name = standard_gradings(i)%name
                        final_M50 = standard_gradings(i)%M50
                        final_w_min = standard_gradings(i)%nll_kg * g / 1000.0_dp
                        final_w_max = standard_gradings(i)%nul_kg * g / 1000.0_dp
                        found = .TRUE.
                    END IF
                END IF
            END DO

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
                ld%grading_name = 'No Standard Fit (Target outside NLL-NUL)'
                ld%design_valid = .FALSE.
            END IF
        END IF

        IF (.NOT. in%use_en13383 .OR. .NOT. ld%design_valid) THEN
            IF (is_armor) THEN
                ld%grading_name = 'Custom Interpolated Grading'
            ELSE
                ld%grading_name = 'Custom Interpolated Grading Underlayer'
            END IF

            family = in%custom_family
            CALL to_upper(family)
            family = ADJUSTL(family)
            IF (TRIM(family) /= 'HMA' .AND. TRIM(family) /= 'LMA' .AND. TRIM(family) /= 'CP') THEN
                family = select_custom_family(ld%target_M50_kg)
            END IF

            ratio_nul_nll = interpolate_family_ratio(ld%target_M50_kg, TRIM(family), ratio_note)
            ratio_nul_nll = MAX(ratio_nul_nll, 1.01_dp)

            nll_kg = (2.0_dp * target_mass) / (1.0_dp + ratio_nul_nll)
            nul_kg = ratio_nul_nll * nll_kg

            ld%used_custom_interpolation = .TRUE.
            ld%custom_family = family
            ld%custom_ratio_nul_nll = ratio_nul_nll
            ld%custom_ratio_note = ratio_note

            ld%w_min_kg = nll_kg
            ld%w_max_kg = nul_kg
            ld%w_min_kn = nll_kg * g / 1000.0_dp
            ld%w_max_kn = nul_kg * g / 1000.0_dp
            ld%w_mean_kn = ld%target_W_kN
            ld%m_mean_kg = target_mass
            ld%actual_dn = target_dn
            ld%design_valid = .TRUE.
        END IF

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
        INTEGER :: i, j, iostat, n_valid, choice
        INTEGER :: valid_idx(8)
        CHARACTER(LEN=100) :: str_buf, str_buf2
        CHARACTER(LEN=64) :: uname
        
        report%inputs = in
        report%manual_override_message = ""
        
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
            CALL add_justification(report, "   In this regime, the wave height distribution strictly follows &
                &the **Rayleigh distribution**.")
            CALL add_justification(report, "   Key characteristics:")
            CALL add_justification(report, "     * The ratio H2%/Hs is constant at approximately 1.4.")
            CALL add_justification(report, "     * Wave breaking is limited to whitecapping or direct interaction &
                &with the armor layer.")
            CALL add_justification(report, "     * The spectral shape is standard (JONSWAP/Pierson-Moskowitz), &
                &and energy transfer to low-frequencies is minimal.")
            CALL add_justification(report, "")
            CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
            
            ! Van der Meer
            CALL add_justification(report, "   **A. Van der Meer (2021 Rewritten) [RECOMMENDED]**")
            CALL add_justification(report, "      * **Advantages:** This formula is the modernized industry standard.")
            CALL add_justification(report, "        Van der Meer (2021) rewrote the original formula to use the &
                &spectral period (Tm-1,0),")
            
            CALL add_justification(report, "        eliminating the influence of spectral shape.")
            CALL add_justification(report, "        Van der Meer et al. (2024) confirmed its validity for h/Hm0 > 1.5, &
                &preferring Hm0 over H1/3 for nonlinear waves.")
            CALL add_justification(report, "      * **Physics:** It correctly assumes a Rayleigh distribution of wave heights, &
                &aligning with the actual deep-water statistics.")

            ! Van Gent Mod
            CALL add_justification(report, "   **B. Van Gent Modified (2003)**")
            CALL add_justification(report, "      * **Context:** This formula incorporates the ratio H2%/Hs.")
            CALL add_justification(report, "        In deep water, with H2%/Hs = 1.4, this formula essentially &
                &converges closely with the Van der Meer predictions.")
            CALL add_justification(report, "        However, its specific calibration was focused on the effects &
                &of shallow foreshores.")

            ! Etemad-Shahidi
            CALL add_justification(report, "   **C. Etemad-Shahidi (2020)**")
            CALL add_justification(report, "      * **Comparison:** Etemad-Shahidi (2020) provides a robust formula &
                &validated for both deep and shallow water.")
            CALL add_justification(report, "        It introduces a physical permeability parameter (D_core/D_armor) &
                &to replace the nominal P factor,")
            CALL add_justification(report, "        reducing uncertainty. However, Van der Meer remains the primary &
                &standard for deep water.")

            IF (es%valid .AND. vdm%valid .AND. vdm%Dn50 > 0.0_dp) THEN
                IF (ABS(vdm%Dn50 - es%Dn50) / vdm%Dn50 > 0.10_dp) THEN
                    CALL add_justification(report, "      * **Note on Divergence:** The result deviates from Van der Meer &
                        &here. This typically occurs in the 'transition zone'")
                    CALL add_justification(report, "        of the surf similarity parameter (xi approx 2.0 - 4.5). &
                        &Etemad-Shahidi transitions to Surging physics")
                    CALL add_justification(report, "        earlier (xi > 1.8), predicting higher stability, &
                        &whereas Van der Meer maintains Plunging physics")
                    CALL add_justification(report, "        (lower stability) until a higher critical threshold. &
                        &Van der Meer is more conservative here.")
                ELSE
                    CALL add_justification(report, "        It typically converges with Van der Meer here.")
                END IF
            END IF

            CALL add_justification(report, "")
            CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
            CALL add_justification(report, "   **Use [Van der Meer (2021 Rewritten)]**.")
            CALL add_justification(report, "   It provides the most theoretically consistent result for non-depth-limited waves.")
            
            IF (es%valid .AND. vdm%valid) THEN
                 str_buf = real_str(es%Dn50, 3)
                 str_buf2 = real_str(ABS(vdm%Dn50 - es%Dn50), 3)
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
            CALL add_justification(report, "     * **Spectral Truncation:** The largest waves " // &
                "in the spectrum break on the foreshore.")
            CALL add_justification(report, "     * **Distribution Shift:** The wave height distribution deviates from &
                &Rayleigh; H2%/Hm0 drops below 1.4.")
            CALL add_justification(report, "     * **Shoaling:** Significant shoaling modifies the wave shape before impact, &
                &creating peaked crests and flat troughs.")
            CALL add_justification(report, "")
            CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
            
            ! VdM
            CALL add_justification(report, "   **A. Van der Meer (2021)**")
            CALL add_justification(report, "      * **Advantages:** Van der Meer et al. (2024) extensively re-analyzed &
                &shallow water data and concluded")
            CALL add_justification(report, "        that the rewritten Van der Meer formula (using Tm-1,0) is valid &
                &down to h/Hm0 = 1.5.")
        
            CALL add_justification(report, "        It performs reasonably well, with slightly less reliability in the &
                &1.0 < h/Hm0 < 1.5 range.")
            CALL add_justification(report, "      * **Note:** For nonlinear waves in this zone, using Hm0 is preferred &
                &over H1/3 for nonlinear waves to avoid deviations.")

            ! VG Mod
            CALL add_justification(report, "   **B. Van Gent Modified (2003)**")
            CALL add_justification(report, "      * **Constraint:** This formula explicitly relies on the ratio H2%/Hs. &
                &Research by Van der Meer et al. (2024)")
            CALL add_justification(report, "        highlights that predicting H2% accurately in this transition zone &
                &(where the ratio dips to ~1.2)")
            CALL add_justification(report, "        is notoriously inaccurate without physical modeling. The formula is valid, &
                &but the input uncertainty is high.")

            ! VG Simp
            CALL add_justification(report, "   **C. Van Gent et al. (2003) Simplified**")
            CALL add_justification(report, "      * **Context:** This formula was specifically derived for shallow foreshores.")
            CALL add_justification(report, "        However, Van der Meer et al. (2024) found that the simplified formula &
                &often does not match the data")
            CALL add_justification(report, "        in the surging domain as well as the rewritten Van der Meer formula.")
            
            CALL add_justification(report, "")
            CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
            CALL add_justification(report, "   **Use [Van der Meer (2021 Rewritten)]**.")
            CALL add_justification(report, "   Recent research (2024) confirms its validity in this depth range (h/Hm0 > 1.5), &
                &favoring it over simplified methods")
            CALL add_justification(report, "   due to the uncertainties in predicting H2% " // &
                "required for the Van Gent Modified formula.")

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
                CALL add_justification(report, "     * **Saturation:** Wave height is depth-limited (H ~ 0.5h). &
                    &Increasing offshore energy does not increase load.")
                CALL add_justification(report, "     * **Infragravity Dominance:** Scaravaglione et al. (2025) and VdM (2024) &
                    &note that infragravity waves")
                CALL add_justification(report, "       begin to dominate the spectrum, causing " // &
                    "Tm-1,0 to increase massively (up to 4x).")
                CALL add_justification(report, "     * **Formula Deviation:** Standard formulas fail here because the stability &
                    &curves flatten out (Horizontal Trend).")
                CALL add_justification(report, "")
                CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
               
                CALL add_justification(report, "   **A. Standard Formulas (VdM, VG, Standard ES)**")
                CALL add_justification(report, "      * **Failure Mode:** Scaravaglione et al. (2025) demonstrated that these &
                    &formulas fail to converge here.")
                CALL add_justification(report, "        Using the inflated spectral period results in over-predicted stability &
                    &(for surging) or under-predicted (for plunging).")
                
                CALL add_justification(report, "   **B. Scaravaglione (Modified ES 2025) [RECOMMENDED]**")
                CALL add_justification(report, "      * **Advantages:** This formula explicitly **decouples** the wave steepness &
                    &term (s^0.05) from the structure slope term.")
                CALL add_justification(report, "      * **Physics:** It is calibrated specifically for surging/bore-like waves &
                    &in the surf zone using new coefficients (c_ES,new=3.55).")
                CALL add_justification(report, "      * **Safety Cap Applied:** The formula has a weak dependence on steepness. &
                    &For very long swells, it may predict")
                CALL add_justification(report, "        theoretically high stability (Kd > 10). This system has capped Kd at 5.0 &
                    &to ensure physical realism.")
                
                CALL add_justification(report, "")
                CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
                CALL add_justification(report, "   **Use [Scaravaglione (Modified ES 2025)]**.")
                CALL add_justification(report, "   This represents the state-of-the-art for broken waves in very shallow water, &
                    &correcting the overestimation of damage.")
            ELSE
                report%recommended = vg_simp
                CALL add_justification(report, "### 1. HYDRAULIC CONTEXT: Very Shallow Water (Surf Zone) (Plunging)")
                CALL add_justification(report, "   The structure is in the surf zone, but conditions are **Plunging** (xi < 1.8).")
                CALL add_justification(report, "   The Modified ES formula is only calibrated for surging/bore conditions.")
               
                CALL add_justification(report, "")
                CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
                CALL add_justification(report, "   **Use [Van Gent Simplified (2003)]**.")
                CALL add_justification(report, "   It acts as a robust fallback. Caution is " // &
                    "advised as damage may be underpredicted &
                    &for impermeable structures.")
            END IF
        ELSE
            ! ZONE 4: Swash
            report%recommended = mod_vg
            WRITE(str_buf, '(F0.2)') report%hydro%rel_depth
            CALL add_justification(report, "### 1. HYDRAULIC CONTEXT: Extremely Shallow Water (Swash Zone)")
    
            CALL add_justification(report, "   The structure is located effectively in the **Swash Zone** (h/Hm0 = " // &
                TRIM(str_buf) // " <= 0.5).")
            CALL add_justification(report, "   Key characteristics:")
            CALL add_justification(report, "     * **Aeration:** High air entrainment reduces the effective fluid density &
                &and buoyancy of the rocks.")
            CALL add_justification(report, "     * **Impact:** Wave impact is characterized by a high-velocity turbulent bore.")
            CALL add_justification(report, "     * **Hydrostatics:** The hydrostatic cushioning effect is negligible.")
            CALL add_justification(report, "")
            CALL add_justification(report, "### 2. FORMULA COMPARISON & ANALYSIS")
            
            CALL add_justification(report, "   **A. Standard Van Gent (2003)**")
            CALL add_justification(report, "      * **Disadvantages:** Using the standard coefficient (1.75) here is **Unsafe**.")
            CALL add_justification(report, "        Scaravaglione et al. (2025) showed that stability is significantly lower than &
                &predicted by intermediate-depth formulas")
            CALL add_justification(report, "        due to the lack of buoyancy and intense turbulence.")
            
            CALL add_justification(report, "   **B. Scaravaglione (Modified VG 2025) [RECOMMENDED]**")
            CALL add_justification(report, "      * **Advantages:** This formula uses a " // &
                "recalibrated coefficient (C_VG = 3.3 instead of 1.75).")
            CALL add_justification(report, "      * **Physics:** It explicitly accounts for the increased instability in the swash &
                &zone, correcting the underestimation")
            CALL add_justification(report, "        of damage by the original VG formula in this specific regime.")
            
            CALL add_justification(report, "")
            CALL add_justification(report, "### 3. FINAL JUSTIFICATION")
            CALL add_justification(report, "   **Use [Scaravaglione (Modified VG 2025)]**.")
            CALL add_justification(report, "   It provides the necessary safety margin for " // &
                "swash zone instability where standard formulas fail.")
        END IF

        ! Interactive Manual Selection Logic
        IF (in%manual_selection) THEN
            PRINT *, ""
            PRINT *, REPEAT("=", 80)
            PRINT *, "   MANUAL SELECTION MODE"
            PRINT *, REPEAT("=", 80)

            n_valid = 0
            DO i = 1, SIZE(report%comparison)
                IF (report%comparison(i)%valid) THEN
                    n_valid = n_valid + 1
                    valid_idx(n_valid) = i
                    WRITE(*, '(A, I0, A, A)') "   [", n_valid, "] ", TRIM(report%comparison(i)%name)
                END IF
            END DO

            DO
                PRINT *, ""
                WRITE(*, '(A, I0, A)', ADVANCE='NO') "Enter the number of your preferred formula (1-", n_valid, "): "
                READ(*, '(A)') str_buf
                READ(str_buf, *, IOSTAT=iostat) choice
                IF (iostat == 0 .AND. choice >= 1 .AND. choice <= n_valid) THEN
                    report%recommended = report%comparison(valid_idx(choice))
                    uname = report%recommended%name
                    CALL to_upper(uname)
                    report%manual_override_message = "[MANUAL OVERRIDE] User switched selection to: " // TRIM(uname)
                    EXIT
                ELSE
                    PRINT *, "Invalid selection. Try again."
                END IF
            END DO
        END IF

        ! Calculate Layers
        report%armor_layer = design_layer(in%rho_r * (report%recommended%Dn50**3), report%recommended%Dn50, .TRUE., in)
        report%underlayer = design_layer(report%armor_layer%m_mean_kg / 10.0_dp, &
            (report%armor_layer%m_mean_kg / 10.0_dp / in%rho_r)**(1.0_dp/3.0_dp), .FALSE., in)

    END SUBROUTINE solve

    SUBROUTINE generate_report_file(report, filepath)
        TYPE(FullReport), INTENT(IN) :: report
        CHARACTER(LEN=*), INTENT(IN) :: filepath

        INTEGER :: u, iostatus_val, i, mass_i
        CHARACTER(LEN=:), ALLOCATABLE :: line, detail_val
        CHARACTER(LEN=115) :: separator_long
        CHARACTER(LEN=95) :: separator, eqline
        TYPE(FormulaResult) :: res
        REAL(dp) :: mass, nll_val, nul_val, ell_val, eul_val, rep_m50

        separator = REPEAT('-', 95)
        separator_long = REPEAT('-', 115)
        eqline = REPEAT('=', 95)

        OPEN(NEWUNIT=u, FILE=filepath, STATUS='REPLACE', ACTION='WRITE', IOSTAT=iostatus_val)
        IF (iostatus_val /= 0) THEN
            PRINT *, 'Error writing file.'
            RETURN
        END IF

        CALL emit_line(u, 'ROCK SLOPE STABILITY CALCULATOR')
        CALL emit_line(u, '')
        CALL emit_line(u, eqline)
        CALL emit_line(u, '   1. DESIGN INPUT PARAMETERS')
        CALL emit_line(u, eqline)
        CALL emit_line(u, fmt_row35('PARAMETER', 'VALUE', 'UNIT'))
        CALL emit_line(u, separator)
        CALL emit_line(u, fmt_row35('Significant Wave Height (Hs)', real_str(report%inputs%Hs, 2), 'm'))
        CALL emit_line(u, fmt_row35('Spectral Period (Tm-1,0)', real_str(report%inputs%Tm10, 2), 's'))
        CALL emit_line(u, fmt_row35('Water Depth (h_toe)', real_str(report%inputs%h_toe, 2), 'm'))
        CALL emit_line(u, fmt_row35('Structure Slope', '1:' // real_str(report%inputs%slope_m, 2), '(V:H)'))
        CALL emit_line(u, fmt_row35('Foreshore Slope', '1:' // real_str(report%inputs%foreshore_m, 2), '(V:H)'))
        CALL emit_line(u, fmt_row35('Permeability (P_notional)', real_str(report%inputs%P, 2), '(-)'))
        CALL emit_line(u, fmt_row35('Physical Permeability (Cp)', real_str(report%derived%Cp, 2), '(-)'))
        CALL emit_line(u, fmt_row35('Damage Level (S)', real_str(report%inputs%Sd, 1), '(-)'))
        CALL emit_line(u, separator)
        CALL emit_line(u, '')

        CALL emit_line(u, eqline)
        CALL emit_line(u, '   2. CALCULATED HYDRAULIC PARAMETERS')
        CALL emit_line(u, eqline)
        CALL emit_line(u, fmt_row40('PARAMETER', 'VALUE', 'UNIT'))
        CALL emit_line(u, separator)
        CALL emit_line(u, fmt_row40('Deep Water Wavelength (L0)', real_str(report%hydro%L0, 2), 'm'))
        CALL emit_line(u, fmt_row40('Wavelength at Toe (L_toe)', real_str(report%hydro%L_toe, 2), 'm'))
        CALL emit_line(u, fmt_row40('Wave Celerity at Toe (C)', real_str(report%hydro%C, 2), 'm/s'))
        CALL emit_line(u, fmt_row40('Group Celerity at Toe (Cg)', real_str(report%hydro%Cg, 2), 'm/s'))
        CALL emit_line(u, fmt_row40('Deep Water Steepness (s_m-1,0)', real_str(report%hydro%s_m10, 5), '(-)'))
        CALL emit_line(u, fmt_row40('Local Wave Steepness (s_local)', real_str(report%hydro%s_local, 5), '(-)'))
        CALL emit_line(u, fmt_row40('Surf Similarity (xi_m-1,0)', real_str(report%hydro%xi_m10, 2), '(-)'))
        line = pad_right('Breaker Type (Visual/Physical)', 40) // ' | ' // TRIM(report%hydro%breaker_type) // ' | (-)'
        CALL emit_line(u, line)
        CALL emit_line(u, fmt_row40('Relative Depth (h/Hm0)', real_str(report%hydro%rel_depth, 2), '(-)'))
        CALL emit_line(u, fmt_row40('Relative Depth (h/L0)', real_str(report%hydro%rel_depth * report%hydro%s_m10, 3), '(-)'))
        line = pad_right('Hydraulic Zone', 40) // ' | ' // TRIM(report%hydro%zone_desc)
        CALL emit_line(u, line)
        CALL emit_line(u, separator)
        CALL emit_line(u, '')

        CALL emit_line(u, eqline)
        CALL emit_line(u, '   3. FORMULA SELECTION & JUSTIFICATION')
        CALL emit_line(u, eqline)
        CALL emit_line(u, 'COMPARISON OF RESULTS:')
        line = 'Method                         | Ns (-)   | Dn50 (m)   | M50 (kg)   | Kd_eq (-) | NOTES'
        CALL emit_line(u, line)
        CALL emit_line(u, separator_long)
        DO i = 1, SIZE(report%comparison)
            res = report%comparison(i)
            IF (res%valid) THEN
                mass = report%inputs%rho_r * (res%Dn50**3)
                mass_i = NINT(mass)
                CALL emit_line(u, fmt_formula_row(TRIM(res%name), res%Ns, res%Dn50, mass_i, res%Kd, TRIM(res%note)))
            END IF
        END DO
        CALL emit_line(u, '')
        CALL emit_line(u, 'JUSTIFICATION & ANALYSIS:')
        IF (ALLOCATED(report%justification)) THEN
            DO i = 1, SIZE(report%justification)
                CALL emit_line(u, TRIM(report%justification(i)))
            END DO
        END IF
        CALL emit_line(u, separator)

        IF (LEN_TRIM(report%manual_override_message) > 0) THEN
            CALL emit_line(u, '')
            CALL emit_line(u, TRIM(report%manual_override_message))
        END IF
        CALL emit_line(u, '')

        CALL emit_line(u, eqline)
        IF (report%inputs%use_en13383) THEN
            CALL emit_line(u, '   4. ROCK ARMOUR LAYER DESIGN (EN 13383 Standard)')
        ELSE
            CALL emit_line(u, '   4. ROCK ARMOUR LAYER DESIGN (Custom Grading)')
        END IF
        CALL emit_line(u, eqline)
        CALL emit_line(u, 'PRIMARY ARMOR LAYER')
        CALL emit_line(u, '   Theoretical Required W    : ' // real_str(report%armor_layer%target_W_kN, 2) // ' kN')
        CALL emit_line(u, '   Theoretical Required M50  : ' // int_str(NINT(report%armor_layer%target_M50_kg)) // ' kg')
        CALL emit_line(u, '   Theoretical Required Dn50 : ' // real_str(report%armor_layer%target_Dn_m, 3) // ' m')
        CALL emit_line(u, '----------------------------------------')

        IF (.NOT. report%armor_layer%design_valid .AND. report%inputs%use_en13383) THEN
            CALL emit_line(u, '   [WARNING] No standard EN13383 grading found for this mass.')
        ELSE
            CALL emit_line(u, fmt_detail_line('Adopted rock grading', TRIM(report%armor_layer%grading_name)))
            IF (report%inputs%use_en13383) THEN
                nll_val = report%armor_layer%w_min_kg
                nul_val = report%armor_layer%w_max_kg
                ell_val = 0.7_dp * nll_val
                eul_val = 1.5_dp * nul_val
                rep_m50 = 0.5_dp * (nll_val + nul_val)
                CALL emit_line(u, fmt_detail_line('Representative M50', real_str(rep_m50, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal lower limit (NLL)', real_str(nll_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal upper limit (NUL)', real_str(nul_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme lower limit (ELL)', real_str(ell_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme upper limit (EUL)', real_str(eul_val, 1) // ' kg'))
            ELSE
                nll_val = report%armor_layer%w_min_kg
                nul_val = report%armor_layer%w_max_kg
                ell_val = 0.7_dp * nll_val
                eul_val = 1.5_dp * nul_val
                CALL emit_line(u, fmt_detail_line('Custom family basis', TRIM(report%armor_layer%custom_family)))
                CALL emit_line(u, fmt_detail_line('Interpolated ratio (NUL/NLL)', &
                    real_str(report%armor_layer%custom_ratio_nul_nll, 3)))
                CALL emit_line(u, fmt_detail_line('Interpolation rule', TRIM(report%armor_layer%custom_ratio_note)))
                CALL emit_line(u, fmt_detail_line('Representative M50', real_str(report%armor_layer%m_mean_kg, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal lower limit (NLL)', real_str(nll_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal upper limit (NUL)', real_str(nul_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme lower limit (ELL)', real_str(ell_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme upper limit (EUL)', real_str(eul_val, 1) // ' kg'))
            END IF
            CALL emit_line(u, fmt_detail_line('Nominal Diameter (Dn_rock)', real_str(report%armor_layer%actual_dn, 3) // ' m'))
            CALL emit_line(u, fmt_detail_line('Double Layer Thickness', real_str(report%armor_layer%thickness, 2) // ' m'))
            CALL emit_line(u, fmt_detail_line('Packing Density [rocks/100m2]', real_str(report%armor_layer%packing_density, 2)))
        END IF

        CALL emit_line(u, separator)
        CALL emit_line(u, 'UNDERLAYER (FILTER LAYER)')
        CALL emit_line(u, '   Target Weight (M50 / 10)  : ' // real_str(report%underlayer%target_W_kN, 3) // ' kN')
        CALL emit_line(u, '   Target Mass (M50 / 10)    : ' // real_str(report%underlayer%target_M50_kg, 1) // ' kg')
        CALL emit_line(u, '----------------------------------------')

        IF (.NOT. report%underlayer%design_valid .AND. report%inputs%use_en13383) THEN
            CALL emit_line(u, '   [WARNING] No suitable standard underlayer grading found.')
        ELSE
            CALL emit_line(u, fmt_detail_line('Adopted rock grading', TRIM(report%underlayer%grading_name)))
            IF (report%inputs%use_en13383) THEN
                nll_val = report%underlayer%w_min_kg
                nul_val = report%underlayer%w_max_kg
                ell_val = 0.7_dp * nll_val
                eul_val = 1.5_dp * nul_val
                rep_m50 = 0.5_dp * (nll_val + nul_val)
                CALL emit_line(u, fmt_detail_line('Representative M50', real_str(rep_m50, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal lower limit (NLL)', real_str(nll_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal upper limit (NUL)', real_str(nul_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme lower limit (ELL)', real_str(ell_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme upper limit (EUL)', real_str(eul_val, 1) // ' kg'))
            ELSE
                nll_val = report%underlayer%w_min_kg
                nul_val = report%underlayer%w_max_kg
                ell_val = 0.7_dp * nll_val
                eul_val = 1.5_dp * nul_val
                CALL emit_line(u, fmt_detail_line('Custom family basis', TRIM(report%underlayer%custom_family)))
                CALL emit_line(u, fmt_detail_line('Interpolated ratio (NUL/NLL)', &
                    real_str(report%underlayer%custom_ratio_nul_nll, 3)))
                CALL emit_line(u, fmt_detail_line('Interpolation rule', TRIM(report%underlayer%custom_ratio_note)))
                CALL emit_line(u, fmt_detail_line('Representative M50', real_str(report%underlayer%m_mean_kg, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal lower limit (NLL)', real_str(nll_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Nominal upper limit (NUL)', real_str(nul_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme lower limit (ELL)', real_str(ell_val, 1) // ' kg'))
                CALL emit_line(u, fmt_detail_line('Extreme upper limit (EUL)', real_str(eul_val, 1) // ' kg'))
            END IF
            CALL emit_line(u, fmt_detail_line('Nominal Diameter (Dn_rock)', real_str(report%underlayer%actual_dn, 3) // ' m'))
            CALL emit_line(u, fmt_detail_line('Double Layer Thickness', real_str(report%underlayer%thickness, 2) // ' m'))
            CALL emit_line(u, fmt_detail_line('Packing Density [rocks/100m2]', real_str(report%underlayer%packing_density, 2)))
        END IF
        CALL emit_line(u, eqline)
        CLOSE(u)
        WRITE(*, '(A)') ''
        WRITE(*, '(A)') '[System] Report saved to ' // TRIM(filepath)
    END SUBROUTINE generate_report_file


END PROGRAM RockSlopeCalculator