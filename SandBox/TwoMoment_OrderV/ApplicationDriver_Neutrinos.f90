PROGRAM ApplicationDriver_Neutrinos

  USE KindModule, ONLY: &
    DP, Zero, One, Two, SqrtTiny
  USE UnitsModule, ONLY: &
    Kilometer, &
    Millisecond, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF, uDF
  USE Euler_UtilitiesModule_NonRelativistic, ONLY: &
    ComputeFromConserved_Euler_NonRelativistic
  USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplySlopeLimiter_Euler_NonRelativistic_TABLE
  USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
    ApplyPositivityLimiter_Euler_NonRelativistic_TABLE
  USE RadiationFieldsModule, ONLY: &
    uCR, uPR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE TwoMoment_UtilitiesModule_OrderV, ONLY: &
    ComputeFromConserved_TwoMoment
  USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
    ApplySlopeLimiter_TwoMoment
  USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
    ApplyPositivityLimiter_TwoMoment
  USE TwoMoment_DiscretizationModule_Collisions_Neutrinos_OrderV, ONLY: &
    ComputeIncrement_TwoMoment_Implicit
  USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
    Update_IMEX_RK
  USE InitializationModule_Neutrinos, ONLY: &
    InitializeFields

  IMPLICIT NONE

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: TimeSteppingScheme
  CHARACTER(32) :: CoordinateSystem
  CHARACTER(64) :: EosTableName
  CHARACTER(64) :: OpacityTableName_AbEm
  CHARACTER(64) :: OpacityTableName_Iso
  CHARACTER(64) :: OpacityTableName_NES
  CHARACTER(64) :: OpacityTableName_Pair
  LOGICAL       :: EvolveEuler
  LOGICAL       :: UseSlopeLimiter_Euler
  LOGICAL       :: UseSlopeLimiter_TwoMoment
  LOGICAL       :: UsePositivityLimiter_Euler
  LOGICAL       :: UsePositivityLimiter_TwoMoment
  LOGICAL       :: FixedTimeStep
  INTEGER       :: nSpecies
  INTEGER       :: nNodes
  INTEGER       :: nE, bcE, nX(3), bcX(3)
  INTEGER       :: iCycle, iCycleD, iCycleW
  REAL(DP)      :: xL(3), xR(3), ZoomX(3) = One
  REAL(DP)      :: eL, eR, ZoomE = One
  REAL(DP)      :: t, dt, dt_CFL, dt_FXD, t_end

  ProgramName = 'Relaxation'

  CoordinateSystem = 'CARTESIAN'

  EosTableName          = 'wl-EOS-SFHo-15-25-50.h5'
  OpacityTableName_AbEm = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5'
  OpacityTableName_Iso  = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5'
  OpacityTableName_NES  = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5'
  OpacityTableName_Pair = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5'

  FixedTimeStep = .FALSE.

  SELECT CASE( TRIM( ProgramName ) )

    CASE( 'Relaxation' )

      nSpecies = 2
      nNodes   = 2

      nX  = [ 1, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      bcX = [ 0, 0, 0 ]

      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 0
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'BackwardEuler'

      t_end = 1.0d0 * Millisecond

      FixedTimeStep = .TRUE.
      dt_FXD        = 1.0d-3 * Millisecond
      iCycleD       = 1
      iCycleW       = 10

      EvolveEuler                    = .FALSE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .FALSE.
      UsePositivityLimiter_TwoMoment = .FALSE.

    CASE( 'Deleptonization' )

      nSpecies = 2
      nNodes   = 2

      nX  = [ 64, 1, 1 ]
      xL  = [ - 50.0_DP, 0.0_DP, 0.0_DP ] * Kilometer
      xR  = [ + 50.0_DP, 1.0_DP, 1.0_DP ] * Kilometer
      bcX = [ 31, 1, 1 ]

      nE    = 16
      eL    = 0.0d0 * MeV
      eR    = 3.0d2 * MeV
      bcE   = 10
      ZoomE = 1.266038160710160_DP

      TimeSteppingScheme = 'IMEX_PDARS'

      t_end = 1.0d0 * Millisecond

      iCycleD = 1
      iCycleW = 10

      EvolveEuler                    = .TRUE.
      UseSlopeLimiter_Euler          = .FALSE.
      UseSlopeLimiter_TwoMoment      = .FALSE.
      UsePositivityLimiter_Euler     = .TRUE.
      UsePositivityLimiter_TwoMoment = .TRUE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A6,A,A)') '', 'Unknown Program Name: ', TRIM( ProgramName )
      WRITE(*,*)
      STOP

  END SELECT

  CALL InitializeDriver

  CALL InitializeFields

  ! --- Apply Slope Limiter to Initial Data ---

  CALL ApplySlopeLimiter_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplySlopeLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  ! --- Apply Positivity Limiter to Initial Data ---

  CALL ApplyPositivityLimiter_Euler_NonRelativistic_TABLE &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uDF )

  CALL ApplyPositivityLimiter_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGE, uGF, uCF, uCR )

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  ! --- Evolve ---

  t = Zero

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, ' to t = ', t_end / Millisecond
  WRITE(*,*)

  iCycle = 0
  DO WHILE( t < t_end )

    iCycle = iCycle + 1

    IF( FixedTimeStep )THEN

      dt = dt_FXD

    ELSE

      dt_CFL = 0.3_DP * MINVAL( (xR-xL)/DBLE(nX) ) / ( Two*DBLE(nNodes-1)+One )

      dt = dt_CFL

    END IF

    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN

      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A5,ES12.6E2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = '    ,  t / Millisecond, &
          '', 'dt = '   , dt / Millisecond

    END IF

    CALL Update_IMEX_RK &
           ( dt, uGE, uGF, uCF, uCR, ComputeIncrement_TwoMoment_Implicit )

    t = t + dt

    IF( MOD( iCycle, iCycleW ) == 0 )THEN

      CALL ComputeFromConserved_Euler_NonRelativistic &
             ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

      CALL ComputeFromConserved_TwoMoment &
             ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR )

      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .TRUE., &
               WriteFF_Option = .TRUE., &
               WriteRF_Option = .TRUE. )

    END IF

  END DO

  CALL ComputeFromConserved_Euler_NonRelativistic &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, uCF, uPF, uAF )

  CALL ComputeFromConserved_TwoMoment &
         ( iZ_B0, iZ_E0, iZ_B1, iZ_E1, uGF, uCF, uCR, uPR )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )

  CALL FinalizeDriver

CONTAINS


  SUBROUTINE InitializeDriver

    USE TwoMoment_TimersModule_OrderV, ONLY: &
      InitializeTimers
    USE ProgramInitializationModule, ONLY: &
      InitializeProgram
    USE ReferenceElementModuleX, ONLY: &
      InitializeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      InitializeReferenceElementX_Lagrange
    USE GeometryComputationModule, ONLY: &
      ComputeGeometryX
    USE ReferenceElementModuleE, ONLY: &
      InitializeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      InitializeReferenceElementE_Lagrange
    USE GeometryComputationModuleE, ONLY: &
      ComputeGeometryE
    USE ReferenceElementModuleZ, ONLY: &
      InitializeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      InitializeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      InitializeReferenceElement_Lagrange
    USE EquationOfStateModule_TABLE, ONLY: &
      InitializeEquationOfState_TABLE, &
      MinD, MaxD, MinT, MaxT, MinY, MaxY
    USE OpacityModule_TABLE, ONLY: &
      InitializeOpacities_TABLE
    USE TwoMoment_ClosureModule, ONLY: &
      InitializeClosure_TwoMoment
    USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
      InitializeSlopeLimiter_Euler_NonRelativistic_TABLE
    USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
      InitializePositivityLimiter_Euler_NonRelativistic_TABLE
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      InitializeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
      InitializeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
      InitializePositivityLimiter_TwoMoment
    USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
      Initialize_IMEX_RK

    CALL InitializeTimers

    CALL InitializeProgram &
           ( ProgramName_Option &
               = TRIM( ProgramName ), &
             nX_Option &
               = nX, &
             swX_Option &
               = [ 1, 0, 0 ], &
             bcX_Option &
               = bcX, &
             xL_Option &
               = xL, &
             xR_Option &
               = xR, &
             zoomX_Option &
               = zoomX, &
             nE_Option &
               = nE, &
             swE_Option &
               = 1, &
             bcE_Option &
               = bcE, &
             eL_Option &
               = eL, &
             eR_Option &
               = eR, &
             zoomE_Option &
               = zoomE, &
             nNodes_Option &
               = nNodes, &
             CoordinateSystem_Option &
               = TRIM( CoordinateSystem ), &
             ActivateUnits_Option &
               = .TRUE., &
             nSpecies_Option &
               = nSpecies, &
             BasicInitialization_Option &
               = .TRUE. )

    ! --- Position Space Reference Element and Geometry ---

    CALL InitializeReferenceElementX

    CALL InitializeReferenceElementX_Lagrange

    CALL ComputeGeometryX &
           ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

    ! --- Energy Space Reference Element and Geometry ---

    CALL InitializeReferenceElementE

    CALL InitializeReferenceElementE_Lagrange

    CALL ComputeGeometryE &
           ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

    ! --- Phase Space Reference Element ---

    CALL InitializeReferenceElementZ

    CALL InitializeReferenceElement

    CALL InitializeReferenceElement_Lagrange

    ! --- Initialize Equation of State ---

    CALL InitializeEquationOfState_TABLE &
           ( EquationOfStateTableName_Option &
               = EosTableName, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Opacities ---

    CALL InitializeOpacities_TABLE &
           ( OpacityTableName_EmAb_Option &
               = TRIM( OpacityTableName_AbEm ), &
             OpacityTableName_Iso_Option  &
               = TRIM( OpacityTableName_Iso ), &
             OpacityTableName_NES_Option &
               = TRIM( OpacityTableName_NES ), &
             OpacityTableName_Pair_Option &
               = TRIM( OpacityTableName_Pair ), &
             EquationOfStateTableName_Option &
               = TRIM( EosTableName ), &
             Verbose_Option = .TRUE. )

    ! --- Initialize Moment Closure ---

    CALL InitializeClosure_TwoMoment

    ! --- Initialize Slope Limiter (Euler) ---

    CALL InitializeSlopeLimiter_Euler_NonRelativistic_TABLE &
           ( BetaTVD_Option &
               = 1.75_DP, &
             SlopeTolerance_Option &
               = 1.0d-6, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter_Euler, &
             UseTroubledCellIndicator_Option &
               = .FALSE., &
             LimiterThresholdParameter_Option &
               = Zero )

    ! --- Initialize Positivity Limiter (Euler) ---

    CALL InitializePositivityLimiter_Euler_NonRelativistic_TABLE &
           ( UsePositivityLimiter_Option &
               = UsePositivityLimiter_Euler, &
             Verbose_Option &
               = .TRUE., &
             Min_1_Option &
               = ( One + EPSILON( One ) ) * MinD, &
             Min_2_Option &
               = ( One + EPSILON( One ) ) * MinT, &
             Min_3_Option &
               = ( One + EPSILON( One ) ) * MinY, &
             Max_1_Option &
               = ( One - EPSILON( One ) ) * MaxD, &
             Max_2_Option &
               = ( One - EPSILON( One ) ) * MaxT, &
             Max_3_Option &
               = ( One - EPSILON( One ) ) * MaxY )

    ! --- Initialize Troubled Cell Indicator (Two-Moment) ---

    CALL InitializeTroubledCellIndicator_TwoMoment &
           ( UseTroubledCellIndicator_Option &
               = .FALSE., &
             C_TCI_Option &
               = Zero, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Slope Limiter (Two-Moment) ---

    CALL InitializeSlopeLimiter_TwoMoment &
           ( BetaTVD_Option &
               = 1.75_DP, &
             UseSlopeLimiter_Option &
               = UseSlopeLimiter_TwoMoment, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Positivity Limiter (Two-Moment) ---

    CALL InitializePositivityLimiter_TwoMoment &
           ( Min_1_Option &
               = SqrtTiny, &
             Min_2_Option &
               = SqrtTiny, &
             UsePositivityLimiter_Option &
               = UsePositivityLimiter_TwoMoment, &
             Verbose_Option &
               = .TRUE. )

    ! --- Initialize Time Stepper ---

    CALL Initialize_IMEX_RK &
           ( TRIM( TimeSteppingScheme ), EvolveEuler_Option = EvolveEuler )

  END SUBROUTINE InitializeDriver


  SUBROUTINE FinalizeDriver

    USE TwoMoment_TimeSteppingModule_OrderV, ONLY: &
      Finalize_IMEX_RK
    USE EquationOfStateModule_TABLE, ONLY: &
      FinalizeEquationOfState_TABLE
    USE OpacityModule_TABLE, ONLY: &
      FinalizeOpacities_TABLE
    USE Euler_SlopeLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE
    USE Euler_PositivityLimiterModule_NonRelativistic_TABLE, ONLY: &
      FinalizePositivityLimiter_Euler_NonRelativistic_TABLE
    USE TwoMoment_TroubledCellIndicatorModule, ONLY: &
      FinalizeTroubledCellIndicator_TwoMoment
    USE TwoMoment_SlopeLimiterModule_OrderV, ONLY: &
      FinalizeSlopeLimiter_TwoMoment
    USE TwoMoment_PositivityLimiterModule_OrderV, ONLY: &
      FinalizePositivityLimiter_TwoMoment
    USE ReferenceElementModuleX, ONLY: &
      FinalizeReferenceElementX
    USE ReferenceElementModuleX_Lagrange, ONLY: &
      FinalizeReferenceElementX_Lagrange
    USE ReferenceElementModuleE, ONLY: &
      FinalizeReferenceElementE
    USE ReferenceElementModuleE_Lagrange, ONLY: &
      FinalizeReferenceElementE_Lagrange
    USE ReferenceElementModuleZ, ONLY: &
      FinalizeReferenceElementZ
    USE ReferenceElementModule, ONLY: &
      FinalizeReferenceElement
    USE ReferenceElementModule_Lagrange, ONLY: &
      FinalizeReferenceElement_Lagrange
    USE ProgramInitializationModule, ONLY: &
      FinalizeProgram
    USE TwoMoment_TimersModule_OrderV, ONLY: &
      FinalizeTimers

    CALL Finalize_IMEX_RK

    CALL FinalizeEquationOfState_TABLE

    CALL FinalizeOpacities_TABLE

    CALL FinalizeSlopeLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizePositivityLimiter_Euler_NonRelativistic_TABLE

    CALL FinalizeTroubledCellIndicator_TwoMoment

    CALL FinalizeSlopeLimiter_TwoMoment

    CALL FinalizePositivityLimiter_TwoMoment

    CALL FinalizeReferenceElementX

    CALL FinalizeReferenceElementX_Lagrange

    CALL FinalizeReferenceElementE

    CALL FinalizeReferenceElementE_Lagrange

    CALL FinalizeReferenceElementZ

    CALL FinalizeReferenceElement

    CALL FinalizeReferenceElement_Lagrange

    CALL FinalizeProgram

    CALL FinalizeTimers

  END SUBROUTINE FinalizeDriver


END PROGRAM ApplicationDriver_Neutrinos
