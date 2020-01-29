MODULE InitializationModule_Relativistic

  USE KindModule, ONLY: &
    DP, Zero, SqrtTiny, Half, One, Two, Three, Pi, Four, TwoPi, FourPi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nX, nNodesX, &
    nDOFX, &
    iX_B0, iX_B1, iX_E0, iX_E1
  USE ReferenceElementModuleX, ONLY: &
    NodeNumberTableX
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  USE FluidFieldsModule, ONLY: &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    uAF, iAF_P
  USE EquationOfStateModule_IDEAL, ONLY: &
    Gamma_IDEAL
  USE Euler_UtilitiesModule_Relativistic, ONLY: &
    ComputeConserved_Euler_Relativistic
  USE UnitsModule, ONLY: &
    Meter, SpeedOfLight, Kilometer, SolarMass, Kilogram, Second, Joule, &
    Gram, Centimeter, Erg
  USE UtilitiesModule, ONLY: &
    Locate, &
    NodeNumberX

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeFields_Relativistic

CONTAINS


  SUBROUTINE InitializeFields_Relativistic &
               ( AdvectionProfile_Option, &
                 RiemannProblemName_Option, &
                 nDetCells_Option, Eblast_Option, &
                 MassPNS_Option, ShockRadius_Option, &
                 AccretionRate_Option, MachNumber_Option, &
                 ApplyPerturbation_Option, PerturbationOrder_Option, &
                 PerturbationAmplitude_Option, &
                 rPerturbationInner_Option, rPerturbationOuter_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: AdvectionProfile_Option
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: RiemannProblemName_Option
    INTEGER,          INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Eblast_Option
    REAL(DP),         INTENT(in), OPTIONAL :: MassPNS_Option
    REAL(DP),         INTENT(in), OPTIONAL :: ShockRadius_Option
    REAL(DP),         INTENT(in), OPTIONAL :: AccretionRate_Option
    REAL(DP),         INTENT(in), OPTIONAL :: MachNumber_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplyPerturbation_Option
    INTEGER,          INTENT(in), OPTIONAL :: PerturbationOrder_Option
    REAL(DP),         INTENT(in), OPTIONAL :: PerturbationAmplitude_Option
    REAL(DP),         INTENT(in), OPTIONAL :: rPerturbationInner_Option
    REAL(DP),         INTENT(in), OPTIONAL :: rPerturbationOuter_Option

    CHARACTER(LEN=64) :: AdvectionProfile = 'SineWave'
    CHARACTER(LEN=64) :: RiemannProblemName = 'Sod'

    ! --- Sedov-Taylor Blast Wave (Defaults) ---
    INTEGER  :: nDetCells = 1
    REAL(DP) :: Eblast    = 1.0d-3

    ! --- Standing Accretion Shock (Defaults) ---
    REAL(DP) :: MassPNS               = 1.4_DP * SolarMass
    REAL(DP) :: ShockRadius           = 180.0_DP * Kilometer
    REAL(DP) :: AccretionRate         = 0.3_DP * SolarMass / Second
    REAL(DP) :: MachNumber            = 10.0_DP
    LOGICAL  :: ApplyPerturbation     = .FALSE.
    INTEGER  :: PerturbationOrder     = 0
    REAL(DP) :: PerturbationAmplitude = 0.0_DP
    REAL(DP) :: rPerturbationInner    = 0.0_DP
    REAL(DP) :: rPerturbationOuter    = 0.0_DP

    IF( PRESENT( AdvectionProfile_Option ) ) &
      AdvectionProfile = TRIM( AdvectionProfile_Option )

    IF( PRESENT( RiemannProblemName_Option ) ) &
      RiemannProblemName = TRIM( RiemannProblemName_Option )

    IF( PRESENT( nDetCells_Option ) ) &
      nDetCells = nDetCells_Option
    IF( PRESENT( Eblast_Option ) ) &
      Eblast = Eblast_Option

    IF( PRESENT( MassPNS_Option ) ) &
      MassPNS = MassPNS_Option
    IF( PRESENT( ShockRadius_Option ) ) &
      ShockRadius = ShockRadius_Option
    IF( PRESENT( AccretionRate_Option ) ) &
      AccretionRate = AccretionRate_Option
    IF( PRESENT( MachNumber_Option ) ) &
      MachNumber = MachNumber_Option
    IF( PRESENT( ApplyPerturbation_Option ) ) &
      ApplyPerturbation = ApplyPerturbation_Option
    IF( PRESENT( PerturbationOrder_Option ) ) &
      PerturbationOrder = PerturbationOrder_Option
    IF( PRESENT( PerturbationAmplitude_Option ) ) &
      PerturbationAmplitude = PerturbationAmplitude_Option
    IF( PRESENT( rPerturbationInner_Option ) ) &
      rPerturbationInner = rPerturbationInner_Option
    IF( PRESENT( rPerturbationOuter_Option ) ) &
      rPerturbationOuter = rPerturbationOuter_Option

    WRITE(*,*)
    WRITE(*,'(A,A)') '    INFO: ', TRIM( ProgramName )

    SELECT CASE ( TRIM( ProgramName ) )

      CASE( 'Polynomial' )

        CALL InitializeFields_Polynomial

      CASE( 'Advection' )

        CALL InitializeFields_Advection &
               ( TRIM( AdvectionProfile ) )

      CASE( 'Advection2D' )

        CALL InitializeFields_Advection2D &
               ( TRIM( AdvectionProfile ) )

      CASE( 'RiemannProblem' )

        CALL InitializeFields_RiemannProblem &
               ( TRIM( RiemannProblemName ), &
                 nDetCells_Option = nDetCells, &
                 Eblast_Option    = Eblast )

      CASE( 'RiemannProblem2D' )

        CALL InitializeFields_RiemannProblem2D &
               ( TRIM( RiemannProblemName ) )

      CASE( 'RiemannProblemSpherical' )

        CALL InitializeFields_RiemannProblemSpherical &
               ( TRIM( RiemannProblemName ) )

      CASE( 'SedovTaylorBlastWave' )

        CALL InitializeFields_SedovTaylorBlastWave &
               ( nDetCells, Eblast )

      CASE( 'KelvinHelmholtzInstability' )

         CALL InitializeFields_KelvinHelmholtzInstability

      CASE( 'StandingAccretionShock' )

        CALL InitializeFields_StandingAccretionShock &
               ( MassPNS, ShockRadius, AccretionRate, MachNumber, &
                 ApplyPerturbation, PerturbationOrder, PerturbationAmplitude, &
                 rPerturbationInner, rPerturbationOuter )

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
        WRITE(*,'(A)')     'Stopping...'
        STOP

    END SELECT

  END SUBROUTINE InitializeFields_Relativistic


  SUBROUTINE InitializeFields_Polynomial

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1, a0, a1, a2, a3, a4, a5

    a0 = 1.0_DP
    a1 = 0.1_DP
    a2 = 4.0_DP
    a3 = -0.4_DP
    a4 = -1.0_DP

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        uPF(iNodeX,iX1,iX2,iX3,iPF_D ) &
          = a0*X1**0 + a1*X1**1 + a2*X1**2 &
              + a3*X1**3 + a4*X1**4
        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = One
        uPF(iNodeX,iX1,iX2,iX3,iPF_E ) = One / ( Gamma_IDEAL - One )

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Polynomial


  SUBROUTINE InitializeFields_Advection( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          CASE( 'SineWave' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE( 'TopHat' )

            IF( X1 .GT. 0.45 .AND. X1 .LT. 0.55 )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 2.0_DP

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP

            END IF

              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWave'
            WRITE(*,'(A)') '  TopHat'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection


  SUBROUTINE InitializeFields_Advection2D( AdvectionProfile )

    CHARACTER(LEN=*), INTENT(in) :: AdvectionProfile

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Advection Profile: ', TRIM( AdvectionProfile )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        SELECT CASE( TRIM( AdvectionProfile ) )

          CASE( 'SineWaveX1' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X1 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE( 'SineWaveX2' )

            uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = One + 0.1_DP * SIN( TwoPi * X2 )
            uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.1_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
            uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = 1.0_DP
            uPF(iNodeX,iX1,iX2,iX3,iPF_E )  &
              = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for AdvectionProfile: ', AdvectionProfile
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') '  SineWaveX1'
            WRITE(*,'(A)') '  SineWaveX2'
            WRITE(*,*)
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_Advection2D


  SUBROUTINE InitializeFields_RiemannProblem &
               ( RiemannProblemName, &
                 nDetCells_Option, Eblast_Option )

    CHARACTER(LEN=*), INTENT(in)           :: RiemannProblemName
    INTEGER,          INTENT(in), OPTIONAL :: nDetCells_Option
    REAL(DP),         INTENT(in), OPTIONAL :: Eblast_Option

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1
    REAL(DP) :: X1

    INTEGER  :: nDetCells
    REAL(DP) :: Eblast, X_D

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'Sod' )

            IF( X1 .LE. Half )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.125_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF

          CASE( 'MBProblem1' )

            IF( X1 .LE. Half )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.9_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 10.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF

          CASE( 'MBProblem4' )

            IF( X1 .LE. Half )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d3
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d-2
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF

          CASE( 'PerturbedShockTube' )

            IF( X1 .LE. Half )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 5.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 50.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  &
                = 2.0_DP + 0.3_DP * SIN( 50.0_DP * X1 )
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 5.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF

          CASE( 'ShockReflection' )

            IF( X1 .LE. One )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.99999_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.01_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF


          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for RiemannProblemName: ', RiemannProblemName
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') &
              "  'Sod' - &
              Sod's shock tube"
            WRITE(*,'(A)') &
              "  'MBProblem1' - &
              Mignone & Bodo (2005) MNRAS, 364, 126, Problem 1"
            WRITE(*,'(A)') &
              "  'MBProblem4' - &
              Mignone & Bodo (2005) MNRAS, 364, 126, Problem 4"
            WRITE(*,'(A)') &
              "  'PerturbedShockTube' - &
              Del Zanna & Bucciantini (2002) AA, 390, 1177, &
              Sinusoidal density perturbation"
            WRITE(*,'(A)') &
              "  'ShockReflection' - &
              Del Zanna & Bucciantini (2002) AA, 390, 1177, &
              Planar shock reflection"
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblem



  SUBROUTINE InitializeFields_RiemannProblem2D( RiemannProblemName )

    CHARACTER(LEN=*), INTENT(in) :: RiemannProblemName

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1, iNodeX2
    REAL(DP)      :: X1, X2

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', '2D Riemann Problem Name: ', TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'DzB2002' )

            ! --- SW ---
            IF( X1 .LE. Half .AND. X2 .LE. Half )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.5_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ! --- NW ---
            ELSE IF( X1 .LE. Half .AND. X2 .GT. Half )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.99_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ! --- NE ---
            ELSE IF( X1 .GT. Half .AND. X2 .GT. Half )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.01_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ! --- SE ---
            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.99_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF


          CASE DEFAULT

            WRITE(*,*)
            WRITE(*,'(A,A)') &
              'Invalid choice for RiemannProblemName: ', &
                RiemannProblemName
            WRITE(*,'(A)') 'Valid choices:'
            WRITE(*,'(A)') &
              "  'DzB2002' - &
              Del-Zanna & Bucciantini, 2D Riemann Problem"
            WRITE(*,'(A)') 'Stopping...'
            STOP

        END SELECT

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO


  END SUBROUTINE InitializeFields_RiemannProblem2D


  SUBROUTINE InitializeFields_RiemannProblemSpherical( RiemannProblemName )

    CHARACTER(LEN=*), INTENT(in) :: RiemannProblemName

    INTEGER       :: iX1, iX2, iX3
    INTEGER       :: iNodeX, iNodeX1
    REAL(DP)      :: X1

    WRITE(*,*)
    WRITE(*,'(A4,A,A)') &
      '', 'Spherical Riemann Problem Name: ', &
        TRIM( RiemannProblemName )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        SELECT CASE ( TRIM( RiemannProblemName ) )

          CASE( 'SphericalSod' )

            IF( X1 <= One )THEN

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            ELSE

              uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 0.125_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
              uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 0.1_DP
              uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
                = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

            END IF

         CASE DEFAULT

            WRITE(*,*)
            WRITE(*,*) &
              'Invalid choice for RiemannProblemName: ', &
              RiemannProblemName
            WRITE(*,*) 'Valid choices:'
            WRITE(*,*) &
              "'SphericalSod' - ", &
              "Spherical Sod's shock tube"
            STOP

          END SELECT

        END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_RiemannProblemSpherical


  SUBROUTINE InitializeFields_SedovTaylorBlastWave( nDetCells, Eblast )

    INTEGER,  INTENT(in) :: nDetCells
    REAL(DP), INTENT(in) :: Eblast

    INTEGER  :: iX1, iX2, iX3, iNodeX1, iNodeX
    REAL(DP) :: X1, X_D

    X_D = DBLE( nDetCells ) * MeshX(1) % Width(1)
    WRITE(*,*)
    WRITE(*,'(A,I4.4)')      '     nDetCells:              ', nDetCells
    WRITE(*,'(A,ES23.16E3)') '     Initial blast radius:   ', X_D
    WRITE(*,'(A,ES23.16E3)') '     Blast energy:           ', Eblast
    WRITE(*,'(A,ES23.16E3)') '     Initial blast pressure: ', &
                                     ( Gamma_IDEAL - One ) &
                                       * Eblast / ( FourPi / Three * X_D**3 )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 <= X_D)THEN

          uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
            = Eblast / ( FourPi / Three * X_D**3 )
          uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
            = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D)  = 1.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0_DP
          uPF(iNodeX,iX1,iX2,iX3,iPF_E)  &
            = 1.0d-5
          uAF(iNodeX,iX1,iX2,iX3,iAF_P)  &
            = ( Gamma_IDEAL - One ) * uPF(iNodeX,iX1,iX2,iX3,iPF_E)

        END IF

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_SedovTaylorBlastWave


  ! --- Relativistic 2D Kelvin-Helmholtz instability a la
  !     Beckwith & Stone (2011), ApjS, 193, 6 (typo in Eq. (63)) ---
  SUBROUTINE InitializeFields_KelvinHelmholtzInstability

    INTEGER  :: iX1, iX2, iX3
    INTEGER  :: iNodeX, iNodeX1, iNodeX2
    REAL(DP) :: X1, X2
    REAL(DP) :: rho0, rho1
    REAL(DP) :: Vshear, a, X2_Offset, sigma, A0

    rho0 = 0.505d0
    rho1 = 0.495d0

    Vshear    = 0.5d0
    a         = 0.01d0
    X2_Offset = 0.5d0
    sigma     = 0.1d0

    A0 = 0.1d0

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      DO iNodeX = 1, nDOFX

        iNodeX1 = NodeNumberTableX(1,iNodeX)
        iNodeX2 = NodeNumberTableX(2,iNodeX)

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
        X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

        ! --- Top ---
        IF( X2 .GT. 0.0d0 )THEN
          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = rho0 + rho1 * TANH( ( X2 - X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = Vshear      * TANH( ( X2 - X2_Offset ) / a )

          ! --- This is where the typo is. The following expression is
          !     taken from Radice & Rezzolla, 2012, AA, 547, A26, Eq. (48) ---
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = A0 * Vshear * SIN( 2.0d0 * Pi * X1 ) &
                * EXP( -( ( X2 - X2_Offset ) / sigma )**2 )

        ! --- Bottom ---
        ELSE
          uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
            = rho0 - rho1 * TANH( ( X2 + X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V1) &
            = -Vshear     * TANH( ( X2 + X2_Offset ) / a )
          uPF(iNodeX,iX1,iX2,iX3,iPF_V2) &
            = -A0 * Vshear * SIN( 2.0d0 * Pi * X1 ) &
                * EXP( -( ( X2 + X2_Offset ) / sigma )**2 )

         END IF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = 0.0d0
        uAF(iNodeX,iX1,iX2,iX3,iAF_P)  = 1.0d0
        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = uAF(iNodeX,iX1,iX2,iX3,iAF_P) / ( Gamma_IDEAL - One )

      END DO

      CALL ComputeConserved_Euler_Relativistic &
             ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_V1), &
               uPF(:,iX1,iX2,iX3,iPF_V2), uPF(:,iX1,iX2,iX3,iPF_V3), &
               uPF(:,iX1,iX2,iX3,iPF_E ), uPF(:,iX1,iX2,iX3,iPF_Ne), &
               uCF(:,iX1,iX2,iX3,iCF_D ), uCF(:,iX1,iX2,iX3,iCF_S1), &
               uCF(:,iX1,iX2,iX3,iCF_S2), uCF(:,iX1,iX2,iX3,iCF_S3), &
               uCF(:,iX1,iX2,iX3,iCF_E ), uCF(:,iX1,iX2,iX3,iCF_Ne), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_11), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_22), &
               uGF(:,iX1,iX2,iX3,iGF_Gm_dd_33), &
               uAF(:,iX1,iX2,iX3,iAF_P) )

    END DO
    END DO
    END DO


  END SUBROUTINE InitializeFields_KelvinHelmholtzInstability


  SUBROUTINE InitializeFields_StandingAccretionShock &
    ( MassPNS, ShockRadius, AccretionRate, MachNumber, &
      ApplyPerturbation, PerturbationOrder, PerturbationAmplitude, &
      rPerturbationInner, rPerturbationOuter )

    REAL(DP), INTENT(in) :: MassPNS, ShockRadius, AccretionRate, MachNumber
    LOGICAL,  INTENT(in) :: ApplyPerturbation
    INTEGER,  INTENT(in) :: PerturbationOrder
    REAL(DP), INTENT(in) :: PerturbationAmplitude
    REAL(DP), INTENT(in) :: rPerturbationInner
    REAL(DP), INTENT(in) :: rPerturbationOuter

    INTEGER  :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3, iNodeX
    INTEGER  :: iX1_1, iX1_2, iNodeX1_1, iNodeX1_2
    REAL(DP) :: X1_1, X1_2, D_1, D_2, V_1, V_2, P_2
    REAL(DP) :: Alpha, Psi, V0, VSq, W
    REAL(DP) :: X1, X2, dX1, PolytropicConstant, MassConstant
    REAL(DP) :: D(1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: V(1:nNodesX(1),iX_B1(1):iX_E1(1))
    REAL(DP) :: P(1:nNodesX(1),iX_B1(1):iX_E1(1))
    LOGICAL  :: FirstPreShockElement = .FALSE.

    WRITE(*,*)
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Shock radius:   ', ShockRadius / Kilometer, ' km'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'PNS Mass:       ', MassPNS / SolarMass, ' Msun'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Accretion Rate: ', AccretionRate / ( SolarMass / Second ), &
      ' Msun/s'
    WRITE(*,'(6x,A,ES9.2E3)') &
      'Mach number:    ', MachNumber
    WRITE(*,*)
    WRITE(*,'(6x,A,L)') &
      'Apply Perturbation: ', ApplyPerturbation
    WRITE(*,'(6x,A,I1)') &
      'Perturbation order: ', PerturbationOrder
    WRITE(*,'(6x,A,ES9.2E3)') &
      'Perturbation amplitude: ', PerturbationAmplitude
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Inner radius of perturbation: ', rPerturbationInner / Kilometer, ' km'
    WRITE(*,'(6x,A,ES9.2E3,A)') &
      'Outer radius of perturbation: ', rPerturbationOuter / Kilometer, ' km'

    !  --- Locate first element of un-shocked fluid ---

    DO iX1 = iX_B1(1), iX_E1(1)
      DO iNodeX1 = 1, nNodesX(1)

        dX1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 ) - X1
        X1  = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. ShockRadius ) CYCLE

        IF( X1 .GT. ShockRadius .AND. .NOT. FirstPreShockElement )THEN

          iX1_1     = iX1
          iNodeX1_1 = iNodeX1
          X1_1      = X1
          X1_2      = X1 - dX1

          IF( iNodeX1_1 .EQ. 1 )THEN

            iX1_2     = iX1_1 - 1
            iNodeX1_2 = nNodesX(1)

          ELSE

            iX1_2     = iX1_1
            iNodeX1_2 = iNodeX1_1

          END IF

          FirstPreShockElement = .TRUE.

        END IF

      END DO
    END DO

    ! --- Compute fields, pre-shock ---

    DO iX1 = iX_E1(1), iX1_1, -1
      DO iNodeX1 = nNodesX(1), 1, -1

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .LE. ShockRadius ) CYCLE

        Alpha = LapseFunction  ( X1, MassPNS )
        Psi   = ConformalFactor( X1, MassPNS )

        V(iNodeX1,iX1) &
          = -Psi**(-2) * SpeedOfLight * SQRT( One - Alpha**2 )

        D(iNodeX1,iX1) &
          = Psi**(-6) * AccretionRate &
              / ( FourPi * X1**2 * ABS( V(iNodeX1,iX1) ) )

        VSq = Psi**4 * V(iNodeX1,iX1)**2

        P(iNodeX1,iX1) &
          = D(iNodeX1,iX1) * VSq &
              / ( Gamma_IDEAL * MachNumber**2 ) &
              / ( One - ( VSq / SpeedOfLight**2 ) &
              / ( MachNumber**2 * ( Gamma_IDEAL - One ) ) )

      END DO
    END DO

    ! --- Apply jump conditions ---

    D_1 = D(iNodeX1_1,iX1_1)
    V_1 = V(iNodeX1_1,iX1_1)

    CALL ApplyJumpConditions &
           ( iX1_1, iNodeX1_1, X1_1, D_1, V_1, &
             iX1_2, iNodeX1_2, X1_2, &
             D_2, V_2, P_2, MassPNS, PolytropicConstant )

    WRITE(*,*)
    WRITE(*,'(6x,A)') 'Shock location:'
    WRITE(*,'(8x,A)') 'Pre-shock:'
    WRITE(*,'(10x,A,I4.4)')       'iX1     = ', iX1_1
    WRITE(*,'(10x,A,I2.2)')       'iNodeX1 = ', iNodeX1_1
    WRITE(*,'(10x,A,ES13.6E3,A)') 'X1      = ', X1_1 / Kilometer, ' km'
    WRITE(*,'(8x,A)') 'Post-shock:'
    WRITE(*,'(10x,A,I4.4)')       'iX1     = ', iX1_2
    WRITE(*,'(10x,A,I2.2)')       'iNodeX1 = ', iNodeX1_2
    WRITE(*,'(10x,A,ES13.6E3,A)') 'X1      = ', X1_2 / Kilometer, ' km'
    WRITE(*,*)
    WRITE(*,'(6x,A,ES13.6E3)') &
      'Compression Ratio LOG10(D_2/D_1) = ', LOG( D_2 / D_1 ) / LOG( 1.0d1 )
    WRITE(*,*)

    ! --- Compute fields, post-shock ---

    Alpha = LapseFunction  ( X1_1, MassPNS )
    Psi   = ConformalFactor( X1_1, MassPNS )
    W     = LorentzFactor( Psi, V_1 )

    MassConstant = Psi**6 * Alpha * X1_1**2 * D_1 * W * V_1

    V0 = V_2

    DO iX1 = iX1_2, iX_B1(1), -1
      DO iNodeX1 = nNodesX(1), 1, -1

        X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )

        IF( X1 .GT. ShockRadius ) CYCLE

        Alpha = LapseFunction  ( X1, MassPNS )
        Psi   = ConformalFactor( X1, MassPNS )

        CALL NewtonRaphson_PostShockVelocity &
               ( Alpha, Psi, MassConstant, PolytropicConstant, &
                 MassPNS, AccretionRate, X1, V0  )

        V(iNodeX1,iX1) = V0

        W = LorentzFactor( Psi, V0 )

        D(iNodeX1,iX1) &
          = MassConstant / ( Psi**6 * Alpha * X1**2  * W * V0 )

        P(iNodeX1,iX1) &
          = PolytropicConstant * D(iNodeX1,iX1)**( Gamma_IDEAL )

      END DO
    END DO

    ! --- Map to 3D domain ---

    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)

        iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )

        IF( ApplyPerturbation )THEN

          X1 = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
          X2 = NodeCoordinate( MeshX(2), iX2, iNodeX2 )

          IF( X1 .GE. rPerturbationInner .AND. X1 .LE. rPerturbationOuter )THEN

            IF( PerturbationOrder .EQ. 0 ) &
              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D(iNodeX1,iX1) * ( One + PerturbationAmplitude )

            IF( PerturbationOrder .EQ. 1 ) &
              uPF(iNodeX,iX1,iX2,iX3,iPF_D) &
                = D(iNodeX1,iX1) * ( One + PerturbationAmplitude * COS( X2 ) )

          ELSE

            uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D(iNodeX1,iX1)

          END IF

        ELSE

          uPF(iNodeX,iX1,iX2,iX3,iPF_D) = D(iNodeX1,iX1)

        END IF

        uPF(iNodeX,iX1,iX2,iX3,iPF_V1) = V(iNodeX1,iX1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_V2) = Zero
        uPF(iNodeX,iX1,iX2,iX3,iPF_V3) = Zero
        uAF(iNodeX,iX1,iX2,iX3,iAF_P ) = P(iNodeX1,iX1)
        uPF(iNodeX,iX1,iX2,iX3,iPF_E) &
          = uAF(iNodeX,iX1,iX2,iX3,iAF_P ) / ( Gamma_IDEAL - One )

        CALL ComputeConserved_Euler_Relativistic &
               ( uPF(iNodeX,iX1,iX2,iX3,iPF_D ),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V1),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V2),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_V3),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_E ),       &
                 uPF(iNodeX,iX1,iX2,iX3,iPF_Ne),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_D ),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S1),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S2),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_S3),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_E ),       &
                 uCF(iNodeX,iX1,iX2,iX3,iCF_Ne),       &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                 uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33), &
                 uAF(iNodeX,iX1,iX2,iX3,iAF_P ) )

      END DO
      END DO
      END DO
    END DO
    END DO
    END DO

  END SUBROUTINE InitializeFields_StandingAccretionShock


  ! --- Auxiliary utilities for standing accretion shock problem ---

  SUBROUTINE ApplyJumpConditions &
    ( iX1_1, iNodeX1_1, X1_1, D_1, V_1, &
      iX1_2, iNodeX1_2, X1_2, &
      D_2, V_2, P_2, MassPNS, PolytropicConstant )

    INTEGER,  INTENT(in)  :: iX1_1, iNodeX1_1, iX1_2, iNodeX1_2
    REAL(DP), INTENT(in)  :: X1_1, X1_2, D_1, V_1, MassPNS
    REAL(DP), INTENT(out) :: D_2, V_2, P_2, PolytropicConstant

    REAL(DP) :: Alpha, Psi
    REAL(DP) :: C1, C2, C3, a0, a1, a2, a3, a4, X1
    REAL(DP) :: W

    REAL(DP), PARAMETER :: ShockTolerance = 0.1_DP
    LOGICAL             :: FoundShockVelocity = .FALSE.

    ! --- Constants from three jump conditions ---

    Alpha = LapseFunction  ( X1_1, MassPNS )
    Psi   = ConformalFactor( X1_1, MassPNS )

    C1 = D_1 * V_1 / Alpha

    C2 = D_1 * SpeedOfLight**2 / Alpha**2 * ( V_1 / SpeedOfLight )**2

    C3 = D_1 * SpeedOfLight**2 / Alpha**2 * V_1

    ! --- Five constants for post-shock fluid-velocity ---

    a4 = Psi**8 &
          * One / ( Gamma_IDEAL - One )**2 * C3**2 / SpeedOfLight**6
    a3 = -Two * Psi**8 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One )**2 * C2 * C3 / SpeedOfLight**4
    a2 = Psi**4 &
          / SpeedOfLight**2 * ( Psi**4 * Gamma_IDEAL**2 &
          / ( Gamma_IDEAL - One )**2 * C2**2 + Two * One &
          / ( Gamma_IDEAL - One ) &
          * C3**2 / SpeedOfLight**2 + C1**2 * SpeedOfLight**2 )
    a1 = -Two * Psi**4 &
          * Gamma_IDEAL / ( Gamma_IDEAL - One ) * C2 * C3 / SpeedOfLight**2
    a0 = One / SpeedOfLight**2 * ( C3**2 - C1**2 * SpeedOfLight**4 )

    ! --- Newton-Raphson method for post-shock fluid-velocity ---

    V_2 = Two * V_1

    ! --- Ensure that shocked velocity is obtained ---

    FoundShockVelocity = .FALSE.
    DO WHILE( .NOT. FoundShockVelocity )

      V_2 = Half * V_2
      CALL NewtonRaphson_JumpConditions( a0, a1, a2, a3, a4, V_2 )

      IF( ABS( V_2 - V_1 ) / ABS( V_1 ) .GT. ShockTolerance ) &
        FoundShockVelocity = .TRUE.

    END DO

    ! --- Post-shock density, velocity, pressure, and polytropic constant ---

    Psi = ConformalFactor( X1_2, MassPNS )
    W   = LorentzFactor( Psi, V_2 )

    D_2 = ABS( C1 ) * SQRT( One / V_2**2 - Psi**4 / SpeedOfLight**2 )

    P_2 = ( Gamma_IDEAL - One ) / Gamma_IDEAL &
            * ( C3 - D_2 * SpeedOfLight**2 * W**2 * V_2 ) / ( W**2 * V_2 )

    PolytropicConstant = P_2 / D_2**( Gamma_IDEAL )

  END SUBROUTINE ApplyJumpConditions


  SUBROUTINE NewtonRaphson_JumpConditions( a0, a1, a2, a3, a4, V )

    REAL(DP), INTENT(in)    :: a0, a1, a2, a3, a4
    REAL(DP), INTENT(inout) :: V

    REAL(DP) :: f, df, dV
    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION

    INTEGER,  PARAMETER :: MAX_ITER = 10
    REAL(DP), PARAMETER :: TOLERANCE = 1.0d-15

    CONVERGED = .FALSE.
    ITERATION = 0
    DO WHILE( .NOT. CONVERGED .AND. ITERATION .LT. MAX_ITER )

      ITERATION = ITERATION + 1

      f  = a4 * V**4 + a3 * V**3 + a2 * V**2 + a1 * V + a0
      df = Four * a4 * V**3 + Three * a3 * V**2 + Two * a2 * V + a1

      dV = -f / df
      V = V + dV

      IF( ABS( dV / V ) .LT. TOLERANCE )THEN
        CONVERGED = .TRUE.
      END IF

    END DO

  END SUBROUTINE NewtonRaphson_JumpConditions


  SUBROUTINE NewtonRaphson_PostShockVelocity &
    ( Alpha, Psi, MassConstant, PolytropicConstant, &
      MassPNS, AccretionRate, X1, V )

    REAL(DP), INTENT(in)    :: Alpha, Psi, MassConstant, &
                               PolytropicConstant, MassPNS, AccretionRate, X1
    REAL(DP), INTENT(inout) :: V

    REAL(DP) :: f, df, dV, W
    INTEGER  :: ITERATION
    LOGICAL  :: CONVERGED

    INTEGER,  PARAMETER :: MAX_ITER = 20
    REAL(DP), PARAMETER :: TOLERANCE = 1.0d-15

    CONVERGED = .FALSE.
    ITERATION = 0
    DO WHILE( .NOT. CONVERGED .AND. ITERATION .LT. MAX_ITER )

      ITERATION = ITERATION + 1

      W = LorentzFactor( Psi, V )

      f  = Gamma_IDEAL / ( Gamma_IDEAL - One ) &
             * PolytropicConstant / SpeedOfLight**2 * ( MassConstant &
             / ( Psi**6 * Alpha * X1**2 * W * V ) )**( Gamma_IDEAL - One ) &
             - One / ( Alpha * W ) + One

      df = -Gamma_IDEAL * PolytropicConstant / SpeedOfLight**2 &
             * ( MassConstant &
                 / ( Psi**6 * Alpha * X1**2 * W * V ) )**( Gamma_IDEAL - One ) &
                 * ( Psi**4 * V / SpeedOfLight**2 * W**2 + One / V ) &
                 + W / Alpha * Psi**4 * V / SpeedOfLight**2

      dV = -f / df
      V = V + dV

      IF( ABS( dV / V ) .LT. TOLERANCE ) &
        CONVERGED = .TRUE.

    END DO

  END SUBROUTINE NewtonRaphson_PostShockVelocity


  REAL(DP) FUNCTION LapseFunction( R, M )

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    LapseFunction = ABS( ( MAX( ABS( R ), SqrtTiny ) - Half * M ) &
                       / ( MAX( ABS( R ), SqrtTiny ) + Half * M ) )

    RETURN
  END FUNCTION LapseFunction


  REAL(DP) FUNCTION ConformalFactor( R, M )

    REAL(DP), INTENT(in) :: R, M

    ! --- Schwarzschild Metric in Isotropic Coordinates ---

    ConformalFactor = One + Half * M / MAX( ABS( R ), SqrtTiny )

    RETURN
  END FUNCTION ConformalFactor


  REAL(DP) FUNCTION LorentzFactor( Psi, V )

    REAL(DP), INTENT(in) :: Psi, V

    LorentzFactor = One / SQRT( One - Psi**4 * ( V / SpeedOfLight )**2 )

    RETURN
  END FUNCTION LorentzFactor


END MODULE InitializationModule_Relativistic
