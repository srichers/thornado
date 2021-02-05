MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    One, &
    Two
  USE ProgramHeaderModule, ONLY: &
    nDOFX, &
    nDimsX, &
    nNodes, &
    nNodesX, &
    bcX
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply, &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    WeightsX_q, &
    NodeNumberTableX
  USE UtilitiesModule, ONLY: &
    MinModB, &
    NodeNumberX
  USE PolynomialBasisModule_Legendre, ONLY: &
    P_X1, &
    P_X2, &
    P_X3, &
    IndPX_Q, &
    MassPX
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Fluid, &
    MapModalToNodal_Fluid, &
    Kij_X, &
    Pij_X
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_Alpha, &
    iGF_Beta_1, &
    iGF_Beta_2, &
    iGF_Beta_3, &
    iGF_SqrtGm
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E, &
    iCF_Ne, &
    iDF_TCI
  USE Euler_BoundaryConditionsModule, ONLY: &
    ApplyInnerBC_Euler, &
    ApplyOuterBC_Euler, &
    iApplyBC_Euler_Both, &
    ApplyBoundaryConditions_Euler
  USE Euler_CharacteristicDecompositionModule_Relativistic_IDEAL, ONLY: &
    ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL
  USE Euler_DiscontinuityDetectionModule, ONLY: &
    InitializeTroubledCellIndicator_Euler, &
    FinalizeTroubledCellIndicator_Euler, &
    UseTroubledCellIndicator, &
    LimiterThreshold, &
    DetectTroubledCells_Euler
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_SlopeLimiter, &
    Timer_Euler_CharacteristicDecomposition

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: FinalizeSlopeLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: ApplySlopeLimiter_Euler_Relativistic_IDEAL

  LOGICAL      :: UseSlopeLimiter
  LOGICAL      :: UseCharacteristicLimiting
  LOGICAL      :: UseConservativeCorrection
  CHARACTER(4) :: SlopeLimiterMethod

  ! --- WENO Limiter ---
  REAL(DP), ALLOCATABLE :: OrthonormalBasis(:,:,:)
  REAL(DP), ALLOCATABLE :: VandermondeMatrix(:,:)

  ! --- TVD Limiter ---
  REAL(DP) :: BetaTVD, BetaTVB
  REAL(DP) :: SlopeTolerance
  REAL(DP) :: LimiterThresholdParameter
  REAL(DP) :: I_6x6(1:6,1:6)

  REAL(DP), ALLOCATABLE :: LegendreX(:,:)

#if defined(THORNADO_OMP_OL)
  !$OMP DECLARE TARGET &
  !$OMP   ( UseSlopeLimiter, UseCharacteristicLimiting, &
  !$OMP     UseConservativeCorrection, &
  !$OMP     BetaTVD, BetaTVB, SlopeTolerance, I_6x6, LegendreX )
#elif defined(THORNADO_OACC)
  !$ACC DECLARE CREATE &
  !$ACC   ( UseSlopeLimiter, UseCharacteristicLimiting, &
  !$ACC     UseConservativeCorrection, &
  !$ACC     BetaTVD, BetaTVB, SlopeTolerance, I_6x6, LegendreX )
#endif


CONTAINS


  SUBROUTINE InitializeSlopeLimiter_Euler_Relativistic_IDEAL &
    ( UseSlopeLimiter_Option,           &
      SlopeLimiterMethod_Option,        &
      BetaTVD_Option,                   &
      BetaTVB_Option,                   &
      SlopeTolerance_Option,            &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option,  &
      LimiterThresholdParameter_Option, &
      UseConservativeCorrection_Option, &
      Verbose_Option )

    REAL(DP), INTENT(in),     OPTIONAL :: &
      BetaTVD_Option, BetaTVB_Option, &
      SlopeTolerance_Option
    LOGICAL,  INTENT(in),     OPTIONAL :: &
      UseSlopeLimiter_Option, &
      UseCharacteristicLimiting_Option, &
      UseTroubledCellIndicator_Option,  &
      UseConservativeCorrection_Option, &
      Verbose_Option
    REAL(DP), INTENT(in),     OPTIONAL :: &
      LimiterThresholdParameter_Option
    CHARACTER(*), INTENT(in), OPTIONAL :: &
      SlopeLimiterMethod_Option

    INTEGER :: i, iPol, iNX, iNX1, iNX2, iNX3
    LOGICAL :: Verbose

    UseSlopeLimiter = .TRUE.
    IF( PRESENT( UseSlopeLimiter_Option ) ) &
      UseSlopeLimiter = UseSlopeLimiter_Option

    SlopeLimiterMethod = 'TVD'
    IF( PRESENT( SlopeLimiterMethod_Option ) ) &
      SlopeLimiterMethod = SlopeLimiterMethod_Option

    BetaTVD = One
    IF( PRESENT( BetaTVD_Option ) ) &
      BetaTVD = BetaTVD_Option

    BetaTVB = Zero
    IF( PRESENT( BetaTVB_Option ) ) &
      BetaTVB = BetaTVB_Option

    SlopeTolerance = 1.0e-3_DP
    IF( PRESENT( SlopeTolerance_Option ) ) &
      SlopeTolerance = SlopeTolerance_Option

    UseCharacteristicLimiting = .FALSE.
    IF( PRESENT( UseCharacteristicLimiting_Option ) ) &
      UseCharacteristicLimiting = UseCharacteristicLimiting_Option

    UseTroubledCellIndicator = .TRUE.
    IF( PRESENT( UseTroubledCellIndicator_Option ) ) &
      UseTroubledCellIndicator = UseTroubledCellIndicator_Option

    LimiterThresholdParameter = 0.03_DP
    IF( PRESENT( LimiterThresholdParameter_Option ) ) &
      LimiterThresholdParameter = LimiterThresholdParameter_Option
    LimiterThreshold = LimiterThresholdParameter * Two**( nNodes - 2 )

    UseConservativeCorrection = .TRUE.
    IF( PRESENT( UseConservativeCorrection_Option ) ) &
      UseConservativeCorrection = UseConservativeCorrection_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Slope Limiter (Euler, Relativistic, IDEAL)'
      WRITE(*,'(A)') &
        '    ------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseSlopeLimiter: ' , &
        UseSlopeLimiter
      WRITE(*,*)
      WRITE(*,'(A4,A27,A)')         '', 'SlopeLimiterMethod: ', &
        TRIM( SlopeLimiterMethod )
      WRITE(*,*)

      IF( TRIM( SlopeLimiterMethod ) .EQ. 'TVD' )THEN

        WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'BetaTVD: ' , &
          BetaTVD
        WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'BetaTVB: ' , &
          BetaTVB
        WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'SlopeTolerance: ' , &
          SlopeTolerance
        WRITE(*,*)

      END IF

      WRITE(*,'(A4,A27,L1)'       ) '', 'UseCharacteristicLimiting: ' , &
        UseCharacteristicLimiting
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseTroubledCellIndicator: ' , &
        UseTroubledCellIndicator
      WRITE(*,'(A4,A27,ES10.3E3)' ) '', 'LimiterThreshold: ' , &
        LimiterThreshold
      WRITE(*,*)
      WRITE(*,'(A4,A27,L1)'       ) '', 'UseConservativeCorrection: ' , &
        UseConservativeCorrection
    END IF

    I_6x6 = Zero
    DO i = 1, 6
      I_6x6(i,i) = One
    END DO

    ALLOCATE( LegendreX(1:nDOFX,1:nDOFX) )

    DO iPol = 1, nDOFX ! Only need for iPol = 2,3,4 (FIXME)

      DO iNX3 = 1, nNodesX(3)
      DO iNX2 = 1, nNodesX(2)
      DO iNX1 = 1, nNodesX(1)

        iNX = NodeNumberX( iNX1, iNX2, iNX3 )

        LegendreX(iNX,iPol) &
          = P_X1  (IndPX_Q(1,iPol)) % P( MeshX(1) % Nodes(iNX1) ) &
            * P_X2(IndPX_Q(2,iPol)) % P( MeshX(2) % Nodes(iNX2) ) &
            * P_X3(IndPX_Q(3,iPol)) % P( MeshX(3) % Nodes(iNX3) )

      END DO
      END DO
      END DO

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET UPDATE TO &
    !$OMP   ( UseSlopeLimiter, UseCharacteristicLimiting, &
    !$OMP     UseConservativeCorrection, &
    !$OMP     BetaTVD, BetaTVB, SlopeTolerance, I_6x6, LegendreX )
#elif defined(THORNADO_OACC)
    !$ACC UPDATE DEVICE &
    !$ACC   ( UseSlopeLimiter, UseCharacteristicLimiting, &
    !$ACC     UseConservativeCorrection, &
    !$ACC     BetaTVD, BetaTVB, SlopeTolerance, I_6x6, LegendreX )
#endif

    CALL InitializeTroubledCellIndicator_Euler

    IF( TRIM( SlopeLimiterMethod ).EQ. 'WENO' )THEN

      CALL InitializeSlopeLimiter_Euler_WENO

    END IF

  END SUBROUTINE InitializeSlopeLimiter_Euler_Relativistic_IDEAL


!!$  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL &
!!$    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )
!!$
!!$    INTEGER,  INTENT(in)           :: &
!!$      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
!!$    REAL(DP), INTENT(in)           :: &
!!$      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
!!$    REAL(DP), INTENT(inout)        :: &
!!$      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
!!$    REAL(DP), INTENT(inout)        :: &
!!$      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
!!$    LOGICAL,  INTENT(in), OPTIONAL :: &
!!$      SuppressBC_Option
!!$    INTEGER,  INTENT(in), OPTIONAL :: &
!!$      iApplyBC_Option(3)
!!$
!!$    SELECT CASE( TRIM( SlopeLimiterMethod ) )
!!$
!!$      CASE( 'TVD' )
!!$
!!$        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
!!$               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
!!$                 SuppressBC_Option, iApplyBC_Option )
!!$
!!$      CASE( 'WENO' )
!!$
!!$        IF( UseCharacteristicLimiting )THEN
!!$
!!$          CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO_Characteristic &
!!$                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
!!$                   SuppressBC_Option, iApplyBC_Option )
!!$
!!$        ELSE
!!$
!!$          CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO_ComponentWise &
!!$                 ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
!!$                   SuppressBC_Option, iApplyBC_Option )
!!$
!!$        END IF
!!$
!!$      CASE DEFAULT
!!$
!!$        CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD &
!!$               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
!!$                 SuppressBC_Option, iApplyBC_Option )
!!$
!!$    END SELECT
!!$
!!$  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    IF( .NOT. UseCharacteristicLimiting )THEN

      CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD_Componentwise &
             ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
               SuppressBC_Option, iApplyBC_Option )

    ELSE

      CALL ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD_Characteristic &
               ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, &
                 SuppressBC_Option, iApplyBC_Option )
    END IF

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD_Componentwise &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    LOGICAL  :: LimitedCell(nCF,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    LOGICAL  :: SuppressBC
    LOGICAL  :: ExcludeInnerGhostCell(3), ExcludeOuterGhostCell(3)
    INTEGER  :: iNX, iX1, iX2, iX3, iCF
    INTEGER  :: iApplyBC(3)
    REAL(DP) :: SlopeDifference(nCF)
    REAL(DP) :: a(nDimsX), b(nDimsX), c(nDimsX)
    REAL(DP) :: dU(nCF,nDimsX)

    INTEGER  :: nX_B, nX_E, nCF_B, nCF_E

    REAL(DP) :: SqrtGm(1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                     iX_B0(2):iX_E0(2), &
                                     iX_B0(3):iX_E0(3))
    REAL(DP) :: Vol   (              iX_B0(1):iX_E0(1), &
                                     iX_B0(2):iX_E0(2), &
                                     iX_B0(3):iX_E0(3))
    REAL(DP) :: U_X   (1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                     iX_B0(2):iX_E0(2), &
                                     iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K   (1:nCF        ,iX_B0(1):iX_E0(1), &
                                     iX_B0(2):iX_E0(2), &
                                     iX_B0(3):iX_E0(3))

    REAL(DP) :: U_N   (1:nDOFX,1:nCF,iX_B1(1):iX_E1(1), &
                                     iX_B1(2):iX_E1(2), &
                                     iX_B1(3):iX_E1(3))
    REAL(DP) :: U_M   (1:nDOFX,1:nCF,iX_B1(1):iX_E1(1), &
                                     iX_B1(2):iX_E1(2), &
                                     iX_B1(3):iX_E1(3))

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    nX_B  = PRODUCT( iX_E0 - iX_B0 + 1 )
    nX_E  = PRODUCT( iX_E1 - iX_B1 + 1 )
    nCF_B = nCF * nX_B
    nCF_E = nCF * nX_E

    iApplyBC = iApplyBC_Euler_Both
    IF( PRESENT( iApplyBC_Option ) ) &
       iApplyBC = iApplyBC_Option

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    ExcludeInnerGhostCell = .FALSE.
    ExcludeOuterGhostCell = .FALSE.
    IF( ApplyInnerBC_Euler( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeInnerGhostCell(1) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeOuterGhostCell(1) = .TRUE.
    IF( ApplyInnerBC_Euler( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeInnerGhostCell(2) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeOuterGhostCell(2) = .TRUE.
    IF( ApplyInnerBC_Euler( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeInnerGhostCell(3) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeOuterGhostCell(3) = .TRUE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$OMP             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$OMP MAP( alloc: SqrtGm, Vol, U_X, U_K, U_N, U_M, LimitedCell )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$ACC             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$ACC CREATE(     SqrtGm, Vol, U_X, U_K, U_N, U_M, LimitedCell )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_N, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_N(iNX,iCF,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, SqrtGm, U_X, G, U_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGm(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

      U_X(iNX,iCF_D ,iX1,iX2,iX3) &
        = U_N(iNX,iCF_D ,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_S1,iX1,iX2,iX3) &
        = U_N(iNX,iCF_S1,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_S2,iX1,iX2,iX3) &
        = U_N(iNX,iCF_S2,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_S3,iX1,iX2,iX3) &
        = U_N(iNX,iCF_S3,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_E ,iX1,iX2,iX3) &
        = U_N(iNX,iCF_E ,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)
      U_X(iNX,iCF_Ne,iX1,iX2,iX3) &
        = U_N(iNX,iCF_Ne,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    ! --- Compute volumes of compute cells ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nX_B , One, SqrtGm, nDOFX, &
             WeightsX_q, 1, Zero, Vol, 1 )

    ! --- Compute cell integrals ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_B, One, U_X   , nDOFX, &
             WeightsX_q, 1, Zero, U_K, 1 )

    ! --- Form cell averages ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K, Vol )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      U_K(iCF,iX1,iX2,iX3) = U_K(iCF,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO

    ! --- Map Nodal to Modal ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_E, nDOFX, One, Kij_X, nDOFX, &
             U_N, nDOFX, Zero, U_M, nDOFX )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_M, MassPX )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_M(iNX,iCF,iX1,iX2,iX3) = U_M(iNX,iCF,iX1,iX2,iX3) * MassPX(iNX)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( dU, SlopeDifference, a, b, c )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, D, U_M, dX1, dX2, dX3, LimitedCell, &
    !$ACC          ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$ACC PRIVATE( dU, SlopeDifference, a, b, c )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4) &
    !$OMP PRIVATE( dU, SlopeDifference )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      a(1) = U_M(2,iCF,iX1,iX2,iX3)

      IF     ( iX1 .EQ. iX_B0(1) .AND. ExcludeInnerGhostCell(1) )THEN

        c(1) = BetaTVD * ( U_M(1,iCF,iX1+1,iX2,iX3) &
                         - U_M(1,iCF,iX1  ,iX2,iX3) )

        b(1) = c(1)

      ELSE IF( iX1 .EQ. iX_E0(1) .AND. ExcludeOuterGhostCell(1) )THEN

        b(1) = BetaTVD * ( U_M(1,iCF,iX1  ,iX2,iX3) &
                         - U_M(1,iCF,iX1-1,iX2,iX3) )
        c(1) = b(1)


      ELSE

        b(1) = BetaTVD * ( U_M(1,iCF,iX1  ,iX2,iX3) &
                         - U_M(1,iCF,iX1-1,iX2,iX3) )
        c(1) = BetaTVD * ( U_M(1,iCF,iX1+1,iX2,iX3) &
                         - U_M(1,iCF,iX1  ,iX2,iX3) )

      END IF

      dU(iCF,1) = MinModB( a(1), b(1), c(1), dX1(iX1), BetaTVB )

      IF( nDimsX .GT. 1 )THEN

        a(2) = U_M(3,iCF,iX1,iX2,iX3)

        IF     ( iX2 .EQ. iX_B0(2) .AND. ExcludeInnerGhostCell(2) )THEN

          c(2) = BetaTVD * ( U_M(1,iCF,iX1,iX2+1,iX3) &
                           - U_M(1,iCF,iX1,iX2  ,iX3) )

          b(2) = c(2)

        ELSE IF( iX2 .EQ. iX_E0(2) .AND. ExcludeOuterGhostCell(2) )THEN

          b(2) = BetaTVD * ( U_M(1,iCF,iX1,iX2  ,iX3) &
                           - U_M(1,iCF,iX1,iX2-1,iX3) )
          c(2) = b(2)


        ELSE

          b(2) = BetaTVD * ( U_M(1,iCF,iX1,iX2  ,iX3) &
                           - U_M(1,iCF,iX1,iX2-1,iX3) )
          c(2) = BetaTVD * ( U_M(1,iCF,iX1,iX2+1,iX3) &
                           - U_M(1,iCF,iX1,iX2  ,iX3) )

        END IF

        dU(iCF,2) = MinModB( a(2), b(2), c(2), dX2(iX2), BetaTVB )

      END IF

      IF( nDimsX .GT. 2 )THEN

        a(3) = U_M(4,iCF,iX1,iX2,iX3)

        IF     ( iX3 .EQ. iX_B0(3) .AND. ExcludeInnerGhostCell(3) )THEN

          c(3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3+1) &
                           - U_M(1,iCF,iX1,iX2,iX3  ) )

          b(3) = c(3)

        ELSE IF( iX3 .EQ. iX_E0(3) .AND. ExcludeOuterGhostCell(3) )THEN

          b(3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3  ) &
                        - U_M(1,iCF,iX1,iX2,iX3-1) )
          c(3) = b(3)


        ELSE

          b(3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3  ) &
                           - U_M(1,iCF,iX1,iX2,iX3-1) )
          c(3) = BetaTVD * ( U_M(1,iCF,iX1,iX2,iX3+1) &
                           - U_M(1,iCF,iX1,iX2,iX3  ) )

        END IF

        dU(iCF,3) = MinModB( a(3), b(3), c(3), dX3(iX3), BetaTVB )

      END IF

      ! --- Compare Limited Slopes to Original Slopes ---

      SlopeDifference(iCF) = ABS( U_M(2,iCF,iX1,iX2,iX3) - dU(iCF,1) )

      IF( nDimsX .GT. 1 ) &
        SlopeDifference(iCF) &
          = MAX( SlopeDifference(iCF), &
                 ABS( U_M(3,iCF,iX1,iX2,iX3) - dU(iCF,2) ) )

      IF( nDimsX .GT. 2 ) &
        SlopeDifference(iCF) &
          = MAX( SlopeDifference(iCF), &
                 ABS( U_M(4,iCF,iX1,iX2,iX3) - dU(iCF,3) ) )

      ! --- Replace Slopes and Discard High-Order Components ---
      ! --- if Limited Slopes Deviate too Much from Original ---

      LimitedCell(iCF,iX1,iX2,iX3) = .FALSE.

      IF( SlopeDifference(iCF) &
            .GT. SlopeTolerance * ABS( U_M(1,iCF,iX1,iX2,iX3) ) )THEN

        DO iNX = 2, nDOFX

          U_M(iNX,iCF,iX1,iX2,iX3) = Zero

        END DO

        U_M(2,iCF,iX1,iX2,iX3) = dU(iCF,1)

        IF( nDimsX .GT. 1 ) U_M(3,iCF,iX1,iX2,iX3) = dU(iCF,2)

        IF( nDimsX .GT. 2 ) U_M(4,iCF,iX1,iX2,iX3) = dU(iCF,3)

        LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

      END IF

    END DO
    END DO
    END DO
    END DO

    ! --- Apply Conservative Correction ---

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    ! --- Map Modal to Nodal ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_E, nDOFX, One, Pij_X, nDOFX, &
             U_M, nDOFX, Zero, U_N, nDOFX )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U, U_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U(iNX,iX1,iX2,iX3,iCF) = U_N(iNX,iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G, D, dX1, dX2, dX3, &
    !$OMP               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$OMP               SqrtGm, Vol, U_X, U_K, U_N, U_M, LimitedCell )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G, D, dX1, dX2, dX3, &
    !$ACC               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$ACC               SqrtGm, Vol, U_X, U_K, U_N, U_M, LimitedCell )
#endif

    END ASSOCIATE

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD_Componentwise


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD_Characteristic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    LOGICAL  :: LimitedCell(nCF,iX_B0(1):iX_E0(1), &
                                iX_B0(2):iX_E0(2), &
                                iX_B0(3):iX_E0(3))
    LOGICAL  :: SuppressBC
    LOGICAL  :: ExcludeInnerGhostCell(3), ExcludeOuterGhostCell(3)
    INTEGER  :: iNX, iX1, iX2, iX3, iCF, jCF, iGF
    INTEGER  :: iApplyBC(3)
    REAL(DP) :: SlopeDifference(nCF)
    REAL(DP) :: a(nCF,nDimsX), b(nCF,nDimsX), c(nCF,nDimsX)
    REAL(DP) :: dU(nCF,nDimsX), dU_C(nCF)

    INTEGER  :: nX_B, nX_E, nCF_B, nCF_E, nGF_B

    REAL(DP) :: SqrtGm (1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: Vol    (              iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: U_X    (1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K    (1:nCF        ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: G_X    (1:nDOFX,1:8  ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: G_K    (1:8          ,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: U_N    (1:nDOFX,1:nCF,iX_B1(1):iX_E1(1), &
                                      iX_B1(2):iX_E1(2), &
                                      iX_B1(3):iX_E1(3))
    REAL(DP) :: U_M    (1:nDOFX,1:nCF,iX_B1(1):iX_E1(1), &
                                      iX_B1(2):iX_E1(2), &
                                      iX_B1(3):iX_E1(3))

    REAL(DP) :: R_X1   (1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X1(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X2   (1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X2(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: R_X3   (1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))
    REAL(DP) :: invR_X3(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                      iX_B0(2):iX_E0(2), &
                                      iX_B0(3):iX_E0(3))

    REAL(DP) :: R(1:nCF,1:nCF), invR(1:nCF,1:nCF)
    REAL(DP) :: GK(8), UK(nCF)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    ASSOCIATE( dX1 => MeshX(1) % Width, &
               dX2 => MeshX(2) % Width, &
               dX3 => MeshX(3) % Width )

    nX_B  = PRODUCT( iX_E0 - iX_B0 + 1 )
    nX_E  = PRODUCT( iX_E1 - iX_B1 + 1 )
    nCF_B = nCF * nX_B
    nCF_E = nCF * nX_E
    nGF_B = 8   * nX_B

    iApplyBC = iApplyBC_Euler_Both
    IF( PRESENT( iApplyBC_Option ) ) &
       iApplyBC = iApplyBC_Option

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    ExcludeInnerGhostCell = .FALSE.
    ExcludeOuterGhostCell = .FALSE.
    IF( ApplyInnerBC_Euler( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeInnerGhostCell(1) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(1) ) .AND. bcX(1) .NE. 1 ) &
      ExcludeOuterGhostCell(1) = .TRUE.
    IF( ApplyInnerBC_Euler( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeInnerGhostCell(2) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(2) ) .AND. bcX(2) .NE. 1 ) &
      ExcludeOuterGhostCell(2) = .TRUE.
    IF( ApplyInnerBC_Euler( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeInnerGhostCell(3) = .TRUE.
    IF( ApplyOuterBC_Euler( iApplyBC(3) ) .AND. bcX(3) .NE. 1 ) &
      ExcludeOuterGhostCell(3) = .TRUE.

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$OMP             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$OMP MAP( alloc: SqrtGm, Vol, U_X, U_K, U_N, U_M, G_X, G_K, LimitedCell, &
    !$OMP             R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, dX1, dX2, dX3, &
    !$ACC             ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$ACC CREATE(     SqrtGm, Vol, U_X, U_K, U_N, U_M, G_X, G_K, LimitedCell, &
    !$ACC             R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#endif

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_N, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_N(iNX,iCF,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, SqrtGm, G, G_X )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGm(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

      G_X(iNX,1,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_11)
      G_X(iNX,2,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_22)
      G_X(iNX,3,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Gm_dd_33)
      G_X(iNX,4,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)
      G_X(iNX,5,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Alpha)
      G_X(iNX,6,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Beta_1)
      G_X(iNX,7,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Beta_2)
      G_X(iNX,8,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_Beta_3)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B0, iX_E0, U_X, U_N, SqrtGm )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_X(iNX,iCF,iX1,iX2,iX3) &
        = U_N(iNX,iCF,iX1,iX2,iX3) * SqrtGm(iNX,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

    ! --- Compute volumes of compute cells ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nX_B , One, SqrtGm, nDOFX, &
             WeightsX_q, 1, Zero, Vol, 1 )

    ! --- Compute cell integrals ---

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nCF_B, One, U_X   , nDOFX, &
             WeightsX_q, 1, Zero, U_K, 1 )

    CALL MatrixVectorMultiply &
           ( 'T', nDOFX, nGF_B, One, G_X   , nDOFX, &
             WeightsX_q, 1, Zero, G_K, 1 )

    ! --- Form cell averages ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, U_K, Vol, G_K )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      U_K(iCF_D ,iX1,iX2,iX3) = U_K(iCF_D ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      U_K(iCF_S1,iX1,iX2,iX3) = U_K(iCF_S1,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      U_K(iCF_S2,iX1,iX2,iX3) = U_K(iCF_S2,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      U_K(iCF_S3,iX1,iX2,iX3) = U_K(iCF_S3,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      U_K(iCF_E ,iX1,iX2,iX3) = U_K(iCF_E ,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      U_K(iCF_Ne,iX1,iX2,iX3) = U_K(iCF_Ne,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

      G_K(1,iX1,iX2,iX3) = G_K(1,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      G_K(2,iX1,iX2,iX3) = G_K(2,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      G_K(3,iX1,iX2,iX3) = G_K(3,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      G_K(4,iX1,iX2,iX3) = G_K(4,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      G_K(5,iX1,iX2,iX3) = G_K(5,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      G_K(6,iX1,iX2,iX3) = G_K(6,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      G_K(7,iX1,iX2,iX3) = G_K(7,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)
      G_K(8,iX1,iX2,iX3) = G_K(8,iX1,iX2,iX3) / Vol(iX1,iX2,iX3)

    END DO
    END DO
    END DO

    ! --- Map Nodal to Modal ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_E, nDOFX, One, Kij_X, nDOFX, &
             U_N, nDOFX, Zero, U_M, nDOFX )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U_M, MassPX )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_M(iNX,iCF,iX1,iX2,iX3) = U_M(iNX,iCF,iX1,iX2,iX3) * MassPX(iNX)

    END DO
    END DO
    END DO
    END DO
    END DO



! Isolate all Matrix-Vector multiplies into one block

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, D, G_K, U_M, &
    !$ACC          R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 ) &
    !$ACC PRIVATE( R, invR, UK, GK )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      DO iGF = 1, 8

        GK(iGF) = G_K(iGF,iX1,iX2,iX3)

      END DO

      DO iCF = 1, nCF

        UK(iCF) = U_M(1,iCF,iX1,iX2,iX3)

      END DO

      CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
             ( 1, GK, UK, R, invR )

      DO iCF = 1, nCF
      DO jCF = 1, nCF

        R_X1   (iCF,jCF,iX1,iX2,iX3) = R   (iCF,jCF)
        invR_X1(iCF,jCF,iX1,iX2,iX3) = invR(iCF,jCF)

      END DO
      END DO

      IF( nDimsX .GT. 1 )THEN

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 2, GK, UK, R, invR )

        DO iCF = 1, nCF
        DO jCF = 1, nCF

          R_X2   (iCF,jCF,iX1,iX2,iX3) = R   (iCF,jCF)
          invR_X2(iCF,jCF,iX1,iX2,iX3) = invR(iCF,jCF)

        END DO
        END DO

      END IF

      IF( nDimsX .GT. 2 )THEN

        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 3, GK, UK, R, invR )

        DO iCF = 1, nCF
        DO jCF = 1, nCF

          R_X3   (iCF,jCF,iX1,iX2,iX3) = R   (iCF,jCF)
          invR_X3(iCF,jCF,iX1,iX2,iX3) = invR(iCF,jCF)

        END DO
        END DO

      END IF

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, D, LimitedCell, G_K, U_M, dX1, dX2, dX3, &
    !$ACC          R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3, &
    !$ACC          ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$ACC PRIVATE( a, b, c, dU, dU_C, SlopeDifference )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE

      DO iCF = 1, nCF

        a(iCF,1) = Zero
        b(iCF,1) = Zero
        c(iCF,1) = Zero

        DO jCF = 1, nCF

          a(iCF,1) &
            = a(iCF,1) + invR_X1(iCF,jCF,iX1,iX2,iX3) * U_M(2,jCF,iX1,iX2,iX3)

          b(iCF,1) &
            = BetaTVD * b(iCF,1) &
                + invR_X1(iCF,jCF,iX1,iX2,iX3) * ( U_M(1,jCF,iX1  ,iX2,iX3) &
                                                 - U_M(1,jCF,iX1-1,iX2,iX3) )

          c(iCF,1) &
            = BetaTVD * c(iCF,1) &
                + invR_X1(iCF,jCF,iX1,iX2,iX3) * ( U_M(1,jCF,iX1+1,iX2,iX3) &
                                                 - U_M(1,jCF,iX1  ,iX2,iX3) )

        END DO

      END DO

      IF( iX1 .EQ. iX_B0(1) .AND. ExcludeInnerGhostCell(1) )THEN

        DO iCF = 1, nCF

          b(iCF,1) = c(iCF,1)

        END DO

      END IF

      IF( iX1 .EQ. iX_E0(1) .AND. ExcludeOuterGhostCell(1) )THEN

        DO iCF = 1, nCF

          c(iCF,1) = b(iCF,1)

        END DO

      END IF

      DO iCF = 1, nCF

        dU_C(iCF) = MinModB( a(iCF,1), b(iCF,1), c(iCF,1), dX1(iX1), BetaTVB )

        dU(iCF,1) = Zero

      END DO

      DO iCF = 1, nCF
      DO jCF = 1, nCF

        dU(iCF,1) = dU(iCF,1) + R_X1(iCF,jCF,iX1,iX2,iX3) * dU_C(jCF)

      END DO
      END DO

      IF( nDimsX .GT. 1 )THEN

        DO iCF = 1, nCF

          a(iCF,2) = Zero
          b(iCF,2) = Zero
          c(iCF,2) = Zero

          DO jCF = 1, nCF

            a(iCF,2) &
              = a(iCF,2) + invR_X2(iCF,jCF,iX1,iX2,iX3) * U_M(2,jCF,iX1,iX2,iX3)

            b(iCF,2) &
              = BetaTVD * b(iCF,2) &
                  + invR_X2(iCF,jCF,iX1,iX2,iX3) * ( U_M(1,jCF,iX1,iX2  ,iX3) &
                                                   - U_M(1,jCF,iX1,iX2-1,iX3) )

            c(iCF,2) &
              = BetaTVD * c(iCF,2) &
                  + invR_X2(iCF,jCF,iX1,iX2,iX3) * ( U_M(1,jCF,iX1,iX2+1,iX3) &
                                                   - U_M(1,jCF,iX1,iX2  ,iX3) )

          END DO

        END DO

        IF( iX2 .EQ. iX_B0(2) .AND. ExcludeInnerGhostCell(2) )THEN

          DO iCF = 1, nCF

            b(iCF,2) = c(iCF,2)

          END DO

        END IF

        IF( iX2 .EQ. iX_E0(2) .AND. ExcludeOuterGhostCell(2) )THEN

          DO iCF = 1, nCF

            c(iCF,2) = b(iCF,2)

          END DO

        END IF

        DO iCF = 1, nCF

          dU_C(iCF) = MinModB( a(iCF,2), b(iCF,2), c(iCF,2), dX2(iX2), BetaTVB )

          dU(iCF,2) = Zero

        END DO

        DO iCF = 1, nCF
        DO jCF = 1, nCF

          dU(iCF,2) = dU(iCF,2) + R_X2(iCF,jCF,iX1,iX2,iX3) * dU_C(jCF)

        END DO
        END DO

      END IF

      IF( nDimsX .GT. 2 )THEN

        DO iCF = 1, nCF

          a(iCF,3) = Zero
          b(iCF,3) = Zero
          c(iCF,3) = Zero

          DO jCF = 1, nCF

            a(iCF,3) &
              = a(iCF,3) + invR_X3(iCF,jCF,iX1,iX2,iX3) * U_M(2,jCF,iX1,iX2,iX3)

            b(iCF,3) &
              = BetaTVD * b(iCF,3) &
                  + invR_X3(iCF,jCF,iX1,iX2,iX3) * ( U_M(1,jCF,iX1,iX2,iX3  ) &
                                                   - U_M(1,jCF,iX1,iX2,iX3-1) )

            c(iCF,3) &
              = BetaTVD * c(iCF,3) &
                  + invR_X3(iCF,jCF,iX1,iX2,iX3) * ( U_M(1,jCF,iX1,iX2,iX3+1) &
                                                   - U_M(1,jCF,iX1,iX2,iX3  ) )

          END DO

        END DO

        IF( iX3 .EQ. iX_B0(3) .AND. ExcludeInnerGhostCell(3) )THEN

          DO iCF = 1, nCF

            b(iCF,3) = c(iCF,3)

          END DO

        END IF

        IF( iX3 .EQ. iX_E0(3) .AND. ExcludeOuterGhostCell(3) )THEN

          DO iCF = 1, nCF

            c(iCF,3) = b(iCF,3)

          END DO

        END IF

        DO iCF = 1, nCF

          dU_C(iCF) = MinModB( a(iCF,3), b(iCF,3), c(iCF,3), dX3(iX3), BetaTVB )

          dU(iCF,3) = Zero

        END DO

        DO iCF = 1, nCF
        DO jCF = 1, nCF

          dU(iCF,3) = dU(iCF,3) + R_X3(iCF,jCF,iX1,iX2,iX3) * dU_C(jCF)

        END DO
        END DO

      END IF

    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( iX_B0, iX_E0, D, LimitedCell, G_K, U_M, dX1, dX2, dX3, &
    !$ACC          R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3, &
    !$ACC          ExcludeInnerGhostCell, ExcludeOuterGhostCell ) &
    !$ACC PRIVATE( a, b, c, dU, dU_C, SlopeDifference )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( D(1,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) CYCLE
      ! --- Compare Limited Slopes to Original Slopes ---

      DO iCF = 1, nCF

        SlopeDifference(iCF) = ABS( U_M(2,iCF,iX1,iX2,iX3) - dU(iCF,1) )

        IF( nDimsX .GT. 1 ) &
          SlopeDifference(iCF) &
            = MAX( SlopeDifference(iCF), &
                   ABS( U_M(3,iCF,iX1,iX2,iX3) - dU(iCF,2) ) )

        IF( nDimsX .GT. 2 ) &
          SlopeDifference(iCF) &
            = MAX( SlopeDifference(iCF), &
                   ABS( U_M(4,iCF,iX1,iX2,iX3) - dU(iCF,3) ) )

        ! --- Replace Slopes and Discard High-Order Components ---
        ! --- if Limited Slopes Deviate too Much from Original ---

        LimitedCell(iCF,iX1,iX2,iX3) = .FALSE.

        IF( SlopeDifference(iCF) &
              .GT. SlopeTolerance * ABS( U_M(1,iCF,iX1,iX2,iX3) ) )THEN

          DO iNX = 2, nDOFX

            U_M(iNX,iCF,iX1,iX2,iX3) = Zero

          END DO

          U_M(2,iCF,iX1,iX2,iX3) = dU(iCF,1)

          IF( nDimsX .GT. 1 ) U_M(3,iCF,iX1,iX2,iX3) = dU(iCF,2)

          IF( nDimsX .GT. 2 ) U_M(4,iCF,iX1,iX2,iX3) = dU(iCF,3)

          LimitedCell(iCF,iX1,iX2,iX3) = .TRUE.

        END IF

      END DO

    END DO
    END DO
    END DO

    ! --- Apply Conservative Correction ---

    CALL ApplyConservativeCorrection &
           ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    ! --- Map Modal to Nodal ---

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nDOFX, nCF_E, nDOFX, One, Pij_X, nDOFX, &
             U_M, nDOFX, Zero, U_N, nDOFX )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( iX_B1, iX_E1, U, U_N )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B1(3), iX_E1(3)
    DO iX2 = iX_B1(2), iX_E1(2)
    DO iX1 = iX_B1(1), iX_E1(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U(iNX,iX1,iX2,iX3,iCF) = U_N(iNX,iCF,iX1,iX2,iX3)

    END DO
    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, G, D, dX1, dX2, dX3, &
    !$OMP               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$OMP               SqrtGm, Vol, U_X, U_K, U_N, U_M, G_X, G_K, LimitedCell,&
    !$OMP               R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, G, D, dX1, dX2, dX3, &
    !$ACC               ExcludeInnerGhostCell, ExcludeOuterGhostCell, &
    !$ACC               SqrtGm, Vol, U_X, U_K, U_N, U_M, G_X, G_K, LimitedCell,&
    !$ACC               R_X1, invR_X1, R_X2, invR_X2, R_X3, invR_X3 )
#endif

    END ASSOCIATE

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_TVD_Characteristic


  SUBROUTINE InitializeSlopeLimiter_Euler_WENO

    INTEGER  :: iNodeX, iNodeX1, iNodeX2, iGridPt, iPoly, k, nP
    REAL(DP) :: OrthonormalBasis1D(1:nNodes,0:nNodes-1,0:nNodes-1)
    REAL(DP) :: eta

    IF( nDimsX .GT. 2 )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
         'WENO slope-limiter not implemented for nDimsX > 2. Stopping...'
      STOP

    END  IF

    IF( nNodes .GT. 3 )THEN

      WRITE(*,*)
      WRITE(*,'(A)') &
         'WENO slope-limiter not implemented for nNodes > 3. Stopping...'
      STOP

    END  IF

    IF( nNodes .EQ. 2 )THEN

      DO iNodeX = 1, nNodes

        eta = MeshX(1) % Nodes( iNodeX )

        ! --- 0th derivative ---
        OrthonormalBasis1D(iNodeX,0,0) &
          = One
        OrthonormalBasis1D(iNodeX,1,0) &
          = SQRT(12.0_DP) * eta

        ! --- 1st derivative ---
        OrthonormalBasis1D(iNodeX,0,1) &
          = Zero
        OrthonormalBasis1D(iNodeX,1,1) &
          = SQRT(12.0_DP)

      END DO

    ELSE IF( nNodes .EQ. 3 )THEN

      DO iNodeX = 1, nNodes

        eta = MeshX(1) % Nodes( iNodeX )

        ! --- 0th derivative ---
        OrthonormalBasis1D(iNodeX,0,0) &
          = One
        OrthonormalBasis1D(iNodeX,1,0) &
          = SQRT(12.0_DP) * eta
        OrthonormalBasis1D(iNodeX,2,0) &
          = SQRT(180.0_DP) * ( eta**2 - One / 12.0_DP )

        ! --- 1st derivative ---
        OrthonormalBasis1D(iNodeX,0,1) &
          = Zero
        OrthonormalBasis1D(iNodeX,1,1) &
          = SQRT(12.0_DP)
        OrthonormalBasis1D(iNodeX,2,1) &
          = SQRT(180.0_DP) * Two * eta

        ! --- 2nd derivative ---
        OrthonormalBasis1D(iNodeX,0,2) &
          = Zero
        OrthonormalBasis1D(iNodeX,1,2) &
          = Zero
        OrthonormalBasis1D(iNodeX,2,2) &
          = SQRT(180.0_DP) * Two

      END DO

    END IF

    k  = nNodes - 1
    nP = UpperLimit_pSpace( k )

    ALLOCATE( OrthonormalBasis(1:nDOFX,0:nDOFX-1,0:nP) )

    IF     ( nDimsX .EQ. 1 )THEN

      OrthonormalBasis = OrthonormalBasis1D

    ELSE IF( nDimsX .EQ. 2 )THEN

      OrthonormalBasis = Zero

      IF( k .EQ. 1 )THEN

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)

          ! --- 0th derivative ---
          OrthonormalBasis(iNodeX,0,0) = OrthonormalBasis1D(iNodeX1,0,0)
          OrthonormalBasis(iNodeX,1,0) = OrthonormalBasis1D(iNodeX1,1,0)
          OrthonormalBasis(iNodeX,2,0) = OrthonormalBasis1D(iNodeX2,1,0)

          ! --- 1st eta1 derivative ---
          OrthonormalBasis(iNodeX,0,1) = OrthonormalBasis1D(iNodeX1,0,1)
          OrthonormalBasis(iNodeX,1,1) = OrthonormalBasis1D(iNodeX1,1,1)
          OrthonormalBasis(iNodeX,2,1) = Zero

          ! --- 1st eta2 derivative ---
          OrthonormalBasis(iNodeX,0,2) = OrthonormalBasis1D(iNodeX2,0,1)
          OrthonormalBasis(iNodeX,1,2) = Zero
          OrthonormalBasis(iNodeX,2,2) = OrthonormalBasis1D(iNodeX2,1,1)

        END DO

      ELSE IF( k .EQ. 2 )THEN

        DO iNodeX = 1, nDOFX

          iNodeX1 = NodeNumberTableX(1,iNodeX)
          iNodeX2 = NodeNumberTableX(2,iNodeX)

          ! --- 0th derivative ---
          OrthonormalBasis(iNodeX,0,0) = OrthonormalBasis1D(iNodeX1,0,0)
          OrthonormalBasis(iNodeX,1,0) = OrthonormalBasis1D(iNodeX1,1,0)
          OrthonormalBasis(iNodeX,2,0) = OrthonormalBasis1D(iNodeX2,1,0)
          OrthonormalBasis(iNodeX,3,0) = OrthonormalBasis1D(iNodeX1,2,0)
          OrthonormalBasis(iNodeX,4,0) = OrthonormalBasis1D(iNodeX1,1,0) &
                                           * OrthonormalBasis1D(iNodeX2,1,0)
          OrthonormalBasis(iNodeX,5,0) = OrthonormalBasis1D(iNodeX2,2,0)

          ! --- 1st eta1 derivative ---
          OrthonormalBasis(iNodeX,0,1) = OrthonormalBasis1D(iNodeX1,0,1)
          OrthonormalBasis(iNodeX,1,1) = OrthonormalBasis1D(iNodeX1,1,1)
          OrthonormalBasis(iNodeX,2,1) = Zero
          OrthonormalBasis(iNodeX,3,1) = OrthonormalBasis1D(iNodeX1,2,1)
          OrthonormalBasis(iNodeX,4,1) = OrthonormalBasis1D(iNodeX2,1,0) &
                                           * OrthonormalBasis1D(iNodeX1,1,1)

          OrthonormalBasis(iNodeX,5,1) = Zero

          ! --- 1st eta2 derivative ---
          OrthonormalBasis(iNodeX,0,2) = OrthonormalBasis1D(iNodeX2,0,1)
          OrthonormalBasis(iNodeX,1,2) = Zero
          OrthonormalBasis(iNodeX,2,2) = OrthonormalBasis1D(iNodeX2,1,1)
          OrthonormalBasis(iNodeX,3,2) = Zero
          OrthonormalBasis(iNodeX,4,2) = OrthonormalBasis1D(iNodeX1,1,0) &
                                           * OrthonormalBasis1D(iNodeX2,1,1)
          OrthonormalBasis(iNodeX,5,2) = OrthonormalBasis1D(iNodeX2,2,1)

          ! --- 2nd eta1 derivative ---
          OrthonormalBasis(iNodeX,0,3) = OrthonormalBasis1D(iNodeX1,0,2)
          OrthonormalBasis(iNodeX,1,3) = OrthonormalBasis1D(iNodeX1,1,2)
          OrthonormalBasis(iNodeX,2,3) = Zero
          OrthonormalBasis(iNodeX,3,3) = OrthonormalBasis1D(iNodeX1,2,2)
          OrthonormalBasis(iNodeX,4,3) = OrthonormalBasis1D(iNodeX2,1,0) &
                                           * OrthonormalBasis1D(iNodeX1,1,2)
          OrthonormalBasis(iNodeX,5,3) = Zero

          ! --- 1st mixed derivative ---
          OrthonormalBasis(iNodeX,0,4) = Zero
          OrthonormalBasis(iNodeX,1,4) = Zero
          OrthonormalBasis(iNodeX,2,4) = Zero
          OrthonormalBasis(iNodeX,3,4) = Zero
          OrthonormalBasis(iNodeX,4,4) = OrthonormalBasis1D(iNodeX1,1,1) &
                                           * OrthonormalBasis1D(iNodeX2,1,1)
          OrthonormalBasis(iNodeX,5,4) = Zero

          ! --- 2nd eta2 derivative ---
          OrthonormalBasis(iNodeX,0,5) = OrthonormalBasis1D(iNodeX2,0,2)
          OrthonormalBasis(iNodeX,1,5) = Zero
          OrthonormalBasis(iNodeX,2,5) = OrthonormalBasis1D(iNodeX2,1,2)
          OrthonormalBasis(iNodeX,3,5) = Zero
          OrthonormalBasis(iNodeX,4,5) = OrthonormalBasis1D(iNodeX1,1,0) &
                                           * OrthonormalBasis1D(iNodeX2,1,2)
          OrthonormalBasis(iNodeX,5,5) = OrthonormalBasis1D(iNodeX2,2,2)

        END DO

      END IF ! Order of accuracy

    END IF ! Spatial dimensions

    ALLOCATE( VandermondeMatrix(0:nDOFX-1,1:nDOFX) )

    DO iPoly = 0, nDOFX-1

      DO iGridPt = 1, nDOFX

        VandermondeMatrix(iPoly,iGridPt) &
          = WeightsX_q(iGridPt) * OrthonormalBasis(iGridPt,iPoly,0)

      END DO

    END DO

  END SUBROUTINE InitializeSlopeLimiter_Euler_WENO


  SUBROUTINE FinalizeSlopeLimiter_Euler_Relativistic_IDEAL

    DEALLOCATE( LegendreX )

    IF( UseTroubledCellIndicator ) &
      CALL FinalizeTroubledCellIndicator_Euler

    IF( TRIM( SlopeLimiterMethod ) .EQ. 'WENO' )THEN

      DEALLOCATE( VandermondeMatrix )
      DEALLOCATE( OrthonormalBasis )

    END IF

  END SUBROUTINE FinalizeSlopeLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO_ComponentWise &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    LOGICAL  :: SuppressBC
    INTEGER  :: iX1, iX2, iX3, iGridPt, iCF, ell
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: U_M(0:nDOFX-1,0:2*nDimsX,nCF)
    REAL(DP) :: UU(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),nCF)
    LOGICAL  :: LimitedCell(iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

    ! --- WENO Limiter ---

    INTEGER  :: nPspace, k
    REAL(DP) :: q    (nDOFX,0:nNodes-1  ,nCF)
    REAL(DP) :: LinearWeights       (2*(nNodes-1))

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    k = nNodes - 1

    ! --- Get linear weights ---

    DO ell = 1, 2*k

      IF( MOD( ell, 2 ) .NE. 0 )THEN

        LinearWeights(ell) = 0.01_DP

      ELSE

        LinearWeights(ell) = 0.99_DP

      END IF

    END DO

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    LimitedCell = .FALSE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ALL( D(:,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) ) CYCLE

      LimitedCell(iX1,iX2,iX3) = .TRUE.

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      DO iCF = 1, nCF

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1,iX2,iX3,iCF), &
                 U_M(:,0,iCF) )

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1-1,iX2,iX3,iCF), &
                 U_M(:,1,iCF) )

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1+1,iX2,iX3,iCF), &
                 U_M(:,2,iCF) )

        IF( nDimsX .GT. 1 )THEN

          CALL MapNodalToModal_Fluid_WENO &
                 ( VandermondeMatrix, &
                   U(:,iX1,iX2-1,iX3,iCF), &
                   U_M(:,3,iCF) )

          CALL MapNodalToModal_Fluid_WENO &
                 ( VandermondeMatrix, &
                   U(:,iX1,iX2+1,iX3,iCF), &
                   U_M(:,4,iCF) )

        END IF

        ! --- Step 1.1 ---

        DO ell = 0, k ! Loop over basis polynomials

          nPspace = UpperLimit_pSpace( ell )

          DO iGridPt = 1, nDOFX

            q(iGridPt,ell,iCF) &
              = SUM( U_M(0:nPspace,0,iCF) &
                       * OrthonormalBasis(iGridPt,0:nPspace,0) )

          END DO

        END DO

        ! --- Steps 1.2-1.5 ---

        CALL ApplySlopeLimiter_WENO_Scalar &
               ( k, LinearWeights, q(:,:,iCF), U_M(:,:,iCF), &
                 UU(:,iX1,iX2,iX3,iCF) )

      END DO ! End of loop over fields

    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( LimitedCell(iX1,iX2,iX3) )THEN

        DO iCF = 1, nCF

          U(:,iX1,iX2,iX3,iCF) = UU(:,iX1,iX2,iX3,iCF)

        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO_ComponentWise


  SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO_Characteristic &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U, D, SuppressBC_Option, iApplyBC_Option )

    INTEGER,  INTENT(in)           :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)           :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout)        :: &
      D(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    LOGICAL,  INTENT(in), OPTIONAL :: &
      SuppressBC_Option
    INTEGER,  INTENT(in), OPTIONAL :: &
      iApplyBC_Option(3)

    LOGICAL  :: SuppressBC
    INTEGER  :: iX1, iX2, iX3, iNodeX, iGridPt, iCF, iCell, ell
    REAL(DP) :: dX1, dX2, dX3
    REAL(DP) :: U_M   (0:nDOFX-1,0:2*nDimsX,nCF)
    REAL(DP) :: U_M_X1(0:nDOFX-1,0:2*nDimsX,nCF)
    REAL(DP) :: U_M_X2(0:nDOFX-1,0:2*nDimsX,nCF)
    REAL(DP) :: UU(nDOFX,iX_B0(1):iX_E0(1), &
                         iX_B0(2):iX_E0(2), &
                         iX_B0(3):iX_E0(3),nCF)
    LOGICAL  :: LimitedCell(iX_B0(1):iX_E0(1), &
                            iX_B0(2):iX_E0(2), &
                            iX_B0(3):iX_E0(3))

    REAL(DP) :: q  (nDOFX,0:nNodes-1,nCF)
    REAL(DP) :: qX1(nDOFX,0:nNodes-1,nCF)
    REAL(DP) :: qX2(nDOFX,0:nNodes-1,nCF)
    REAL(DP) :: LinearWeights(2*(nNodes-1))
    INTEGER  :: nPspace, k

    ! --- Characteristic limiting ---
    REAL(DP) :: R_X1(nCF,nCF), invR_X1(nCF,nCF)
    REAL(DP) :: R_X2(nCF,nCF), invR_X2(nCF,nCF)
    REAL(DP) :: U_X1(nDOFX,nCF), U_X2(nDOFX,nCF)
    REAL(DP) :: U_K(nCF), G_K(nGF)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UseSlopeLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_SlopeLimiter )

    k = nNodes - 1

    ! --- Get linear weights ---
    DO iNodeX = 1, 2*k

      IF( MOD( iNodeX, 2 ) .NE. 0 )THEN

        LinearWeights(iNodeX) = 0.01_DP

      ELSE

        LinearWeights(iNodeX) = 0.99_DP

      END IF

    END DO

    SuppressBC = .FALSE.
    IF( PRESENT( SuppressBC_Option ) ) &
      SuppressBC = SuppressBC_Option

    IF( .NOT. SuppressBC ) &
      CALL ApplyBoundaryConditions_Euler &
             ( iX_B0, iX_E0, iX_B1, iX_E1, U )

    CALL DetectTroubledCells_Euler &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, D )

    LimitedCell = .FALSE.

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( ALL( D(:,iX1,iX2,iX3,iDF_TCI) .LT. LimiterThreshold ) ) CYCLE

      LimitedCell(iX1,iX2,iX3) = .TRUE.

      dX1 = MeshX(1) % Width(iX1)
      dX2 = MeshX(2) % Width(iX2)
      dX3 = MeshX(3) % Width(iX3)

      ! --- Compute cell averages (ignore geometry) ---

      DO iCF = 1, nCF

        U_K(iCF) = SUM( WeightsX_q * U(:,iX1,iX2,iX3,iCF) )

      END DO

      ! --- Step 2.1.1 ---

      G_K(iGF_Gm_dd_11) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_Gm_dd_11) )
      G_K(iGF_Gm_dd_22) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_Gm_dd_22) )
      G_K(iGF_Gm_dd_33) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_Gm_dd_33) )
      G_K(iGF_SqrtGm) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_SqrtGm  ) )
      G_K(iGF_Alpha) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_Alpha   ) )
      G_K(iGF_Beta_1) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_Beta_1  ) )
      G_K(iGF_Beta_2) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_Beta_2  ) )
      G_K(iGF_Beta_3) &
        = DOT_PRODUCT( WeightsX_q, G(:,iX1,iX2,iX3,iGF_Beta_3  ) )

      CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
             ( 1, G_K, U_K, R_X1, invR_X1 )

      IF( nDimsX .GT. 1 ) &
        CALL ComputeCharacteristicDecomposition_Euler_Relativistic_IDEAL &
               ( 2, G_K, U_K, R_X2, invR_X2 )

      DO iCF = 1, nCF

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1,iX2,iX3,iCF), &
                 U_M(:,0,iCF) )

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1-1,iX2,iX3,iCF), &
                 U_M(:,1,iCF) )

        CALL MapNodalToModal_Fluid_WENO &
               ( VandermondeMatrix, &
                 U(:,iX1+1,iX2,iX3,iCF), &
                 U_M(:,2,iCF) )

        IF( nDimsX .GT. 1 )THEN

          CALL MapNodalToModal_Fluid_WENO &
                 ( VandermondeMatrix, &
                   U(:,iX1,iX2-1,iX3,iCF), &
                   U_M(:,3,iCF) )

          CALL MapNodalToModal_Fluid_WENO &
                 ( VandermondeMatrix, &
                   U(:,iX1,iX2+1,iX3,iCF), &
                   U_M(:,4,iCF) )

        END IF

      END DO

      ! --- Step 2.1.0 (compute "associated polynomial vectors") ---

      DO iCF = 1, nCF

        DO ell = 0, k

          nPspace = UpperLimit_pSpace( ell )

          DO iGridPt = 1, nDOFX

            q(iGridPt,ell,iCF) &
              = SUM( U_M(0:nPspace,0,iCF) &
                       * OrthonormalBasis(iGridPt,0:nPspace,0) )

          END DO

        END DO

      END DO

      ! --- Step 2.1.2 ---

      DO ell = 0, k

        DO iGridPt = 1, nDOFX

          qX1(iGridPt,ell,:) &
            = MATMUL( invR_X1, q(iGridPt,ell,:) )

        END DO

      END DO

      IF( nDimsX .GT. 1 )THEN

        DO ell = 0, k

          DO iGridPt = 1, nDOFX

            qX2(iGridPt,ell,:) &
              = MATMUL( invR_X2, q(iGridPt,ell,:) )

          END DO

        END DO

      END IF

      DO iCell = 0, 2*nDimsX

        DO iNodeX = 0, nDOFX-1

          U_M_X1(iNodeX,iCell,:) = MATMUL( invR_X1, U_M(iNodeX,iCell,:) )

        END DO

      END DO

      IF( nDimsX .GT. 1 )THEN

        DO iCell = 0, 2*nDimsX

          DO iNodeX = 0, nDOFX-1

            U_M_X2(iNodeX,iCell,:) = MATMUL( invR_X2, U_M(iNodeX,iCell,:) )

          END DO

        END DO

      END IF

      ! --- Step 2.1.3 ---

      DO iCF = 1, nCF

        CALL ApplySlopeLimiter_WENO_Scalar &
               ( k, LinearWeights, qX1(:,:,iCF), U_M_X1(:,:,iCF), &
                 U_X1(:,iCF) )

      END DO

      IF( nDimsX .GT. 1 )THEN

        DO iCF = 1, nCF

          CALL ApplySlopeLimiter_WENO_Scalar &
               ( k, LinearWeights, qX2(:,:,iCF), U_M_X2(:,:,iCF), &
                 U_X2(:,iCF) )

        END DO

      END IF

      ! --- Step 2.1.4 (X1) ---

      DO iNodeX = 1, nDOFX

        U_X1(iNodeX,:) = MATMUL( R_X1, U_X1(iNodeX,:) )

      END DO

      ! --- Step 2.2 (X1) ---

      UU(:,iX1,iX2,iX3,:) = U_X1 / DBLE( nDimsX )

      ! --- Step 2.1.4 (X2) ---

      IF( nDimsX .GT. 1 )THEN

        DO iNodeX = 1, nDOFX

          U_X2(iNodeX,:) = MATMUL( R_X2, U_X2(iNodeX,:) )

        END DO

        ! --- Step 2.2 (X2)---

        UU(:,iX1,iX2,iX3,:) = UU(:,iX1,iX2,iX3,:) + U_X2 / DBLE( nDimsX )

      END IF

    END DO
    END DO
    END DO

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( LimitedCell(iX1,iX2,iX3) )THEN

        DO iCF = 1, nCF

          U(:,iX1,iX2,iX3,iCF) = UU(:,iX1,iX2,iX3,iCF)

        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_SlopeLimiter )

  END SUBROUTINE ApplySlopeLimiter_Euler_Relativistic_IDEAL_WENO_Characteristic


  SUBROUTINE ApplySlopeLimiter_WENO_Scalar( k, LinearWeights, q, U_M, UU )

    INTEGER,  INTENT(in)  :: k
    REAL(DP), INTENT(in)  :: LinearWeights(2*k), &
                             q(nDOFX,0:k), &
                             U_M(0:nDOFX-1,0:2*nDimsX)
    REAL(DP), INTENT(out) :: UU(nDOFX)

    INTEGER  :: iGridPt
    REAL(DP) :: SmoothnessIndicators(2*k)
    REAL(DP) :: NonLinearWeights    (2*k)
    REAL(DP) :: pCoeffs(2*k,nNodes)
    REAL(DP) :: p(nDOFX,2*k)

    ! --- 2nd Order ---

    ! --- Step 1.2 ---

    CALL ComputeCoefficients_Order2( LinearWeights, pCoeffs )

    DO iGridPt = 1, nDOFX

      p(iGridPt,1:2) = MATMUL( pCoeffs(1:2,1:2), q(iGridPt,0:1) )

    END DO

    ! --- Step 1.3 ---

    ! --- Hard-code smoothness indicator beta_01 ---

    IF     ( nDimsX .EQ. 1 )THEN

      SmoothnessIndicators(1) &
        = OrthonormalBasis(1,1,1)**2 &
            * MIN( U_M(1,1)**2, &
                   U_M(1,2)**2 )

    ELSE IF( nDimsX .EQ. 2 )THEN

      SmoothnessIndicators(1) &
        = OrthonormalBasis(1,1,1)**2 &
            * MIN( U_M(1,1)**2 + U_M(2,1)**2, &
                   U_M(1,2)**2 + U_M(2,2)**2, &
                   U_M(1,3)**2 + U_M(2,3)**2, &
                   U_M(1,4)**2 + U_M(2,4)**2 )

    END IF

    CALL ComputeSmoothnessIndicator_Order2 &
           ( LinearWeights(2), &
             OrthonormalBasis(1,1,1), &
             U_M(1:nDimsX,0), &
             SmoothnessIndicators(2) )

    ! --- Step 1.4 ---

    CALL ComputeNonLinearWeights &
           ( SmoothnessIndicators(1:2), &
             LinearWeights       (1:2), &
             NonLinearWeights    (1:2) )

    ! --- 3rd order ---

    IF( nNodes .GT. 2 )THEN

      ! --- Step 1.2 ---

      CALL ComputeCoefficients_Order3 &
             ( LinearWeights, NonLinearWeights, pCoeffs )

      DO iGridPt = 1, nDOFX

        p(iGridPt,3:4) &
          = MATMUL( pCoeffs(3:4,1:3), q(iGridPt,0:2) )

      END DO

      ! --- Step 1.3 ---

      CALL ComputeSmoothnessIndicators_Order3 &
             ( LinearWeights(2), NonLinearWeights(2), &
               U_M(:,0), pCoeffs(4,:), SmoothnessIndicators )

      ! --- Step 1.4 ---

      CALL ComputeNonLinearWeights &
             ( SmoothnessIndicators(3:4), &
               LinearWeights       (3:4), &
               NonLinearWeights    (3:4) )

    END IF

    ! --- Step 1.5 ---

    DO iGridPt = 1, nDOFX

      UU(iGridPt) &
        = SUM( NonLinearWeights(2*k-1:2*k) &
                 * p(iGridPt,2*k-1:2*k) )

    END DO

  END SUBROUTINE ApplySlopeLimiter_WENO_Scalar


  INTEGER FUNCTION UpperLimit_pSpace( ell )

    INTEGER, INTENT(in) :: ell

    IF( nDimsX .EQ. 1 )THEN

      UpperLimit_pSpace = ell

    ELSE IF( nDimsX .EQ. 2 )THEN

      UpperLimit_pSpace = ( ell + 1 ) * ( ell + 2 ) / 2 - 1

    END IF

    RETURN
  END FUNCTION UpperLimit_pSpace


  SUBROUTINE ComputeSmoothnessIndicator_Order2 &
    ( LinearWeight_11, Coeff_v1, U_M, SmoothnessIndicator )

    REAL(DP), INTENT(in)  :: LinearWeight_11, Coeff_v1, U_M(:)
    REAL(DP), INTENT(out) :: SmoothnessIndicator

    ! --- beta_11 ---

    SmoothnessIndicator = ( U_M(1) / LinearWeight_11 * Coeff_v1 )**2

    IF( nDimsX .GT. 1 ) &
      SmoothnessIndicator &
        = SmoothnessIndicator + ( U_M(2) / LinearWeight_11 * Coeff_v1 )**2

  END SUBROUTINE ComputeSmoothnessIndicator_Order2


  SUBROUTINE ComputeSmoothnessIndicators_Order3 &
    ( LinearWeight_11, NonLinearWeight_11, U_M, &
      pCoeffs, SmoothnessIndicators )

    REAL(DP), INTENT(in)  :: LinearWeight_11, NonLinearWeight_11
    REAL(DP), INTENT(in)  :: U_M(0:nDOFX-1), pCoeffs(:)
    REAL(DP), INTENT(out) :: SmoothnessIndicators(:)

    REAL(DP) :: Coeff_v1
    REAL(DP) :: dq_X1   (nDOFX,0:2), &
                dq_X2   (nDOFX,0:2), &
                ddq_X1  (nDOFX,0:2), &
                ddq_X2  (nDOFX,0:2), &
                ddq_X1X2(nDOFX,0:2)
    REAL(DP) :: dP_X1   (nDOFX), &
                dP_X2   (nDOFX), &
                ddP_X1  (nDOFX), &
                ddP_X2  (nDOFX), &
                ddP_X1X2(nDOFX)
    INTEGER  :: iGridPt

    ! --- beta_12 ---

    Coeff_v1 = OrthonormalBasis(1,1,1)

    SmoothnessIndicators(3) &
      = ( Coeff_v1 * NonLinearWeight_11 / LinearWeight_11 * U_M(1) )**2

    IF( nDimsX .GT. 1 ) &
      SmoothnessIndicators(3) &
        = SmoothnessIndicators(3) &
            + ( Coeff_v1 * NonLinearWeight_11 / LinearWeight_11 * U_M(2) )**2

    ! --- beta_22 ---

    IF( nDimsX .EQ. 1 )THEN

      DO iGridPt = 1, nDOFX

        dq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,1) )
        dq_X1(iGridPt,1) = SUM( U_M(0:1) * OrthonormalBasis(iGridPt,0:1,1) )
        dq_X1(iGridPt,2) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,1) )
        dp_X1 (iGridPt)  = DOT_PRODUCT( pCoeffs, dq_X1(iGridPt,:) )

        ddq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,2) )
        ddq_X1(iGridPt,1) = SUM( U_M(0:1) * OrthonormalBasis(iGridPt,0:1,2) )
        ddq_X1(iGridPt,2) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,2) )
        ddp_X1(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X1(iGridPt,:) )

      END DO

      SmoothnessIndicators(4) &
        = SUM( WeightsX_q * ( dp_X1**2 + ddp_X1**2 ) )

    ELSE IF( nDimsX .GT. 1 )THEN

      DO iGridPt = 1, nDOFX

        dq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,1) )
        dq_X1(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,1) )
        dq_X1(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,1) )
        dp_X1(iGridPt)   = DOT_PRODUCT( pCoeffs, dq_X1(iGridPt,:) )

        dq_X2(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,2) )
        dq_X2(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,2) )
        dq_X2(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,2) )
        dp_X2(iGridPt)   = DOT_PRODUCT( pCoeffs, dq_X2(iGridPt,:) )

        ddq_X1(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,3) )
        ddq_X1(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,3) )
        ddq_X1(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,3) )
        ddp_X1(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X1(iGridPt,:) )

        ddq_X1X2(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,4) )
        ddq_X1X2(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,4) )
        ddq_X1X2(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,4) )
        ddp_X1X2(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X1X2(iGridPt,:) )

        ddq_X2(iGridPt,0) = SUM( U_M(0:0) * OrthonormalBasis(iGridPt,0:0,5) )
        ddq_X2(iGridPt,1) = SUM( U_M(0:2) * OrthonormalBasis(iGridPt,0:2,5) )
        ddq_X2(iGridPt,2) = SUM( U_M(0:5) * OrthonormalBasis(iGridPt,0:5,5) )
        ddp_X2(iGridPt)   = DOT_PRODUCT( pCoeffs, ddq_X2(iGridPt,:) )

      END DO

      SmoothnessIndicators(4) &
        = SUM( WeightsX_q &
                 * ( dp_X1**2 + ddp_X1**2 + dp_X2**2 &
                       + ddp_X2**2 + ddp_X1X2**2 ) )

    END IF

  END SUBROUTINE ComputeSmoothnessIndicators_Order3


  SUBROUTINE ComputeCoefficients_Order2 &
    ( LinearWeights, PolynomialCoefficients )

    REAL(DP), INTENT(in)  :: LinearWeights(:)
    REAL(DP), INTENT(out) :: PolynomialCoefficients(:,:)

    PolynomialCoefficients(1,1) = One

    PolynomialCoefficients(1,2) = Zero

    PolynomialCoefficients(2,1) = -LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(2,2) = One / LinearWeights(2)

  END SUBROUTINE ComputeCoefficients_Order2


  SUBROUTINE ComputeCoefficients_Order3 &
    ( LinearWeights, NonLinearWeights, PolynomialCoefficients )

    REAL(DP), INTENT(in)  :: LinearWeights(:), NonLinearWeights(:)
    REAL(DP), INTENT(out) :: PolynomialCoefficients(:,:)

    PolynomialCoefficients(1,1) &
      = One

    PolynomialCoefficients(1,2) &
      = Zero

    PolynomialCoefficients(1,3) &
      = Zero

    PolynomialCoefficients(2,1) &
      = -LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(2,2) &
      = One / LinearWeights(2)

    PolynomialCoefficients(2,3) &
      = Zero

    PolynomialCoefficients(3,1) &
      = NonLinearWeights(1) &
          - NonLinearWeights(2) * LinearWeights(1) / LinearWeights(2)

    PolynomialCoefficients(3,2) &
      = NonLinearWeights(2) / LinearWeights(2)

    PolynomialCoefficients(3,3) &
      = Zero

    PolynomialCoefficients(4,1) &
      = -LinearWeights(3) / LinearWeights(4) * PolynomialCoefficients(3,1)

    PolynomialCoefficients(4,2) &
      = -LinearWeights(3) / LinearWeights(4) * PolynomialCoefficients(3,2)

    PolynomialCoefficients(4,3) &
      = One / LinearWeights(4)

  END SUBROUTINE ComputeCoefficients_Order3


  SUBROUTINE ComputeNonLinearWeights &
    ( SmoothnessIndicators, LinearWeights, NonLinearWeights )

    REAL(DP), INTENT(in)  :: SmoothnessIndicators(2), LinearWeights(2)
    REAL(DP), INTENT(out) :: NonLinearWeights(2)

    REAL(DP) :: Tau, Normalization
    REAL(DP), PARAMETER :: EPS = 1.0d-10

    Tau = ( SmoothnessIndicators(2) - SmoothnessIndicators(1) )**2

    Normalization &
      = SUM( LinearWeights * ( One + Tau / ( EPS + SmoothnessIndicators ) ) )

    NonLinearWeights(1) &
      = LinearWeights(1) * ( One + Tau / ( EPS + SmoothnessIndicators(1) ) ) &
          / Normalization

    NonLinearWeights(2) &
      = LinearWeights(2) * ( One + Tau / ( EPS + SmoothnessIndicators(2) ) ) &
          / Normalization

  END SUBROUTINE ComputeNonLinearWeights


  SUBROUTINE MapNodalToModal_Fluid_WENO( VandermondeMatrix, uN, uM )

    REAL(DP), INTENT(in)  :: VandermondeMatrix(0:nDOFX-1,nDOFX)
    REAL(DP), INTENT(in)  :: uN(nDOFX)
    REAL(DP), INTENT(out) :: uM(0:nDOFX-1)

    uM = MATMUL( VandermondeMatrix, uN )

  END SUBROUTINE MapNodalToModal_Fluid_WENO


  SUBROUTINE ApplyConservativeCorrection &
    ( iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, LimitedCell, U_M )

    INTEGER, INTENT(in)     :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      SqrtGm(1:nDOFX,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      Vol   (        iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3)), &
      U_K   (1:nCF  ,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    LOGICAL, INTENT(in)     :: &
      LimitedCell(1:nCF,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3))
    REAL(DP), INTENT(inout) :: &
      U_M(1:nDOFX,1:nCF,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3))

    INTEGER  :: iNX, iX1, iX2, iX3, iCF, iDimX
    REAL(DP) :: Correction, Term

    IF( .NOT. UseConservativeCorrection ) RETURN

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, &
    !$OMP          LimitedCell, U_M )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, &
    !$ACC          LimitedCell, U_M )
#endif

    ! --- Applies a correction to the 0-th order ---
    ! --- mode to maintain the cell average.     ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( iX_B0, iX_E0, LimitedCell, WeightsX_q, LegendreX, SqrtGm, &
    !$ACC          U_M, Vol, U_K )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      IF( LimitedCell(iCF,iX1,iX2,iX3) )THEN

        Correction = Zero

        DO iDimX = 1, nDimsX

          Term = Zero

          DO iNX = 1, nDOFX

            Term &
              = Term &
                  + WeightsX_q(iNX) &
                      * LegendreX(iNX,iDimX+1) * SqrtGm(iNX,iX1,iX2,iX3)

          END DO

          Correction &
            = Correction &
                + U_M(iDimX+1,iCF,iX1,iX2,iX3) * Term / Vol(iX1,iX2,iX3)

        END DO

        U_M(1,iCF,iX1,iX2,iX3) = U_K(iCF,iX1,iX2,iX3) - Correction

      END IF

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U_M ) &
    !$OMP MAP( release: iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, &
    !$OMP               LimitedCell )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U_M ) &
    !$ACC DELETE(       iX_B0, iX_E0, iX_B1, iX_E1, SqrtGm, Vol, U_K, &
    !$ACC               LimitedCell )
#endif

  END SUBROUTINE ApplyConservativeCorrection


END MODULE Euler_SlopeLimiterModule_Relativistic_IDEAL
