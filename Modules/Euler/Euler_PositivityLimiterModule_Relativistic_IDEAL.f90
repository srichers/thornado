#ifdef THORNADO_DEBUG
#define THORNADO_DEBUG_EULER_POSITIVITYLIMITER
#endif

!> Initialize, finalize, and apply positivity limiter from
!> Qin et al., (2016), JCP, 315, 323
MODULE Euler_PositivityLimiterModule_Relativistic_IDEAL

  USE KindModule, ONLY: &
    DP, &
    Zero, &
    Half, &
    One, &
    SqrtTiny
  USE ProgramHeaderModule, ONLY: &
    nNodesX, &
    nDOFX
  USE LinearAlgebraModule, ONLY: &
    MatrixVectorMultiply, &
    MatrixMatrixMultiply
  USE ReferenceElementModuleX, ONLY: &
    nDOFX_X1, &
    nDOFX_X2, &
    nDOFX_X3, &
    WeightsX_q
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    LX_X1_Dn, &
    LX_X1_Up, &
    LX_X2_Dn, &
    LX_X2_Up, &
    LX_X3_Dn, &
    LX_X3_Up
  USE GeometryFieldsModule, ONLY: &
    nGF, &
    iGF_h_1, &
    iGF_h_2, &
    iGF_h_3, &
    iGF_Gm_dd_11, &
    iGF_Gm_dd_22, &
    iGF_Gm_dd_33, &
    iGF_SqrtGm
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX_FromScaleFactors
  USE FluidFieldsModule, ONLY: &
    nCF, &
    iCF_D, &
    iCF_S1, &
    iCF_S2, &
    iCF_S3, &
    iCF_E
  USE TimersModule_Euler, ONLY: &
    TimersStart_Euler, &
    TimersStop_Euler, &
    Timer_Euler_PositivityLimiter, &
    Timer_Euler_PL_LimitCells, &
    Timer_Euler_PL_CopyIn, &
    Timer_Euler_PL_CopyOut, &
    Timer_Euler_PL_Permute, &
    Timer_Euler_PL_Integrate, &
    Timer_Euler_PL_ErrorCheck
  USE Euler_ErrorModule, ONLY: &
    DescribeError_Euler

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC :: InitializePositivityLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: FinalizePositivityLimiter_Euler_Relativistic_IDEAL
  PUBLIC :: ApplyPositivityLimiter_Euler_Relativistic_IDEAL

  LOGICAL               :: UsePositivityLimiter
  INTEGER, PARAMETER    :: nPS = 7  ! Number of Positive Point Sets
  INTEGER               :: nPP(nPS) ! Number of Positive Points Per Set
  INTEGER               :: nPT      ! Total number of Positive Points
  INTEGER               :: nCF_K
  REAL(DP)              :: Min_1, Min_2
  REAL(DP), ALLOCATABLE :: U_PP(:,:), G_PP(:,:), L_X(:,:)


CONTAINS


  SUBROUTINE InitializePositivityLimiter_Euler_Relativistic_IDEAL &
    ( UsePositivityLimiter_Option, Verbose_Option, Min_1_Option, Min_2_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: UsePositivityLimiter_Option, &
                                      Verbose_Option
    REAL(DP), INTENT(in), OPTIONAL :: Min_1_Option, Min_2_Option

    INTEGER :: iDim, iNX, iOS
    LOGICAL :: Verbose

    UsePositivityLimiter = .TRUE.
    IF( PRESENT( UsePositivityLimiter_Option ) ) &
      UsePositivityLimiter = UsePositivityLimiter_Option

    Verbose = .TRUE.
    IF( PRESENT( Verbose_Option ) ) &
      Verbose = Verbose_Option

    Min_1 = - HUGE( One )
    IF( PRESENT( Min_1_Option ) ) &
      Min_1 = Min_1_Option

    Min_2 = - HUGE( One )
    IF( PRESENT( Min_2_Option ) ) &
      Min_2 = Min_2_Option

    IF( Verbose )THEN
      WRITE(*,*)
      WRITE(*,'(A)') &
        '    INFO: Positivity Limiter (Euler, Relativistic, IDEAL)'
      WRITE(*,'(A)') &
        '    -----------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A6,A,L1)') &
        '', 'UsePositivityLimiter: ', UsePositivityLimiter
      WRITE(*,*)
      WRITE(*,'(A6,A,ES12.4E3)') &
        '', 'Min_1 = ', Min_1
      WRITE(*,'(A6,A,ES12.4E3)') &
        '', 'Min_2 = ', Min_2
    END IF

    nPP(1:nPS) = 0
    nPP(1)     = PRODUCT( nNodesX(1:3) )

    DO iDim = 1, 3

      IF( nNodesX(iDim) > 1 )THEN

        nPP(2*iDim:2*iDim+1) &
          = PRODUCT( nNodesX(1:3), MASK = [1,2,3] .NE. iDim )

      END IF

    END DO

    nPT = SUM( nPP(1:nPS) )

    ALLOCATE( U_PP(nPT,nCF) )
    ALLOCATE( G_PP(nPT,nGF) )

    ALLOCATE( L_X(nPT,nDOFX) )

    L_X = Zero
    DO iNX = 1, nDOFX

      L_X(iNX,iNX) = One

      IF( SUM( nPP(2:3) ) > 0 )THEN

        iOS = nPP(1)
        L_X(iOS+1:iOS+nDOFX_X1,iNX) = LX_X1_Dn(1:nDOFX_X1,iNX)

        iOS = iOS + nPP(2)
        L_X(iOS+1:iOS+nDOFX_X1,iNX) = LX_X1_Up(1:nDOFX_X1,iNX)

      END IF

      IF( SUM( nPP(4:5) ) > 0 )THEN

        iOS = SUM( nPP(1:3) )
        L_X(iOS+1:iOS+nDOFX_X2,iNX) = LX_X2_Dn(1:nDOFX_X2,iNX)

        iOS = iOS + nPP(4)
        L_X(iOS+1:iOS+nDOFX_X2,iNX) = LX_X2_Up(1:nDOFX_X2,iNX)

      END IF

      IF( SUM( nPP(6:7) ) > 0 )THEN

        iOS = SUM( nPP(1:5) )
        L_X(iOS+1:iOS+nDOFX_X3,iNX) = LX_X3_Dn(1:nDOFX_X3,iNX)

        iOS = iOS + nPP(6)
        L_X(iOS+1:iOS+nDOFX_X3,iNX) = LX_X3_Up(1:nDOFX_X3,iNX)

      END IF

    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to: L_X )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(  L_X )
#endif

  END SUBROUTINE InitializePositivityLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE FinalizePositivityLimiter_Euler_Relativistic_IDEAL

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( release: L_X )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC DELETE(       L_X )
#endif

    DEALLOCATE( L_X )
    DEALLOCATE( U_PP )
    DEALLOCATE( G_PP )

  END SUBROUTINE FinalizePositivityLimiter_Euler_Relativistic_IDEAL


  !> Iterate through the entire spatial domain and apply the positivity
  !> limiter from Qin et al., (2016), JCP, 315, 323 to each element.
  !> @param Theta_D minimum value to ensure physical density
  !> @param Theta_q minimum value to ensure physical internal energy-density
  !>        and velocity
  !> @todo Modify to accomodate GR
  SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL_gpu &
           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )
!real(dp)::ut(ndofx,ix_b1(1):ix_e1(1),ix_b1(2):ix_e1(2),ix_b1(3):ix_e1(3),ncf)
!integer::inx,ix1,ix2,ix3,icf
!
!!$ACC UPDATE HOST( U )
!ut=u
!    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL_cpu &
!           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )
!open(100,file='cpu')
!do icf=1,ncf
!do ix3=ix_b1(3),ix_e1(3)
!do ix2=ix_b1(2),ix_e1(2)
!do ix1=ix_b1(1),ix_e1(1)
!do inx=1,ndofx
!write(100,'(ES24.16E3)')u(inx,ix1,ix2,ix3,icf)
!enddo
!enddo
!enddo
!enddo
!enddo
!close(100)
!u=ut
!    CALL ApplyPositivityLimiter_Euler_Relativistic_IDEAL_gpu &
!           ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )
!!$ACC UPDATE HOST( U )
!open(100,file='gpu')
!do icf=1,ncf
!do ix3=ix_b1(3),ix_e1(3)
!do ix2=ix_b1(2),ix_e1(2)
!do ix1=ix_b1(1),ix_e1(1)
!do inx=1,ndofx
!write(100,'(ES24.16E3)')u(inx,ix1,ix2,ix3,icf)
!enddo
!enddo
!enddo
!enddo
!enddo
!close(100)

  END SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL


  SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL_gpu &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iNX, iX1, iX2, iX3, iCF, iPT, nX_K
    REAL(DP) :: Min_D, Min_K, Min_ESq, Theta_D, Theta_P, q

    INTEGER :: iErr(              iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    LOGICAL :: NegativeStates(2  ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: Theta_q(          iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: SqrtGm(1:nDOFX   ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    REAL(DP) :: U_Q(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_P(1:nPT  ,1:nCF,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: U_K(1:nCF        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    ! --- Scale factors ---

    REAL(DP) :: h1Q(1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h2Q(1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h3Q(1:nDOFX      ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    REAL(DP) :: h1P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h2P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: h3P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    ! --- Metric coefficients ---

    REAL(DP) :: g1P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: g2P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))
    REAL(DP) :: g3P(1:nPT        ,iX_B0(1):iX_E0(1), &
                                  iX_B0(2):iX_E0(2), &
                                  iX_B0(3):iX_E0(3))

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

    nX_K  = PRODUCT( iX_E0 - iX_B0 + 1 )
    nCF_K = nCF * nX_K

    CALL TimersStart_Euler( Timer_Euler_PL_CopyIn )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:    iX_B0, iX_E0, G, U ) &
    !$OMP MAP( alloc: iErr, NegativeStates, Theta_q, SqrtGm, &
    !$OMP             U_Q, U_P, U_K, &
    !$OMP             h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#elif defined(THORNADO_OACC)
    !$ACC ENTER DATA &
    !$ACC COPYIN(     iX_B0, iX_E0, G, U ) &
    !$ACC CREATE(     iErr, NegativeStates, Theta_q, SqrtGm, &
    !$ACC             U_Q, U_P, U_K, &
    !$ACC             h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyIn )

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( SqrtGm, h1Q, h2Q, h3Q, G )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iNX = 1, nDOFX

      SqrtGm(iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_SqrtGm)

      h1Q   (iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_h_1)
      h2Q   (iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_h_2)
      h3Q   (iNX,iX1,iX2,iX3) = G(iNX,iX1,iX2,iX3,iGF_h_3)

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
    !$ACC PRESENT( U_Q, U )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(5)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF
    DO iNX = 1, nDOFX

      U_Q(iNX,iCF,iX1,iX2,iX3) = U(iNX,iX1,iX2,iX3,iCF)

    END DO
    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputePointValues( nCF_K, iX_B0, iX_E0, U_Q, U_P )

    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h1Q, h1P )
    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h2Q, h2P )
    CALL ComputePointValues( nX_K , iX_B0, iX_E0, h3Q, h3P )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    CALL TimersStart_Euler( Timer_Euler_PL_Permute )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT(  iX_B0, iX_E0, g1P, g2P, g3P, h1P, h2P, h3P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iPT = 1, nPT

      g1P(iPT,iX1,iX2,iX3) = MAX( h1P(iPT,iX1,iX2,iX3)**2, SqrtTiny )
      g2P(iPT,iX1,iX2,iX3) = MAX( h2P(iPT,iX1,iX2,iX3)**2, SqrtTiny )
      g3P(iPT,iX1,iX2,iX3) = MAX( h3P(iPT,iX1,iX2,iX3)**2, SqrtTiny )

    END DO
    END DO
    END DO
    END DO

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(4)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(4) &
    !$ACC PRESENT( U_K, WeightsX_q, SqrtGm, U_Q )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(4)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)
    DO iCF = 1, nCF

      U_K(iCF,iX1,iX2,iX3) &
        = SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) * U_Q(:,iCF,iX1,iX2,iX3) ) &
          / SUM( WeightsX_q * SqrtGm(:,iX1,iX2,iX3) )

    END DO
    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_Permute )

    ! --- Limit Mass-Density ---

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_D, Min_K, Theta_D )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( U_K, U_Q, U_P, NegativeStates ) &
    !$ACC PRIVATE( Min_D, Min_K, Theta_D )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_D, Min_K, Theta_D )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      Min_D = Min_1 * U_K(iCF_D,iX1,iX2,iX3)

      Min_K = HUGE( 1.0_DP )
      DO iPT = 1, nPT

        Min_K &
          = MIN( Min_K, U_P(iPT,iCF_D,iX1,iX2,iX3) )

      END DO

      NegativeStates(1,iX1,iX2,iX3) = .FALSE.

      IF( Min_K .LT. Min_D )THEN

        NegativeStates(1,iX1,iX2,iX3) = .TRUE.

        Theta_D &
          =   ( U_K(iCF_D,iX1,iX2,iX3) - Min_D ) &
            / ( U_K(iCF_D,iX1,iX2,iX3) - Min_K )

        DO iNX = 1, nDOFX

          U_Q(iNX,iCF_D,iX1,iX2,iX3) &
            = U_K(iCF_D,iX1,iX2,iX3) &
                + Theta_D * ( U_Q(iNX,iCF_D,iX1,iX2,iX3) &
                                - U_K(iCF_D,iX1,iX2,iX3) )

        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    ! --- Recompute point values ---

    CALL TimersStart_Euler( Timer_Euler_PL_Integrate )

    CALL ComputePointValues( nCF_K, iX_B0, iX_E0, U_Q, U_P )

    CALL TimersStop_Euler( Timer_Euler_PL_Integrate )

    ! --- Limit q-function ---

    CALL TimersStart_Euler( Timer_Euler_PL_LimitCells )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_ESq, Theta_P )
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( U_K, U_Q, U_P, g1P, g2P, g3P, Theta_q, &
    !$ACC          iErr, NegativeStates ) &
    !$ACC PRIVATE( Min_ESq, Theta_P )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE( Min_ESq, Theta_P )
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      iErr(iX1,iX2,iX3) = 0

      IF( U_K(iCF_E,iX1,iX2,iX3) .LT. Zero ) iErr(iX1,iX2,iX3) = 01

      Min_ESq = Min_2 * U_K(iCF_E,iX1,iX2,iX3)**2
      iErr(iX1,iX2,iX3) = 0

      Theta_q(iX1,iX2,iX3) = One

      NegativeStates(2,iX1,iX2,iX3) = .FALSE.

      DO iPT = 1, nPT

        q = Computeq &
              ( U_P(iPT,iCF_D ,iX1,iX2,iX3), &
                U_P(iPT,iCF_S1,iX1,iX2,iX3), &
                U_P(iPT,iCF_S2,iX1,iX2,iX3), &
                U_P(iPT,iCF_S3,iX1,iX2,iX3), &
                U_P(iPT,iCF_E ,iX1,iX2,iX3), &
                g1P(iPT      ,iX1,iX2,iX3), &
                g2P(iPT      ,iX1,iX2,iX3), &
                g3P(iPT      ,iX1,iX2,iX3), &
                Min_ESq )

        IF( q .LT. Zero )THEN

          NegativeStates(2,iX1,iX2,iX3) = .TRUE.

          CALL SolveTheta_Bisection &
                 ( U_P(iPT,iCF_D ,iX1,iX2,iX3), &
                   U_P(iPT,iCF_S1,iX1,iX2,iX3), &
                   U_P(iPT,iCF_S2,iX1,iX2,iX3), &
                   U_P(iPT,iCF_S3,iX1,iX2,iX3), &
                   U_P(iPT,iCF_E ,iX1,iX2,iX3), &
                   U_K(    iCF_D ,iX1,iX2,iX3), &
                   U_K(    iCF_S1,iX1,iX2,iX3), &
                   U_K(    iCF_S2,iX1,iX2,iX3), &
                   U_K(    iCF_S3,iX1,iX2,iX3), &
                   U_K(    iCF_E ,iX1,iX2,iX3), &
                   g1P(iPT      ,iX1,iX2,iX3), &
                   g2P(iPT      ,iX1,iX2,iX3), &
                   g3P(iPT      ,iX1,iX2,iX3), &
                   Min_ESq, Theta_P, &
                   iErr(iX1,iX2,iX3) )

          Theta_q(iX1,iX2,iX3) = MIN( Theta_q(iX1,iX2,iX3), Theta_P )

        END IF

      END DO

    END DO
    END DO
    END DO

    ! --- Limit all variables towards cell-average ---

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3)
#elif defined(THORNADO_OACC)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRESENT( U, U_K, U_Q, Theta_q, NegativeStates )
#elif defined(THORNADO_OMP)
    !$OMP PARALLEL DO SIMD COLLAPSE(3)
#endif
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( NegativeStates(2,iX1,iX2,iX3) )THEN

        DO iCF = 1, nCF
        DO iNX = 1, nDOFX
  
          U(iNX,iX1,iX2,iX3,iCF) &
            = U_K(iCF,iX1,iX2,iX3) &
                + Theta_q(iX1,iX2,iX3) * ( U_Q(iNX,iCF,iX1,iX2,iX3) &
                                             - U_K(iCF,iX1,iX2,iX3) )
  
        END DO
        END DO

      ELSE

        DO iCF = 1, nCF
        DO iNX = 1, nDOFX
  
          U(iNX,iX1,iX2,iX3,iCF) = U_Q(iNX,iCF,iX1,iX2,iX3)
  
        END DO
        END DO

      END IF

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PL_LimitCells )

    CALL TimersStart_Euler( Timer_Euler_PL_CopyOut )

#if defined(THORNADO_OMP_OL)
    !$OMP TARGET EXIT DATA &
    !$OMP MAP( from:    U, iErr ) &
    !$OMP MAP( release: NegativeStates, Theta_q, iX_B0, iX_E0, SqrtGm, &
    !$OMP               G, U_Q, U_P, U_K, &
    !$OMP               h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#elif defined(THORNADO_OACC)
    !$ACC EXIT DATA &
    !$ACC COPYOUT(      U, iErr ) &
    !$ACC DELETE(       NegativeStates, Theta_q, iX_B0, iX_E0, SqrtGm, &
    !$ACC               G, U_Q, U_P, U_K, &
    !$ACC               h1Q, h2Q, h3Q, h1P, h2P, h3P, g1P, g2P, g3P )
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_CopyOut )

    CALL TimersStart_Euler( Timer_Euler_PL_ErrorCheck )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( iErr(iX1,iX2,iX3) .NE. 0 ) &
        CALL DescribeError_Euler( iErr(iX1,iX2,iX3) )

    END DO
    END DO
    END DO

#ifdef THORNADO_DEBUG_EULER_POSITIVITYLIMITER
    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      IF( U_K(iCF_E,iX1,iX2,iX3) .LT. Zero )THEN

        CALL DescribeError_Euler( 01 )

      END IF

    END DO
    END DO
    END DO
#endif

    CALL TimersStop_Euler( Timer_Euler_PL_ErrorCheck )

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL_gpu


  SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL_cpu &
    ( iX_B0, iX_E0, iX_B1, iX_E1, G, U )

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(DP), INTENT(in)    :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)

    INTEGER  :: iX1, iX2, iX3, iCF, iP
    REAL(DP) :: Theta_D, Theta_q, Theta_P
    REAL(DP) :: Min_K, Min_D, Min_ESq
    REAL(DP) :: U_q(nDOFX,nCF), G_q(nDOFX,nGF), U_K(nCF), q(nPT), SSq(nPT)
    REAL(DP) :: chi

    ! --- For de-bugging ---
    REAL(DP) :: q_K(nPT)

    IF( nDOFX .EQ. 1 ) RETURN

    IF( .NOT. UsePositivityLimiter ) RETURN

    CALL TimersStart_Euler( Timer_Euler_PositivityLimiter )

    DO iX3 = iX_B0(3), iX_E0(3)
    DO iX2 = iX_B0(2), iX_E0(2)
    DO iX1 = iX_B0(1), iX_E0(1)

      U_q(1:nDOFX,1:nCF) = U(1:nDOFX,iX1,iX2,iX3,1:nCF)
      G_q(1:nDOFX,1:nGF) = G(1:nDOFX,iX1,iX2,iX3,1:nGF)

      ! --- Check if cell-average of density is physical ---

      U_K(iCF_D) = SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                          * U_q(1:nDOFX,iCF_D) ) &
                     / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) )

      IF( U_K(iCF_D) .GT. Zero )THEN

        Min_D = Min_1 * U_K(iCF_D)

        CALL ComputePointValues_cpu( U_q(1:nDOFX,iCF_D), U_PP(1:nPT,iCF_D) )

        Min_K = MINVAL( U_PP(1:nPT,iCF_D) )

        IF( Min_K .LT. Min_D )THEN

          ! --- Limit point-values of density towards cell-average ---

          Theta_D &
            = ( U_K(iCF_D) - Min_D ) / ( U_K(iCF_D) - Min_K )

          U_q(1:nDOFX,iCF_D) &
            = U_K(iCF_D) + Theta_D * ( U_q(1:nDOFX,iCF_D) - U_K(iCF_D) )

        END IF

        ! --- Modify point-values of q if necessary ---

        DO iCF = 1, nCF

          CALL ComputePointValues_cpu( U_q(1:nDOFX,iCF), U_PP(1:nPT,iCF) )

        END DO

        CALL ComputePointValues_cpu( G_q(1:nDOFX,iGF_h_1), G_PP(1:nPT,iGF_h_1) )
        CALL ComputePointValues_cpu( G_q(1:nDOFX,iGF_h_2), G_PP(1:nPT,iGF_h_2) )
        CALL ComputePointValues_cpu( G_q(1:nDOFX,iGF_h_3), G_PP(1:nPT,iGF_h_3) )

        CALL ComputeGeometryX_FromScaleFactors( G_PP(1:nPT,1:nGF) )

        U_K(iCF_E) &
          = ( SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                     * U_q(1:nDOFX,iCF_E) ) &
                / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) ) )

        IF( U_K(iCF_E) .LT. Zero ) &
            CALL DescribeError_Euler( 01 )

        Min_ESq = Min_2 * U_K(iCF_E)**2

        CALL Computeq_cpu &
               ( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), Min_ESq, q(1:nPT) )

        IF( ANY( q(1:nPT) .LT. Zero ) )THEN

          DO iCF = 1, nCF

            U_K(iCF) = SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                              * U_q(1:nDOFX,iCF) ) &
                         / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) )

          END DO

          ! --- Compute cell-average of q ---
          DO iP = 1, nPT

            CALL Computeq_cpu( 1, U_K(1:nCF), G_PP(iP,1:nGF), Min_ESq, q_K(iP) )

          END DO

          ! --- Artificially inject some thermal energy if any q_K < 0 ---
          IF( ANY( q_K(1:nPT) .LT. Zero ) )THEN

            WRITE(*,*)
            WRITE(*,'(A)') 'Adding internal energy'
            WRITE(*,'(A,3I5.4)') 'iX1, iX2, iX3 = ', iX1, iX2, iX3

            ! --- Use maximum value of SSq(1:nPT) to determine how much
            !     internal energy to inject ---
            DO iP = 1, nPT

              SSq(iP) = U_K(iCF_S1)**2 / G_PP(iP,iGF_Gm_dd_11) &
                           + U_K(iCF_S2)**2 / G_PP(iP,iGF_Gm_dd_22) &
                           + U_K(iCF_S3)**2 / G_PP(iP,iGF_Gm_dd_33)

            END DO

            WRITE(*,'(A,ES24.16E3)') 'tau_K (old):', U_K(iCF_E)

            chi = ( Min_2 * U_K(iCF_E) - U_K(iCF_D) &
                      + SQRT( U_K(iCF_D)**2 + MAXVAL( SSq ) &
                                + Min_ESq ) ) / U_K(iCF_E) - One

            U_K(iCF_E) = U_K(iCF_E) * ( One + chi )

            WRITE(*,'(A,ES13.6E3)') &
              'Fractional amount of internal energy added: ', chi

            DO iCF = 1, nCF

              U_q(1:nDOFX,iCF) = U_K(iCF)

            END DO

          ELSE

            ! --- Solve for Theta_q such that all point-values
            !     of q are positive ---

            Theta_q = One

            DO iP = 1, nPT

              IF( q(iP) .LT. Zero )THEN

                CALL SolveTheta_Bisection_cpu &
                  ( U_PP(iP,1:nCF), U_K(1:nCF), G_PP(iP,1:nGF), Min_ESq, &
                    Theta_P, iX1, iX2, iX3, iP )

                Theta_q = MIN( Theta_q, Theta_P )

              END IF

            END DO

            ! --- Limit all variables towards cell-average ---

            DO iCF = 1, nCF

              U_q(1:nDOFX,iCF) &
                = U_K(iCF) + Theta_q * ( U_q(1:nDOFX,iCF) - U_K(iCF) )

            END DO

          END IF

        END IF ! q < 0

        U(1:nDOFX,iX1,iX2,iX3,1:nCF) = U_q(1:nDOFX,1:nCF)

      ELSE

        WRITE(*,'(A)') &
          'Warning: ApplyPositivityLimiter_Euler_Relativistic_IDEAL'
        WRITE(*,'(A)') 'Cell-average of density <= 0'
        WRITE(*,'(A)') 'Setting all values to cell-average'

        DO iCF = 1, nCF

          U_K(iCF) &
            = SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) &
                     * U_q(1:nDOFX,iCF) ) &
              / SUM( WeightsX_q(1:nDOFX) * G_q(1:nDOFX,iGF_SqrtGm) )

          U(1:nDOFX,iX1,iX2,iX3,iCF) = U_K(iCF)

        END DO

      END IF

      ! --- Ensure that all point-values of q are positive after limiting ---

      DO iCF = 1, nCF

        CALL ComputePointValues_cpu( U(1:nDOFX,iX1,iX2,iX3,iCF), U_PP(1:nPT,iCF) )

      END DO

      CALL ComputePointValues_cpu &
             ( G(1:nDOFX,iX1,iX2,iX3,iGF_h_1), G_PP(1:nPT,iGF_h_1) )
      CALL ComputePointValues_cpu &
             ( G(1:nDOFX,iX1,iX2,iX3,iGF_h_2), G_PP(1:nPT,iGF_h_2) )
      CALL ComputePointValues_cpu &
             ( G(1:nDOFX,iX1,iX2,iX3,iGF_h_3), G_PP(1:nPT,iGF_h_3) )

      CALL ComputeGeometryX_FromScaleFactors( G_PP(1:nPT,1:nGF) )

      CALL Computeq_cpu &
             ( nPT, U_PP(1:nPT,1:nCF), G_PP(1:nPT,1:nGF), Min_ESq, q(1:nPT) )

      IF( ANY( q(1:nPT) .LT. Zero ) ) &
          CALL DescribeError_Euler( 04 )

    END DO
    END DO
    END DO

    CALL TimersStop_Euler( Timer_Euler_PositivityLimiter )

  END SUBROUTINE ApplyPositivityLimiter_Euler_Relativistic_IDEAL_cpu


  SUBROUTINE ComputePointValues_cpu( X_Q, X_P )

    REAL(DP), INTENT(in)  :: X_Q(nDOFX)
    REAL(DP), INTENT(out) :: X_P(nPT)

    INTEGER :: iOS

    X_P(1:nDOFX) = X_Q(1:nDOFX)

    IF( SUM( nPP(2:3) ) > 0 )THEN

      ! --- X1 ---

      iOS = nPP(1)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Dn, nDOFX_X1, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X1), 1 )

      iOS = iOS + nPP(2)

      CALL DGEMV &
             ( 'N', nDOFX_X1, nDOFX, One, LX_X1_Up, nDOFX_X1, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X1), 1 )

    END IF

    IF( SUM( nPP(4:5) ) > 0 )THEN

      ! --- X2 ---

      iOS = SUM( nPP(1:3) )

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Dn, nDOFX_X2, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X2), 1 )

      iOS = iOS + nPP(4)

      CALL DGEMV &
             ( 'N', nDOFX_X2, nDOFX, One, LX_X2_Up, nDOFX_X2, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X2), 1 )

    END IF

    IF( SUM( nPP(6:7) ) > 0 )THEN

      ! --- X3 ---

      iOS = SUM( nPP(1:5) )

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Dn, nDOFX_X3, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X3), 1 )

      iOS = iOS + nPP(6)

      CALL DGEMV &
             ( 'N', nDOFX_X3, nDOFX, One, LX_X3_Up, nDOFX_X3, &
               X_Q(1:nDOFX), 1, Zero, X_P(iOS+1:iOS+nDOFX_X3), 1 )

    END IF

  END SUBROUTINE ComputePointValues_cpu


  SUBROUTINE Computeq_cpu( N, U, G, Min_ESq, q )

    INTEGER,  INTENT(in)  :: N
    REAL(DP), INTENT(in)  :: U(N,nCF), G(N,nGF), Min_ESq
    REAL(DP), INTENT(out) :: q(N)

    q = qFun( U(1:N,iCF_D ), &
              U(1:N,iCF_S1), &
              U(1:N,iCF_S2), &
              U(1:N,iCF_S3), &
              U(1:N,iCF_E ), &
              G(1:N,iGF_Gm_dd_11), &
              G(1:N,iGF_Gm_dd_22), &
              G(1:N,iGF_Gm_dd_33), &
              Min_ESq )

  END SUBROUTINE Computeq_cpu


  REAL(DP) ELEMENTAL FUNCTION qFun &
    ( D, S1, S2, S3, tau, Gm11, Gm22, Gm33, Min_ESq )

    REAL(DP), INTENT(in) :: D, S1, S2, S3, tau, Gm11, Gm22, Gm33, Min_ESq

    qFun = tau + D &
             - SQRT( D**2 + ( S1**2 / Gm11 + S2**2 / Gm22 + S3**2 / Gm33 ) &
                       + Min_ESq )

    RETURN
  END FUNCTION qFun


  SUBROUTINE SolveTheta_Bisection_cpu &
    ( U_Q, U_K, G_Q, Min_ESq, Theta_P, iX1, iX2, iX3, iP )

    REAL(DP), INTENT(in)  :: U_Q(nCF), U_K(nCF), G_Q(nGF), Min_ESq
    REAL(DP), INTENT(out) :: Theta_P

    ! --- For de-bugging ---
    INTEGER, INTENT(in)  :: iX1, iX2, iX3, iP
!!$    INTEGER              :: iCF, iGF

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0d-3

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

!!$    WRITE(*,*)
!!$    WRITE(*,'(A,I3)') 'iP = ', iP
!!$    WRITE(*,'(A)') 'U_Q = np.array( ['
!!$    DO iCF = 1, nCF-1
!!$      WRITE(*,'(ES24.16E3,A)') U_Q(iCF), ','
!!$    END DO
!!$    WRITE(*,'(ES24.16E3,A)') U_Q(iCF), ' ] )'
!!$    WRITE(*,*)
!!$    WRITE(*,'(A)') 'U_K = np.array( ['
!!$    DO iCF = 1, nCF-1
!!$      WRITE(*,'(ES24.16E3,A)') U_K(iCF), ','
!!$    END DO
!!$    WRITE(*,'(ES24.16E3,A)') U_K(iCF), ' ] )'
!!$    WRITE(*,*)
!!$    WRITE(*,'(A)') 'G_Q = np.array( ['
!!$    DO iGF = iGF_Gm_dd_11, iGF_Gm_dd_33-1
!!$      WRITE(*,'(ES24.16E3,A)') G_Q(iCF), ','
!!$    END DO
!!$    WRITE(*,'(ES24.16E3,A)') G_Q(iGF), ' ] )'
!!$    WRITE(*,*)

    x_a = Zero
    f_a = qFun &
            ( x_a * U_Q(iCF_D)  + ( One - x_a ) * U_K(iCF_D),         &
              x_a * U_Q(iCF_S1) + ( One - x_a ) * U_K(iCF_S1),        &
              x_a * U_Q(iCF_S2) + ( One - x_a ) * U_K(iCF_S2),        &
              x_a * U_Q(iCF_S3) + ( One - x_a ) * U_K(iCF_S3),        &
              x_a * U_Q(iCF_E)  + ( One - x_a ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33), Min_ESq )

    x_b = One
    f_b = qFun &
            ( x_b * U_Q(iCF_D)  + ( One - x_b ) * U_K(iCF_D),         &
              x_b * U_Q(iCF_S1) + ( One - x_b ) * U_K(iCF_S1),        &
              x_b * U_Q(iCF_S2) + ( One - x_b ) * U_K(iCF_S2),        &
              x_b * U_Q(iCF_S3) + ( One - x_b ) * U_K(iCF_S3),        &
              x_b * U_Q(iCF_E)  + ( One - x_b ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33), Min_ESq )

    IF( .NOT. f_a * f_b < 0 ) &
      CALL DescribeError_Euler( 02 )

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx
      f_c = qFun &
            ( x_c * U_Q(iCF_D)  + ( One - x_c ) * U_K(iCF_D),         &
              x_c * U_Q(iCF_S1) + ( One - x_c ) * U_K(iCF_S1),        &
              x_c * U_Q(iCF_S2) + ( One - x_c ) * U_K(iCF_S2),        &
              x_c * U_Q(iCF_S3) + ( One - x_c ) * U_K(iCF_S3),        &
              x_c * U_Q(iCF_E)  + ( One - x_c ) * U_K(iCF_E),         &
              G_Q(iGF_Gm_dd_11), G_Q(iGF_Gm_dd_22), G_Q(iGF_Gm_dd_33), Min_ESq )

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx .LT. dx_min ) CONVERGED = .TRUE.

      IF( ITERATION .GT. MAX_IT .AND. .NOT. CONVERGED )THEN

        WRITE(*,'(6x,A)') &
          'SolveTheta_Bisection (ApplyPositivityLimiter_Euler_Relativistic_IDEAL):'
        WRITE(*,'(A8,A,I3.3)') &
          '', 'ITERATION ', ITERATION
        WRITE(*,'(A8,A,4ES15.6e3)') &
          '', 'x_a, x_c, x_b, dx = ', x_a, x_c, x_b, dx
        WRITE(*,'(A8,A,3ES15.6e3)') &
          '', 'f_a, f_c, f_b     = ', f_a, f_c, f_b

        IF( ITERATION .GT. MAX_IT + 3 ) &
          CALL DescribeError_Euler( 03 )

      END IF

    END DO

    Theta_P = x_a

  END SUBROUTINE SolveTheta_Bisection_cpu

  ! --- GPU ---

  SUBROUTINE ComputePointValues( N, iX_B0, iX_E0, U_Q, U_P )

    INTEGER,  INTENT(in)  :: N, iX_B0(3), iX_E0(3)
    REAL(DP), INTENT(in)  :: U_Q(1:nDOFX,1:nCF,iX_B0(1):iX_E0(1), &
                                               iX_B0(2):iX_E0(2), &
                                               iX_B0(3):iX_E0(3) )
    REAL(DP), INTENT(out) :: U_P(1:nPT  ,1:nCF,iX_B0(1):iX_E0(1), &
                                               iX_B0(2):iX_E0(2), &
                                               iX_B0(3):iX_E0(3) )

    CALL MatrixMatrixMultiply &
           ( 'N', 'N', nPT, N, nDOFX, One, L_X, nPT, &
             U_Q, nDOFX, Zero, U_P, nPT )

  END SUBROUTINE ComputePointValues


  FUNCTION Computeq &
    ( D, S1, S2, S3, tau, g1, g2, g3, Min_ESq ) RESULT( q )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in) :: D, S1, S2, S3, tau, g1, g2, g3, Min_ESq

    REAL(DP) :: q

    q = tau + D &
          - SQRT( D**2 + ( S1**2 / g1 + S2**2 / g2 + S3**2 / g3 ) &
                  + Min_ESq )

    RETURN
  END FUNCTION Computeq


  SUBROUTINE SolveTheta_Bisection &
    ( D_P, S1_P, S2_P, S3_P, E_P, &
      D_K, S1_K, S2_K, S3_K, E_K, &
      g1_P, g2_P, g3_P, Min_ESq, Theta_P, iErr )

#if defined(THORNADO_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(THORNADO_OACC)
    !$ACC ROUTINE SEQ
#endif

    REAL(DP), INTENT(in)    :: D_P, S1_P, S2_P, S3_P, E_P, &
                               D_K, S1_K, S2_K, S3_K, E_K, &
                               g1_P, g2_P, g3_P, Min_ESq
    REAL(DP), INTENT(out)   :: Theta_P
    INTEGER , INTENT(inout) :: iErr

    INTEGER,  PARAMETER :: MAX_IT = 19
    REAL(DP), PARAMETER :: dx_min = 1.0e-3_DP

    LOGICAL  :: CONVERGED
    INTEGER  :: ITERATION
    REAL(DP) :: x_a, x_b, x_c, dx
    REAL(DP) :: f_a, f_b, f_c

    x_a = Zero
    f_a = Computeq &
            ( x_a *  D_P + ( One - x_a ) *  D_K, &
              x_a * S1_P + ( One - x_a ) * S1_K, &
              x_a * S2_P + ( One - x_a ) * S2_K, &
              x_a * S3_P + ( One - x_a ) * S3_K, &
              x_a *  E_P + ( One - x_a ) *  E_K, &
              g1_P, g2_P, g3_P, Min_ESq )

    x_b = One
    f_b = Computeq &
            ( x_b *  D_P + ( One - x_b ) *  D_K, &
              x_b * S1_P + ( One - x_b ) * S1_K, &
              x_b * S2_P + ( One - x_b ) * S2_K, &
              x_b * S3_P + ( One - x_b ) * S3_K, &
              x_b *  E_P + ( One - x_b ) *  E_K, &
              g1_P, g2_P, g3_P, Min_ESq )

    IF( .NOT. f_a * f_b < 0 ) iErr = 02

    dx = x_b - x_a

    ITERATION = 0
    CONVERGED = .FALSE.
    DO WHILE ( .NOT. CONVERGED )

      ITERATION = ITERATION + 1

      dx = Half * dx
      x_c = x_a + dx
      f_c = Computeq &
              ( x_c *  D_P + ( One - x_c ) *  D_K, &
                x_c * S1_P + ( One - x_c ) * S1_K, &
                x_c * S2_P + ( One - x_c ) * S2_K, &
                x_c * S3_P + ( One - x_c ) * S3_K, &
                x_c *  E_P + ( One - x_c ) *  E_K, &
                g1_P, g2_P, g3_P, Min_ESq )

      IF( f_a * f_c < Zero )THEN

        x_b = x_c
        f_b = f_c

      ELSE

        x_a = x_c
        f_a = f_c

      END IF

      IF( dx .LT. dx_min ) CONVERGED = .TRUE.

      IF( ITERATION .GT. MAX_IT .AND. .NOT. CONVERGED )THEN

        IF( ITERATION .GT. MAX_IT + 3 ) iErr = 03

      END IF

    END DO

    Theta_P = x_a

  END SUBROUTINE SolveTheta_Bisection


END MODULE Euler_PositivityLimiterModule_Relativistic_IDEAL
