PROGRAM Driver

  USE KindModule, ONLY: &
    DP, Pi, TwoPi, Half, Zero
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kilometer, &
    Millisecond, &
    MeV, &
    Erg
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE GeometryFieldsModule, ONLY: &
    uGF, iGF_SqrtGm, iGF_h_1, iGF_h_2, iGF_h_3
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  USE RadiationFieldsModule, ONLY: &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE ThornadoInitializationModule, ONLY: &
    InitThornado, &       ! --- To be called once
    InitThornado_Patch, & ! --- To be called once per patch
    FreeThornado_Patch    ! --- To be called once per parch
  USE TimeSteppingModule_Flash, ONLY: &
    Update_IMEX_PDARS
  USE SubcellReconstructionModule, ONLY: &
    UpdateSubcellReconstruction, &
    Phi_ijq, ProjectionMatrix
  USE MeshModule, ONLY: &
    MeshX
  USE ReferenceElementModuleX, ONLY: &
    NodesX1, WeightsX1, &
    NodesX2, WeightsX2, &
    NodesX3, WeightsX3
  USE ProgramHeaderModule, ONLY: &
    nNodesX

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER  :: i, iDOFX, iX1, iQ, iS, jS
  INTEGER  :: mpierr
  INTEGER, PARAMETER  :: nN = 2
  INTEGER  :: nX(3)
  REAL(DP) :: wTime
  REAL(DP) :: dt
  REAL(DP) :: Volume_Si, SQRTGamma_iq, Inte_Phi_ij
  REAL(DP) :: New_ProjectionMatrix(nN,nN)
  REAL(DP) :: Subcell_Lo, Subcell_Width

  nX = [12, 1, 1]

  CALL MPI_INIT( mpierr )

  wTime = MPI_WTIME( )

  CALL InitThornado &
         ( nNodes = nN, nDimsX = 1, nE = 10, swE = 0, &
           eL_MeV = 0.0d0, eR_MeV = 1.0d2, zoomE = 1.0_DP )

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A,A32,ES8.2E2)') '', 'InitThornado: ', wTime
  WRITE(*,*)

  wTime = MPI_WTIME( )

  WRITE(*,*) 'In Cartesian coordinate:'
  WRITE(*,'(A25,4ES12.3)') 'ProjectionMatrix', ProjectionMatrix

  WRITE(*,*)
  WRITE(*,*) 'In Spherical polar coordinate: '
  DO i = 1, 1

    CALL InitThornado_Patch &
           ( nX    = nX, &
             swX   = [ 0, 0, 0 ], &
             xL    = [ 00.0_DP            , 00.0_DP, 00.0_DP ], &
             xR    = [ 12.0_DP * Kilometer, Pi     , TwoPi ], &
             nSpecies = 1, CoordinateSystem_Option = 'spherical' )

    DO iX1 = 1, nX(1)

      WRITE(*,'(A25,I5)') 'iX1', iX1

      CALL UpdateSubcellReconstruction( [iX1, 1, 1] )

      WRITE(*,'(A25,4ES12.3)') 'ProjectionMatrix', ProjectionMatrix

      New_ProjectionMatrix = Zero

      Subcell_Width = MeshX(1) % Width(iX1) / nNodesX(1)

      DO iS = 1, nN

        ! the low bry of the subcell
        Subcell_Lo = MeshX(1) % Center(iX1) + (iS - 2) * MeshX(1) % Width(iX1) * Half

        ! --- integral for volume Si ---
        Volume_Si = 0.0_DP
        DO iQ = 1, nNodesX(1)
          ! quadrature point location
          SQRTGamma_iq = Subcell_Lo &
            + Subcell_Width * Half + NodesX1(iQ) * Subcell_Width
          ! sqrt gamma, which = 1 for cartesian
          SQRTGamma_iq = SQRTGamma_iq**2
          ! integral with sqrt gamma, Volume_Si = 1 for cartesian
          Volume_Si = Volume_Si + WeightsX1(iQ) * SQRTGamma_iq
        END DO

        DO jS = 1, nN

          ! --- integral for Phi_j at ith subcell ---
          Inte_Phi_ij = 0.0_DP
          DO iQ = 1, nN
            ! quadrature point location
            SQRTGamma_iq = Subcell_Lo &
            + Subcell_Width * Half + NodesX1(iQ) * Subcell_Width
            ! sqrt gamma
            SQRTGamma_iq = SQRTGamma_iq**2
            ! integral wt sqrt gamma
            Inte_Phi_ij = Inte_Phi_ij &
              + WeightsX1(iQ) * SQRTGamma_iq * Phi_ijq(iS,jS,iQ) 
          END DO
  
          New_ProjectionMatrix(iS,jS) = Inte_Phi_ij / Volume_Si

        END DO
      END DO

      WRITE(*,'(A25,4ES12.3)') 'New_ProjectionMatrix', New_ProjectionMatrix

    END DO

    CALL FreeThornado_Patch

  END DO

  wTime = MPI_WTIME( ) - wTime

  WRITE(*,*)
  WRITE(*,'(A,A32,ES8.2E2)') '', 'InitThornado_Patch: ', wTime
  WRITE(*,*)

  CALL MPI_FINALIZE( mpierr )

END PROGRAM Driver
