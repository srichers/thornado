PROGRAM WriteProjectionMatrix

  USE KindModule, ONLY: &
    DP
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE Euler_MeshRefinementModule, ONLY: &
    InitializeMeshRefinement_Euler, &
    FinalizeMeshRefinement_Euler

  IMPLICIT NONE

  CHARACTER(64) :: ProgramName, CoordinateSystem

  INTEGER :: nX(3), swX(3), bcX(3), nNodes
  REAL(DP) :: xL(3), xR(3), ZoomX(3)
  LOGICAL :: ActivateUnits

  ProgramName = 'Advection2D'
  CoordinateSystem = 'CARTESIAN'

  nX = [ 32, 32, 1 ]
  swX = [ 1, 1, 0 ]
  bcX = [ 1, 1, 0 ]

  xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
  xR = [ 1.0_DP, 1.0_DP, 1.0_DP ]

  ZoomX = [ 1.0_DP, 1.0_DP, 1.0_DP ]

  nNodes = 3

  ActivateUnits = .FALSE.

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = swX, &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           ActivateUnits_Option &
             = ActivateUnits, &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL InitializeMeshRefinement_Euler

  CALL FinalizeMeshRefinement_Euler

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementX

  CALL FinalizeProgram

END PROGRAM WriteProjectionMatrix
