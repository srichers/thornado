PROGRAM main

  ! --- AMReX Modules ---
  USE amrex_fort_module,                ONLY: &
    amrex_real
  USE amrex_parallel_module,            ONLY: &
    amrex_parallel_ioprocessor, &
    amrex_parallel_communicator

  ! --- thornado Modules ---
  USE MeshModule,                       ONLY: &
    MeshX, &
    DestroyMesh
  USE InputOutputModuleAMReX,           ONLY: &
    ReadCheckpointFile,          &
    WriteFieldsAMReX_Checkpoint, &
    WriteFieldsAMReX_PlotFile
  USE UnitsModule,                      ONLY: &
    Millisecond
  USE TimersModule_Euler,               ONLY: &
    TimeIt_Euler,           &
    InitializeTimers_Euler, &
    FinalizeTimers_Euler

  ! --- Local Modules ---
  USE MF_Euler_UtilitiesModule,         ONLY: &
    MF_ComputeFromConserved, &
    MF_ComputeTimeStep
  USE MF_Euler_SlopeLimiterModule,      ONLY: &
    MF_ApplySlopeLimiter_Euler
  USE MF_Euler_PositivityLimiterModule, ONLY: &
    MF_ApplyPositivityLimiter_Euler
  USE MF_Euler_dgDiscretizationModule,  ONLY: &
    MF_Euler_ComputeIncrement
  USE MF_TimeSteppingModule_SSPRK,      ONLY: &
    MF_UpdateFluid_SSPRK
  USE FinalizationModule,               ONLY: &
    FinalizeProgram
  USE MF_UtilitiesModule,               ONLY: &
    ShowVariableFromMultifab
  USE MyAmrDataModule,                  ONLY: &
    MF_uGF, &
    MF_uCF, &
    MF_uPF, &
    MF_uAF
  USE InitializationModule,             ONLY: &
    InitializeProgram, &
    chk,               &
    wrt
  USE MyAmrModule,                      ONLY: &
    nLevels,   &
    StepNo,    &
    t,         &
    dt,        &
    t_end,     &
    CFL,       &
    t_wrt,     &
    dt_wrt,    &
    t_chk,     &
    dt_chk,    &
    iCycleD,   &
    iCycleW,   &
    iCycleChk, &
    GEOM
  USE UnitsModule,                      ONLY: &
    UnitsDisplay
  USE TimersModule_AMReX_Euler,         ONLY: &
    TimeIt_AMReX_Euler,            &
    InitializeTimers_AMReX_Euler,  &
    FinalizeTimers_AMReX_Euler,    &
    TimersStart_AMReX_Euler,       &
    TimersStop_AMReX_Euler,        &
    Timer_AMReX_Euler_InputOutput, &
    Timer_AMReX_Euler_MPI_Barrier

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER          :: iErr
  REAL(amrex_real) :: Timer_Evolution

  CALL InitializeProgram

  TimeIt_AMReX_Euler = .TRUE.
  CALL InitializeTimers_AMReX_Euler

  TimeIt_Euler = .TRUE.
  CALL InitializeTimers_Euler

  IF( amrex_parallel_ioprocessor() ) &
      Timer_Evolution = MPI_WTIME()

  CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_MPI_Barrier )
  CALL MPI_BARRIER( amrex_parallel_communicator(), iErr )
  CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_MPI_Barrier )

  DO WHILE( ALL( t .LT. t_end ) )

    StepNo = StepNo + 1

    CALL MF_ComputeTimeStep( MF_uGF, MF_uCF, CFL, dt )

    IF( ALL( t + dt .LE. t_end ) )THEN
      t = t + dt
    ELSE
      dt = t_end - [t]
      t  = [t_end]
    END IF

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )
    IF( amrex_parallel_ioprocessor() )THEN
      IF( MOD( StepNo(0), iCycleD ) .EQ. 0 )THEN
        WRITE(*,'(8x,A8,I8.8,A5,ES13.6E3,1x,A,A6,ES13.6E3,1x,A)') &
          'StepNo: ', StepNo(0), ' t = ', t / UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel ), &
          ' dt = ', dt(0) /  UnitsDisplay % TimeUnit, &
          TRIM( UnitsDisplay % TimeLabel )
     END IF
    END IF
    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

    CALL MF_UpdateFluid_SSPRK &
           ( t, dt, MF_uGF, MF_uCF, &
             GEOM, MF_Euler_ComputeIncrement )

    CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )
    IF( iCycleChk .GT. 0 )THEN
      IF( MOD( StepNo(0), iCycleChk ) .EQ. 0 ) &
        chk = .TRUE.
    ELSE
      IF( ALL( t + dt .GT. t_chk ) )THEN
        t_chk = t_chk + dt_chk
        chk   = .TRUE.
      END IF
    END IF

    IF( chk )THEN

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_Checkpoint &
             ( StepNo, nLevels, dt, t, t_wrt, &
               MF_uGF % BA % P, &
               MF_uGF % P, &
               MF_uCF % P, &
               MF_uPF % P, &
               MF_uAF % P )

      chk = .FALSE.

    END IF

    IF( iCycleW .GT. 0 )THEN
      IF( MOD( StepNo(0), iCycleW ) .EQ. 0 ) &
        wrt = .TRUE.
    ELSE
      IF( ALL( t + dt .GT. t_wrt ) )THEN
        t_wrt = t_wrt + dt_wrt
        wrt   = .TRUE.
      END IF
    END IF

    IF( wrt )THEN

      CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

      CALL WriteFieldsAMReX_PlotFile &
             ( t(0), StepNo, &
               MF_uGF_Option = MF_uGF, &
               MF_uCF_Option = MF_uCF, &
               MF_uPF_Option = MF_uPF, &
               MF_uAF_Option = MF_uAF )

      wrt = .FALSE.

    END IF
    CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

  END DO

  ! --- END of evolution ---

  IF( amrex_parallel_ioprocessor() )THEN
    WRITE(*,*)
    WRITE(*,'(2x,A,ES13.6E3,A)') &
      'Total evolution time: ', MPI_WTIME() - Timer_Evolution, ' s'
  END IF

  StepNo = StepNo + 1

  CALL TimersStart_AMReX_Euler( Timer_AMReX_Euler_InputOutput )
  CALL MF_ComputeFromConserved( MF_uGF, MF_uCF, MF_uPF, MF_uAF )

  CALL WriteFieldsAMReX_Checkpoint &
         ( StepNo, nLevels, dt, t, t_wrt, &
           MF_uGF % BA % P, &
           MF_uGF % P, &
           MF_uCF % P, &
           MF_uPF % P, &
           MF_uAF % P )

  CALL WriteFieldsAMReX_PlotFile &
         ( t(0), StepNo,           &
           MF_uGF_Option = MF_uGF, &
           MF_uCF_Option = MF_uCF, &
           MF_uPF_Option = MF_uPF, &
           MF_uAF_Option = MF_uAF )
  CALL TimersStop_AMReX_Euler( Timer_AMReX_Euler_InputOutput )

  ! --- Finalize everything ---

  CALL FinalizeTimers_Euler &
         ( Verbose_Option = amrex_parallel_ioprocessor(), &
           SuppressApplicationDriver_Option = .TRUE. )

  CALL FinalizeTimers_AMReX_Euler

  CALL FinalizeProgram( GEOM, MeshX )

END PROGRAM main
