program optimizer_affinity_fortran_driver
  use iso_c_binding,only:c_int
  use omp_lib
  implicit none
  integer,parameter :: test_threads = 4
  integer,parameter :: test_tasks = 32

  interface
    integer(c_int) function c_test_prepare(expected_threads,fail_bind_rank, &
    & fail_restore_rank) bind(C,name='crest_optimizer_affinity_test_prepare')
      import :: c_int
      integer(c_int),value :: expected_threads,fail_bind_rank,fail_restore_rank
    end function c_test_prepare

    integer(c_int) function c_bind(rank,team_size) &
    & bind(C,name='crest_optimizer_affinity_bind')
      import :: c_int
      integer(c_int),value :: rank,team_size
    end function c_bind

    integer(c_int) function c_restore(rank,team_size) &
    & bind(C,name='crest_optimizer_affinity_restore')
      import :: c_int
      integer(c_int),value :: rank,team_size
    end function c_restore

    integer(c_int) function c_finish() &
    & bind(C,name='crest_optimizer_affinity_finish')
      import :: c_int
    end function c_finish
  end interface

  call omp_set_dynamic(.false.)
  call omp_set_num_threads(test_threads)
  call run_case('normal_first',-1,-1,0,0,0,test_tasks)
  call run_case('bind_injected',2,-1,0,1,0,0)
  call run_case('team_mismatch',-1,-1,-1,test_threads,0,0)
  call run_case('restore_injected',-1,1,0,0,1,test_tasks)
  call run_case('normal_second',-1,-1,0,0,0,test_tasks)
  write(*,'(a)') 'OPTIMIZER_AFFINITY_FORTRAN_ORACLE_PASS'

contains

  subroutine run_case(label,fail_bind_rank,fail_restore_rank,team_delta, &
  & expected_bind_errors,expected_restore_errors,expected_tasks)
    character(len=*),intent(in) :: label
    integer,intent(in) :: fail_bind_rank,fail_restore_rank,team_delta
    integer,intent(in) :: expected_bind_errors,expected_restore_errors
    integer,intent(in) :: expected_tasks
    integer :: rank,team_size,reported_team,i
    integer :: bind_errors,restore_errors,task_count,observed_team
    integer(c_int) :: cstat,finish_status

    cstat = c_test_prepare(int(test_threads,c_int),int(fail_bind_rank,c_int), &
    & int(fail_restore_rank,c_int))
    if (cstat /= 0_c_int) error stop 'Fortran oracle prepare failed'
    bind_errors = 0
    restore_errors = 0
    task_count = 0
    observed_team = 0
    !$omp parallel default(none) &
    !$omp shared(bind_errors,restore_errors,task_count,observed_team,team_delta) &
    !$omp private(rank,team_size,reported_team,cstat,i)
    rank = omp_get_thread_num()
    team_size = omp_get_num_threads()
    reported_team = team_size+team_delta
    !$omp single
    observed_team = team_size
    !$omp end single
    cstat = c_bind(int(rank,c_int),int(reported_team,c_int))
    if (cstat /= 0_c_int) then
      !$omp atomic update
      bind_errors = bind_errors+1
    end if
    !$omp barrier
    !$omp single
    if (bind_errors == 0) then
      do i = 1,test_tasks
        !$omp task shared(task_count)
        !$omp atomic update
        task_count = task_count+1
        !$omp end task
      end do
      !$omp taskwait
    end if
    !$omp end single
    cstat = c_restore(int(rank,c_int),int(team_size,c_int))
    if (cstat /= 0_c_int) then
      !$omp atomic update
      restore_errors = restore_errors+1
    end if
    !$omp end parallel
    finish_status = c_finish()

    if (observed_team /= test_threads) error stop 'Fortran oracle team size mismatch'
    if (bind_errors /= expected_bind_errors) error stop 'Fortran oracle bind count mismatch'
    if (restore_errors /= expected_restore_errors) &
    & error stop 'Fortran oracle restore count mismatch'
    if (task_count /= expected_tasks) error stop 'Fortran oracle task count mismatch'
    if (finish_status /= 0_c_int) error stop 'Fortran oracle finish failed'
    write(*,'(a,a,a,i0,a,i0,a,i0)') 'FORTRAN_CASE_PASS label=',trim(label), &
    & ' bind_errors=',bind_errors,' restore_errors=',restore_errors, &
    & ' tasks=',task_count
  end subroutine run_case
end program optimizer_affinity_fortran_driver
