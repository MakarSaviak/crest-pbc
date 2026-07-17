!> Focused driver for COM-CV tests
program comcv_tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  use test_comcv, only : collect_comcv
  implicit none
  integer :: stat
  type(testsuite_type) :: suite

  stat = 0
  suite = new_testsuite("comcv", collect_comcv)
  write(error_unit, '("# Testing:",1x,a)') suite%name
  call run_testsuite(suite%collect, error_unit, stat)
  if (stat > 0) then
    write(error_unit, '(i0,1x,a)') stat, "test(s) failed!"
    error stop 1
  end if
end program comcv_tester
