program process_option_driver
  use crest_data,only:status_normal
  use mtd_process_scheduler,only:resolve_mtd_process_isolation
  implicit none
  logical :: requested
  integer :: status
  character(len=256) :: message

  call resolve_mtd_process_isolation(requested,status,message)
  write(*,'(a,l1)') 'requested=',requested
  write(*,'(a,i0)') 'status=',status
  if (status /= status_normal) write(*,'(a)') 'message='//trim(message)
end program process_option_driver
