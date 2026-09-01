!================================================================================!
! Explicit interface for multilevel optimization's optional in-memory ingress.
!================================================================================!
module crest_multilevel_interface
  implicit none
  private
  public :: crest_multilevel_oloop

  interface
    subroutine crest_multilevel_oloop(env,ensnam,multilevel_in,input_buffer)
      use crest_data,only:systemdata
      use crest_poststage_ensemble,only:poststage_ensemble
      implicit none
      type(systemdata),intent(inout) :: env
      character(len=*),intent(in) :: ensnam
      logical,intent(in) :: multilevel_in(6)
      type(poststage_ensemble),intent(inout),optional :: input_buffer
    end subroutine crest_multilevel_oloop
  end interface
end module crest_multilevel_interface
