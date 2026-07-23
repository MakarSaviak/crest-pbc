!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!
module modelhessian_module
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use crest_calculator,only:calcdata,constrhess
  use modelhessian_core
  implicit none

  public :: modhes,modhes_free

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!
!
  subroutine modhes(calc,modh,natoms,xyz,at,Hess,pr)
!**********************************************************
!* subroutine modhes
!* create a model Hessian for a given molecule
!*
!* Input:
!*     natoms - number of atoms
!*       xyz  - Cartesian coordinates
!*        at  - atom types as integers
!*      modh  - model Hessian settings (see above)
!*      calc  - calculation settings (for constraints)
!*        pr  - printout selection
!*
!* Output:
!*      Hess  - the (packed) model Hessian
!**********************************************************
    implicit none
    type(calcdata),intent(in) :: calc
    type(mhparam),intent(in) :: modh
    logical,intent(in) :: pr
    integer :: i
    integer :: nhess
    integer,intent(in) :: natoms
    real(wp),intent(in) :: xyz(3,natoms)
    real(wp),intent(out) :: hess((natoms*3)*((natoms*3)+1)/2)
    integer,intent(in) :: at(natoms)

!>  initialize
    nhess = 3*natoms
    Hess = 0.0_wp

    select case (modh%model)
    case (0)
      if (pr) write (stdout,'(a)') "Using Lindh-Hessian (1995)"
      call ddvopt(xyz,natoms,Hess,at,modh)
!> other model hessians currently not tested
    case (1)
      if (pr) write (stdout,'(a)') "Using Lindh-Hessian"
      call mh_lindh_d2(xyz,natoms,Hess,at,modh)
    case (2)
      if (pr) write (stdout,'(a)') "Using Lindh-Hessian (2007)"
      call mh_lindh(xyz,natoms,Hess,at,modh)
    case (3)
      if (pr) write (stdout,'(a)') "Using Swart-Hessian"
      call mh_swart(xyz,natoms,Hess,at,modh)
    end select

!> add user-set constraint contributions to modelhessian
    call constrhess(natoms,at,xyz,calc,Hess)

    return
  end subroutine modhes

!========================================================================================!
  subroutine modhes_free(calc,modh,natoms,xyz,at,freezelist,hess,pr)
    implicit none
    type(calcdata),intent(in) :: calc
    type(mhparam),intent(in) :: modh
    logical,intent(in) :: pr
    integer,intent(in) :: natoms
    real(wp),intent(in) :: xyz(3,natoms)
    integer,intent(in) :: at(natoms)
    logical,intent(in) :: freezelist(natoms)
    real(wp),intent(out) :: hess(:)
    integer :: nfree3,np,i,j,p,ii,jj,q,atom,c,nfull3
    integer,allocatable :: idx(:)
    real(wp),allocatable :: work(:),fullh(:)

    nfree3 = 3*count(.not.freezelist)
    np = nfree3*(nfree3+1)/2
    if (size(hess) /= np) error stop 'modhes_free: inconsistent Hessian size'

    if (modh%model == 0) then
      allocate (work(0:np),source=0.0_wp)
      if (pr) write (stdout,'(a)') 'Using direct mobile/mobile Lindh-Hessian (1995)'
      call ddvopt_free(xyz,natoms,work,at,modh,freezelist)
      hess = work(1:np)
      deallocate (work)
    else
      ! Preserve the selected alternative model exactly. This uncommon fallback
      ! uses a temporary full Hessian, then extracts the mobile/mobile block.
      nfull3 = 3*natoms
      allocate (fullh(nfull3*(nfull3+1)/2),source=0.0_wp)
      call modhes(calc,modh,natoms,xyz,at,fullh,pr)
      call mobile_indices(freezelist,idx)
      call extract_mobile_block(fullh,idx,hess)
      deallocate (idx,fullh)
      return
    end if

    ! Explicit non-freeze constraints are uncommon in the frozen-host workflow.
    ! Keep their exact contribution through a full packed fallback while pure
    ! $fix-atoms jobs avoid that allocation completely.
    if (calc%nconstraints > 0) then
      nfull3 = 3*natoms
      allocate (fullh(nfull3*(nfull3+1)/2),source=0.0_wp)
      call constrhess(natoms,at,xyz,calc,fullh)
      call mobile_indices(freezelist,idx)
      allocate (work(np),source=0.0_wp)
      call extract_mobile_block(fullh,idx,work)
      hess = hess+work
      deallocate (work,idx,fullh)
    end if
  contains
    subroutine mobile_indices(frozen,indices)
      logical,intent(in) :: frozen(:)
      integer,allocatable,intent(out) :: indices(:)
      integer :: ia,ic,k
      allocate (indices(3*count(.not.frozen)))
      k=0
      do ia=1,size(frozen)
        if (.not.frozen(ia)) then
          do ic=1,3
            k=k+1
            indices(k)=3*(ia-1)+ic
          end do
        end if
      end do
    end subroutine mobile_indices

    subroutine extract_mobile_block(full,indices,reduced)
      real(wp),intent(in) :: full(:)
      integer,intent(in) :: indices(:)
      real(wp),intent(out) :: reduced(:)
      integer :: ii0,jj0,ii,jj,k,qq
      k=0
      do ii0=1,size(indices)
        ii=indices(ii0)
        do jj0=1,ii0
          jj=indices(jj0); k=k+1
          qq=max(ii,jj)*(max(ii,jj)-1)/2+min(ii,jj)
          reduced(k)=full(qq)
        end do
      end do
    end subroutine extract_mobile_block
  end subroutine modhes_free

!========================================================================================!
!########################################################################################!
!========================================================================================!
end module modelhessian_module
