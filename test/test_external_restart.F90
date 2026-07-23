!================================================================================!
! This file is part of crest.
! SPDX-License-Identifier: LGPL-3.0-or-later
!================================================================================!

module test_external_restart
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters,only:wp,bohr
  use strucrd,only:coord
  use external_rerank_restart,only:fragment_aware_topology_equal, &
    & resolve_gff_fragment_ids
  implicit none
  private
  public :: collect_external_restart

contains

  subroutine collect_external_restart(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest('external restart fragment topology',test_fragment_topology), &
      new_unittest('external restart internal bond change',test_internal_change), &
      new_unittest('external restart fragment resolution',test_fragment_resolution) &
    ]
  end subroutine collect_external_restart

  subroutine setup_dimers(a,b)
    type(coord),intent(out) :: a,b
    real(wp) :: xa(3,4),xb(3,4)
    integer :: at(4)
    at = 6
    xa = 0.0_wp
    xa(1,:) = [0.0_wp,1.40_wp,6.00_wp,7.40_wp]
    xb = xa
    ! Move the second fragment close enough to create inter-fragment contacts.
    ! Its internal C-C distance stays 1.40 Angstrom.
    xb(1,3:4) = [2.70_wp,4.10_wp]
    a%nat=4; a%at=at; a%xyz=xa/bohr
    b%nat=4; b%at=at; b%xyz=xb/bohr
  end subroutine setup_dimers

  subroutine test_fragment_topology(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: a,b
    character(len=3) :: fragments(2)
    logical :: with_fragments,without_fragments
    fragments = ['1-2','3-4']
    call setup_dimers(a,b)
    with_fragments = fragment_aware_topology_equal(a,b,fragments)
    without_fragments = fragment_aware_topology_equal(a,b)
    call check(error,with_fragments,.true.)
    if (allocated(error)) return
    call check(error,without_fragments,.false.)
  end subroutine test_fragment_topology

  subroutine test_internal_change(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: a,b
    character(len=3) :: fragments(2)
    logical :: same
    fragments = ['1-2','3-4']
    call setup_dimers(a,b)
    ! Break the bond inside fragment 2. Inter-fragment masking must not hide it.
    b%xyz(1,4) = 6.20_wp/bohr
    same = fragment_aware_topology_equal(a,b,fragments)
    call check(error,same,.false.)
  end subroutine test_internal_change

  subroutine test_fragment_resolution(error)
    type(error_type),allocatable,intent(out) :: error
    integer,allocatable :: ids(:)
    integer :: at(6)
    character(len=3) :: fragments(2)
    logical :: ok
    at = [8,1,1,8,1,1]
    fragments = ['1-3','4-6']
    call resolve_gff_fragment_ids(6,at,fragments,ids,ok)
    call check(error,ok,.true.)
    if (allocated(error)) return
    if (any(ids /= [1,1,1,2,2,2])) then
      call test_failed(error,'unexpected resolved GFN-FF fragment IDs')
    end if
  end subroutine test_fragment_resolution

end module test_external_restart
