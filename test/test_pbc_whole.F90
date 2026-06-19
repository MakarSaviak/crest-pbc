module test_pbc_whole
  use iso_fortran_env,only:wp => real64
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use metadynamics_module,only:read_whole_bond_graph,pbc_make_whole
  use iomod,only:remove
  implicit none
  private

  public :: collect_pbc_whole

contains

  subroutine collect_pbc_whole(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("selected lowercase bond parser",test_selected_bond_parser), &
      new_unittest("disconnected graph rejected",test_disconnected_graph), &
      new_unittest("selected PBC make-whole",test_selected_make_whole) &
    ]
  end subroutine collect_pbc_whole

  subroutine test_selected_bond_parser(error)
    type(error_type),allocatable,intent(out) :: error
    logical :: mask(5)
    integer,allocatable :: bonds(:,:)
    integer :: ich,io,nbonds,ncomponents
    character(len=*),parameter :: fname = 'test_pbcwhole_bondlengths'

    mask = [.false.,.false.,.true.,.true.,.true.]
    open (newunit=ich,file=fname,status='replace',action='write')
    write (ich,'(a)') '$constrain'
    write (ich,'(a)') ' distance: 1, 2, 1.0'
    write (ich,'(a)') ' distance: 3, 4, 1.0'
    write (ich,'(a)') ' distance: 4, 5, 1.0'
    write (ich,'(a)') ' DISTANCE: 3, 5, 1.0, 4.0'
    write (ich,'(a)') '$end'
    close (ich)

    call read_whole_bond_graph(fname,5,mask,bonds,nbonds,ncomponents,io)
    call remove(fname)
    call check(error,io,0)
    if (allocated(error)) return
    call check(error,nbonds,2)
    if (allocated(error)) return
    call check(error,ncomponents,1)
    if (allocated(error)) return
    if (any(.not.mask(bonds(1,:))).or.any(.not.mask(bonds(2,:)))) then
      call test_failed(error,'Parser retained an unselected bond endpoint')
    end if
  end subroutine test_selected_bond_parser

  subroutine test_disconnected_graph(error)
    type(error_type),allocatable,intent(out) :: error
    logical :: mask(5)
    integer,allocatable :: bonds(:,:)
    integer :: ich,io,nbonds,ncomponents
    character(len=*),parameter :: fname = 'test_pbcwhole_disconnected'

    mask = [.false.,.false.,.true.,.true.,.true.]
    open (newunit=ich,file=fname,status='replace',action='write')
    write (ich,'(a)') '$constrain'
    write (ich,'(a)') ' distance: 3, 4, 1.0'
    write (ich,'(a)') '$end'
    close (ich)

    call read_whole_bond_graph(fname,5,mask,bonds,nbonds,ncomponents,io)
    call remove(fname)
    call check(error,io,3)
    if (allocated(error)) return
    call check(error,ncomponents,2)
  end subroutine test_disconnected_graph

  subroutine test_selected_make_whole(error)
    type(error_type),allocatable,intent(out) :: error
    integer,parameter :: nat = 5,nbonds = 2
    logical :: mask(nat)
    integer :: bonds(2,nbonds)
    real(wp) :: lat(3,3),whole_xyz(3,nat),split_xyz(3,nat),result(3,nat)
    real(wp),parameter :: tol = 1.0e-12_wp

    mask = [.false.,.false.,.true.,.true.,.true.]
    bonds = reshape([3,4,4,5],shape(bonds))
    lat = 0.0_wp
    lat(1,1) = 10.0_wp
    lat(2,2) = 10.0_wp
    lat(3,3) = 10.0_wp
    whole_xyz = 0.0_wp
    whole_xyz(:,1) = [2.0_wp,2.0_wp,2.0_wp]
    whole_xyz(:,2) = [8.0_wp,8.0_wp,8.0_wp]
    whole_xyz(:,3) = [4.5_wp,5.0_wp,5.0_wp]
    whole_xyz(:,4) = [5.3_wp,5.0_wp,5.0_wp]
    whole_xyz(:,5) = [6.1_wp,5.0_wp,5.0_wp]
    split_xyz = whole_xyz
    split_xyz(:,5) = split_xyz(:,5)-lat(:,1)

    call pbc_make_whole(nat,split_xyz,lat,mask,bonds,nbonds,result)
    if (any(abs(result(:,3:5)-whole_xyz(:,3:5)) > tol)) then
      call test_failed(error,'Selected atoms were not reconstructed correctly')
      return
    end if
    if (any(result(:,1:2) /= split_xyz(:,1:2))) then
      call test_failed(error,'Unselected atoms were modified')
    end if
  end subroutine test_selected_make_whole

end module test_pbc_whole
