!================================================================================!
! This file is part of crest.
!
! SPDX-License-Identifier: LGPL-3.0-or-later
!================================================================================!

module external_rerank_restart
  use crest_parameters,only:wp,stdout
  use crest_data,only:systemdata
  use crest_restartlog,only:restart_data
  use strucrd,only:coord,rdensemble,get_atlist
  implicit none
  private

  public :: validate_external_rerank_seed
  public :: fragment_aware_topology_equal
  public :: resolve_gff_fragment_ids

contains

  subroutine validate_external_rerank_seed(env,rdat,matched_frame)
!*******************************************************************************
!* Validate the externally selected iteration seed against the completed
!* iteration-1 ensemble. Coordinates may differ for every atom (mobile hosts are
!* supported). The invariant is the chemical system and its fragment-aware
!* covalent topology. User-provided GFN-FF fragments are taken from the TOML-
!* reconstructed calculation object and suppress inter-fragment topology edges.
!*******************************************************************************
    type(systemdata),intent(in) :: env
    type(restart_data),intent(in) :: rdat
    integer,intent(out),optional :: matched_frame

    type(coord) :: seed
    type(coord),allocatable :: refs(:)
    character(len=:),allocatable :: fragments(:)
    integer :: nall,i,match
    logical :: have_fragments,ok,ex

    if (len_trim(rdat%last_file) == 0) then
      error stop '**ERROR** external-rerank checkpoint has no iteration-1 archive'
    end if
    inquire(file=trim(rdat%last_file),exist=ex)
    if (.not.ex) then
      write(stdout,'(a,a)') '**ERROR** iteration-1 archive not found: ',trim(rdat%last_file)
      error stop '**ERROR** cannot validate external-rerank seed'
    end if

    call env%ref%to(seed)
    if (allocated(seed%lat)) then
      error stop '**ERROR** external-rerank seed must be a finite structure without lattice vectors'
    end if

    call rdensemble(trim(rdat%last_file),nall,refs)
    if (nall < 1) error stop '**ERROR** iteration-1 archive is empty'

    call get_consistent_gff_fragments(env,fragments,have_fragments)
    match = 0
    do i = 1,nall
      if (refs(i)%nat /= seed%nat) cycle
      if (any(refs(i)%at /= seed%at)) cycle
      if (have_fragments) then
        ok = fragment_aware_topology_equal(refs(i),seed,fragments)
      else
        ok = fragment_aware_topology_equal(refs(i),seed)
      end if
      if (ok) then
        match = i
        exit
      end if
    end do

    if (match == 0) then
      error stop '**ERROR** external-rerank seed topology does not match the iteration-1 ensemble'
    end if

    write(stdout,'(/,1x,a)') 'External-rerank seed validation:'
    write(stdout,'(3x,a,i0)') 'compatible topology reference frame: ',match
    if (have_fragments) then
      write(stdout,'(3x,a,i0)') 'GFN-FF fragment groups used: ',size(fragments)
      write(stdout,'(3x,a)') 'inter-fragment contacts were excluded from the topology comparison'
    else
      write(stdout,'(3x,a)') 'no user GFN-FF fragments: full-system topology compared'
    end if
    write(stdout,'(3x,a)') 'coordinate equality was not required'

    if (present(matched_frame)) matched_frame = match
    call seed%deallocate()
    do i = 1,size(refs)
      call refs(i)%deallocate()
    end do
    deallocate(refs)
    if (allocated(fragments)) deallocate(fragments)
  end subroutine validate_external_rerank_seed

  logical function fragment_aware_topology_equal(mol_a,mol_b,fragments) result(equal)
!*******************************************************************************
!* Compare covalent topology. When fragments are supplied, each user-defined
!* GFN-FF fragment (including implicit fragment 0) is analyzed independently.
!* This is stronger than calculating a full topology and masking cross edges:
!* close inter-fragment contacts cannot steal neighbours during topology setup.
!*******************************************************************************
    type(coord),intent(in) :: mol_a,mol_b
    character(len=*),intent(in),optional :: fragments(:)

    integer :: nat,ntopo,f,nf,nsub,i,k
    integer,allocatable :: topo_a(:),topo_b(:),frag_ids(:),subat(:),indices(:)
    real(wp),allocatable :: xyz_a(:,:),xyz_b(:,:)
    logical :: ok

    equal = .false.
    if (mol_a%nat /= mol_b%nat) return
    if (any(mol_a%at /= mol_b%at)) return

    nat = mol_a%nat
    if (.not.present(fragments)) then
      ntopo = nat*(nat+1)/2
      allocate(topo_a(ntopo),topo_b(ntopo),source=0)
      call quicktopo(nat,mol_a%at,mol_a%xyz,ntopo,topo_a)
      call quicktopo(nat,mol_b%at,mol_b%xyz,ntopo,topo_b)
      equal = all(topo_a == topo_b)
      deallocate(topo_b,topo_a)
      return
    end if

    call resolve_gff_fragment_ids(nat,mol_a%at,fragments,frag_ids,ok)
    if (.not.ok) return
    nf = maxval(frag_ids)
    do f = 0,nf
      nsub = count(frag_ids == f)
      if (nsub <= 1) cycle
      allocate(indices(nsub),subat(nsub),xyz_a(3,nsub),xyz_b(3,nsub))
      k = 0
      do i = 1,nat
        if (frag_ids(i) /= f) cycle
        k = k+1
        indices(k) = i
        subat(k) = mol_a%at(i)
        xyz_a(:,k) = mol_a%xyz(:,i)
        xyz_b(:,k) = mol_b%xyz(:,i)
      end do
      ntopo = nsub*(nsub+1)/2
      allocate(topo_a(ntopo),topo_b(ntopo),source=0)
      call quicktopo(nsub,subat,xyz_a,ntopo,topo_a)
      call quicktopo(nsub,subat,xyz_b,ntopo,topo_b)
      if (any(topo_a /= topo_b)) then
        deallocate(topo_b,topo_a,xyz_b,xyz_a,subat,indices,frag_ids)
        return
      end if
      deallocate(topo_b,topo_a,xyz_b,xyz_a,subat,indices)
    end do

    equal = .true.
    deallocate(frag_ids)
  end function fragment_aware_topology_equal

  subroutine resolve_gff_fragment_ids(nat,at,fragments,frag_ids,ok)
!*******************************************************************************
!* Resolve the TOML `fragments = [...]` strings with the same last-definition-
!* wins and implicit fragment-0 semantics used by gfnff_set_fragments().
!*******************************************************************************
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    character(len=*),intent(in) :: fragments(:)
    integer,allocatable,intent(out) :: frag_ids(:)
    logical,intent(out) :: ok

    logical,allocatable :: selected(:)
    integer :: f

    ok = .false.
    allocate(frag_ids(nat),source=0)
    do f = 1,size(fragments)
      call get_atlist(nat,selected,trim(fragments(f)),at)
      where(selected) frag_ids = f
      deallocate(selected)
    end do
    ok = .true.
  end subroutine resolve_gff_fragment_ids

  subroutine get_consistent_gff_fragments(env,fragments,found)
    type(systemdata),intent(in) :: env
    character(len=:),allocatable,intent(out) :: fragments(:)
    logical,intent(out) :: found

    integer :: i,j,first

    found = .false.
    first = 0
    if (.not.associated(env%calc)) return
    do i = 1,env%calc%ncalculations
      if (.not.allocated(env%calc%calcs(i)%gff_fragments)) cycle
      if (first == 0) then
        fragments = env%calc%calcs(i)%gff_fragments
        first = i
        found = .true.
      else
        if (size(fragments) /= size(env%calc%calcs(i)%gff_fragments)) then
          error stop '**ERROR** inconsistent GFN-FF fragment definitions across calculation levels'
        end if
        do j = 1,size(fragments)
          if (trim(fragments(j)) /= trim(env%calc%calcs(i)%gff_fragments(j))) then
            error stop '**ERROR** inconsistent GFN-FF fragment definitions across calculation levels'
          end if
        end do
      end if
    end do
  end subroutine get_consistent_gff_fragments

end module external_rerank_restart
