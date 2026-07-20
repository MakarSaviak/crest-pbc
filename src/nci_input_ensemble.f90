!================================================================================!
! This file is part of crest.
!
! SPDX-License-Identifier: LGPL-3.0-or-later
!================================================================================!

module nci_input_ensemble
  use crest_parameters, only: wp, stdout, bohr
  use crest_data, only: systemdata
  use crest_calculator, only: calcdata
  use optimize_module, only: optimize_geometry
  use strucrd, only: coord, rdensembleparam, rdensemble, checkcoordtype
  use axis_module, only: axis, cma
  implicit none
  private

  public :: load_nci_input_ensemble

  real(wp), parameter :: transform_tolerance = 1.0e-7_wp
  real(wp), parameter :: fixed_atom_tolerance = 1.0e-8_wp

contains

  subroutine load_nci_input_ensemble(env, prepared_ref, mols, ninputs)
!*******************************************************************************
!* Read a positional multi-frame XYZ for the first NCI MTD stage.
!*
!* The ordinary positional reader has already consumed, transformed and (when
!* requested) preoptimized frame 1.  We therefore reconstruct that same rigid
!* transformation from the raw first frame and apply it unchanged to every
!* other frame.  Aligning every host--guest placement independently would move
!* the fixed host and change the physical meaning of the supplied placements.
!*******************************************************************************
    type(systemdata), intent(inout) :: env
    type(coord), intent(in) :: prepared_ref
    type(coord), allocatable, intent(out) :: mols(:)
    integer, intent(out) :: ninputs

    type(coord) :: transformed_ref
    integer :: filetype, nat, natmax, iframe, iatom
    integer :: ntopo, nfixed, nelec
    integer, allocatable :: nats(:), ats(:,:), topo0(:), topoi(:)
    real(wp), allocatable :: xyz_ang(:,:,:), xyz_bohr(:,:,:)
    real(wp) :: rot(3), avmom, evec(3,3), center(3), det
    real(wp) :: shifted(3), max_fixed_delta, max_transform_delta
    logical :: is_xyz

    ninputs = 1
    if (.not. env%NCI) return
    if (.not. allocated(env%inputcoords)) return

    call checkcoordtype(trim(env%inputcoords), filetype)
    is_xyz = filetype == 2 ! strucrd's Xmol/XYZ file-type identifier
    if (.not. is_xyz) return

    call rdensembleparam(trim(env%inputcoords), natmax, ninputs)
    if (ninputs <= 1) then
      ninputs = 1
      return
    end if

    write (stdout,'(/,1x,a,i0,a)') &
      'Detected ', ninputs, ' positional XYZ frames for the first NCI MTD stage.'

    allocate (nats(ninputs), ats(natmax,ninputs))
    allocate (xyz_ang(3,natmax,ninputs), source=0.0_wp)
    nats = 0
    ats = 0
    call rdensemble(trim(env%inputcoords), natmax, ninputs, nats, ats, xyz_ang)

    nat = nats(1)
    if (nat /= prepared_ref%nat .or. nat /= env%ref%nat) then
      error stop 'Multi-input NCI: frame 1 atom count differs from the parsed input.'
    end if
    do iframe = 1,ninputs
      if (nats(iframe) /= nat) then
        error stop 'Multi-input NCI: every XYZ frame must have the same atom count.'
      end if
      if (any(ats(1:nat,iframe) /= ats(1:nat,1))) then
        error stop 'Multi-input NCI: element identities or atom ordering differ between frames.'
      end if
    end do
    if (any(ats(1:nat,1) /= prepared_ref%at)) then
      error stop 'Multi-input NCI: XYZ elements/order differ from the parsed reference.'
    end if

    call validate_global_selections(env, nat)
    nelec = sum(ats(1:nat,1)) - env%chrg
    if (nelec < 0 .or. mod(nelec-env%uhf,2) /= 0) then
      error stop 'Multi-input NCI: charge and unpaired-electron count are incompatible.'
    end if

    max_fixed_delta = 0.0_wp
    nfixed = 0
    if (allocated(env%calc%freezelist)) then
      nfixed = count(env%calc%freezelist)
      do iframe = 2,ninputs
        do iatom = 1,nat
          if (env%calc%freezelist(iatom)) then
            max_fixed_delta = max(max_fixed_delta, &
              maxval(abs(xyz_ang(:,iatom,iframe)-xyz_ang(:,iatom,1))))
          end if
        end do
      end do
      if (max_fixed_delta > fixed_atom_tolerance) then
        error stop 'Multi-input NCI: frozen-atom coordinates differ between XYZ frames.'
      end if
    end if

    ! Validate covalent topology in the raw input before any optimization.
    ntopo = nat*(nat+1)/2
    allocate (topo0(ntopo), topoi(ntopo))
    call quicktopo(nat, ats(1:nat,1), xyz_ang(:,1:nat,1)/bohr, ntopo, topo0)
    do iframe = 2,ninputs
      call quicktopo(nat, ats(1:nat,iframe), xyz_ang(:,1:nat,iframe)/bohr, &
        ntopo, topoi)
      if (any(topoi /= topo0)) then
        error stop 'Multi-input NCI: supplied XYZ frames have incompatible topology.'
      end if
    end do

    ! Reproduce inputcoords' frame-1 axis transform, then reuse it for all frames.
    call axis(nat, ats(1:nat,1), xyz_ang(:,1:nat,1), rot, avmom, evec)
    call cma(nat, ats(1:nat,1), xyz_ang(:,1:nat,1), center)
    det = determinant3(evec)
    if (det < 0.0_wp) evec(:,1) = -evec(:,1)

    allocate (xyz_bohr(3,nat,ninputs))
    do iframe = 1,ninputs
      do iatom = 1,nat
        shifted = xyz_ang(:,iatom,iframe) - center
        xyz_bohr(:,iatom,iframe) = matmul(transpose(evec),shifted)/bohr
      end do
    end do

    ! `coord` still contains transformed, unoptimized frame 1 after trialOPT.
    call transformed_ref%open('coord')
    max_transform_delta = maxval(abs(xyz_bohr(:,:,1)-transformed_ref%xyz))
    if (max_transform_delta > transform_tolerance) then
      error stop 'Multi-input NCI: failed to reproduce the common frame-1 transform.'
    end if
    call transformed_ref%deallocate()

    allocate (mols(ninputs))
    mols(1) = prepared_ref
    do iframe = 2,ninputs
      mols(iframe)%nat = nat
      mols(iframe)%at = ats(1:nat,iframe)
      mols(iframe)%xyz = xyz_bohr(:,:,iframe)
      mols(iframe)%chrg = env%chrg
      mols(iframe)%uhf = env%uhf
    end do

    if (env%preopt) call preoptimize_additional_inputs(env, mols, topo0)

    write (stdout,'(1x,a,i0)') 'Validated NCI input structures : ',ninputs
    write (stdout,'(1x,a,i0)') 'Frozen atoms checked           : ',nfixed
    if (nfixed > 0) then
      write (stdout,'(1x,a,es12.4,a)') &
        'Maximum frozen-atom deviation : ',max_fixed_delta,' Angstrom'
    end if
    write (stdout,'(1x,a,es12.4,a)') &
      'Common-transform check error   : ',max_transform_delta,' Bohr'

    deallocate (topoi,topo0,xyz_bohr,xyz_ang,ats,nats)
  end subroutine load_nci_input_ensemble

  subroutine validate_global_selections(env, nat)
    type(systemdata), intent(in) :: env
    integer, intent(in) :: nat

    if (.not. allocated(env%includeRMSD)) then
      error stop 'Multi-input NCI: RMSD/COM atom selection is not initialized.'
    end if
    if (size(env%includeRMSD) /= nat) then
      error stop 'Multi-input NCI: RMSD/COM atom selection has the wrong size.'
    end if
    if (count(env%includeRMSD /= 0) == 0) then
      error stop 'Multi-input NCI: RMSD/COM atom selection is empty.'
    end if
    if (allocated(env%calc%freezelist)) then
      if (size(env%calc%freezelist) /= nat) then
        error stop 'Multi-input NCI: fixed-atom selection has the wrong size.'
      end if
    end if
  end subroutine validate_global_selections

  subroutine preoptimize_additional_inputs(env, mols, input_topology)
!*******************************************************************************
!* Apply the same loose initial optimization used by trialOPT to frames 2..N.
!* GFN-FF is deliberately kept at one thread per calculation.  These few input
!* optimizations are sequential to avoid shared calculator/log/scratch state;
!* the expensive N x bias MTD batch remains globally parallel.
!*******************************************************************************
    type(systemdata), intent(inout) :: env
    type(coord), intent(inout) :: mols(:)
    integer, intent(in) :: input_topology(:)

    type(coord) :: molopt
    type(calcdata) :: tmpcalc
    real(wp), allocatable :: grad(:,:)
    integer, allocatable :: topo(:)
    real(wp) :: energy
    integer :: iframe, io, ntopo, outer_threads, inner_threads
    logical :: pr, wr

    ntopo = size(input_topology)
    allocate (grad(3,mols(1)%nat), topo(ntopo))
    pr = .false.
    wr = .false.

    call new_ompautoset(env,'serial',0,outer_threads,inner_threads)
    write (stdout,'(/,1x,a,i0,a)') 'Preoptimizing ',size(mols)-1, &
      ' additional NCI input structures with one GFN-FF core each.'

    do iframe = 2,size(mols)
      tmpcalc = env%calc
      tmpcalc%optlev = -1
      molopt = mols(iframe)
      grad = 0.0_wp
      call optimize_geometry(mols(iframe),molopt,tmpcalc,energy,grad,pr,wr,io)
      if (io /= 0) then
        error stop 'Multi-input NCI: initial optimization failed for an input frame.'
      end if
      call quicktopo(molopt%nat,molopt%at,molopt%xyz,ntopo,topo)
      if (any(topo /= input_topology)) then
        error stop 'Multi-input NCI: input topology changed during initial optimization.'
      end if
      mols(iframe) = molopt
      call molopt%deallocate()
      write (stdout,'(3x,a,i0,a)') 'Input structure ',iframe,' successfully optimized.'
    end do

    call new_ompautoset(env,'max',0,outer_threads,inner_threads)
    deallocate (topo,grad)
  end subroutine preoptimize_additional_inputs

  pure real(wp) function determinant3(matrix)
    real(wp), intent(in) :: matrix(3,3)
    determinant3 = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)) &
      + matrix(1,2)*(matrix(2,3)*matrix(3,1)-matrix(2,1)*matrix(3,3)) &
      + matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
  end function determinant3

end module nci_input_ensemble
