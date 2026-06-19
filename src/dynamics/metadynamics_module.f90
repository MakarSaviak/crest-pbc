!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2023 Philipp Pracht
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

module metadynamics_module

  use crest_parameters
  use ls_rmsd
  use strucrd

  implicit none

  !======================================================================================!
  !--- private module variables and parameters
  private
  integer,parameter,public :: cv_std_mtd = 1
  integer,parameter,public :: cv_rmsd = 2
  integer,parameter :: damp_heaviside = 3
  integer,parameter :: damp_heaviside_cv = 4
  integer,parameter,public :: cv_rmsd_static = 5

  logical,save :: whole_graph_ready = .false.
  integer,save :: cached_whole_nbonds = 0
  logical,allocatable,save :: cached_whole_mask(:)
  integer,allocatable,save :: cached_whole_bonds(:,:)

  !======================================================================================!
  !data object that contains settings and trackers for a single MTD potential
  type :: mtdpot

    integer :: mtdtype = 0
    integer :: nmax = 0
    integer :: ncur = 0

    real(wp) :: kpush = 1.0_wp
    real(wp) :: alpha = 1.0_wp

    !>--- regular MTD, some CV
    real(wp),allocatable :: cv(:) !list of cv values at each timestep
    real(wp),allocatable :: cvgrd(:,:) !Cartesian gradient of CV at step t

    !>--- RMSD MTD
    integer :: cvdump = 0 !xyz dump counter
    real(wp) :: cvdump_fs = 0.0_wp !xyz dump frequency (in fs)
    integer :: cvdumpstep = 0 !xyz dump frequency (in MD steps)
    integer :: maxsave = 0
    character(len=:),allocatable :: biasfile !> specify a file from which the bias is obtained
    logical,allocatable :: atinclude(:) !specify atoms to include in RMSD potentail
    real(wp),allocatable :: cvxyz(:,:,:) !ensemble of CV structures to calculate RMSD from
    logical :: whole = .false. !> repair selected PBC atoms before RMSD-MTD
    integer :: whole_nbonds = 0
    logical,allocatable :: whole_mask(:)
    integer,allocatable :: whole_bonds(:,:) !> selected-atom bond graph for make-whole

    !>--- damping of the MTD potential
    integer :: damptype = 0
    real(wp) :: ramp = -1.0_wp
    real(wp) :: damp = 1.0_wp
    real(wp),allocatable :: damping(:)  !input for snapshot-specific damping

  contains
    procedure :: deallocate => mtd_deallocate
    procedure :: info => mtd_info
  end type mtdpot

  public :: mtdpot
  public :: mtd_ini
  public :: cv_dump
  public :: calc_mtd
  public :: prepare_whole_bond_graph
  public :: read_whole_bond_graph
  public :: pbc_make_whole

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine mtd_ini(mol,pot,tstep,mdlength,pr)
!*************************************
!* subroutine mtd_ini
!* initialize metadynamics settings
!*************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(in) :: tstep !> MD timestep in fs
    real(wp),intent(in) :: mdlength !> MD length in ps
    logical,intent(in) :: pr
    real(wp) :: dum1,dum2
    integer :: i,j
    integer :: idum1,idum2,nall,nat
    logical :: ex
    integer,allocatable :: at(:)

    if (pr) then
      write (stdout,'(1X,17("─"))', advance='no')
      write (stdout,'(1X,"Metadynamics Parameters")', advance='no')
      write (stdout,'(1X,17("─"))')
    end if

    dum1 = anint((mdlength*1000.0_wp)/tstep)
    idum1 = nint(dum1)

    select case (pot%mtdtype)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case (cv_std_mtd) !>--- "standard" MTD
      pot%nmax = idum1
      allocate (pot%cv(idum1),source=0.0_wp)
      allocate (pot%cvgrd(3,mol%nat),source=0.0_wp)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case (cv_rmsd) !>--- RMSD MTD
      if (pot%cvdump_fs <= 0.0_wp) return !> structure dumpstep in fs must be given
      dum2 = max(1.0_wp, (pot%cvdump_fs/tstep))
      dum2 = anint(dum2)
      pot%cvdumpstep = nint(dum2) !> bias structure dump step in MD steps
      dum1 = dum1/dum2
      dum1 = floor(dum1)
      pot%nmax = nint(dum1) !> max number of bias structure dump
      if (pot%maxsave == 0) pot%maxsave = nint(dum1)
      if (allocated(pot%cvxyz)) deallocate (pot%cvxyz)
      allocate (pot%cvxyz(3,mol%nat,pot%maxsave),source=0.0_wp)
      !>--- automatic ramp parameter (acounted for both different MD time steps and CV dump steps)
      !> (should yield damp≈0.5 for cvdumpstep/2, but is at least 0.03 as in xtb)
      if (pot%ramp <= 0.0_wp) then !> only if not set by the user
        pot%ramp = log(3.0_wp)/(0.5_wp*dum2)
        pot%ramp = max(pot%ramp,0.03_wp)
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case (cv_rmsd_static) !> static RMSD MTD / umbrella sampling with RMSD pot
      pot%cvdumpstep = huge(idum1) !> we won't modify the potential, so set the dumstep to "infinity"
      pot%cvdump_fs = huge(dum1)   !> same for _fs version
      if (allocated(pot%cvxyz)) then !> if structures were already loaded, just determine the rest
        nat = size(pot%cvxyz,2)
        nall = size(pot%cvxyz,3)
      else if (allocated(pot%biasfile)) then !> else try to read from file
        inquire (file=pot%biasfile,exist=ex)
        if (.not.ex) return !> if the file is absent, we return
        call rdensembleparam(pot%biasfile,nat,nall)
        if (nat .ne. 0.and.nall .ne. 0) then
          allocate (pot%cvxyz(3,nat,nall),source=0.0_wp)
          allocate (at(nat),source=0)
          call rdensemble(pot%biasfile,nat,nall,at,pot%cvxyz)
          deallocate (at)
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
          !>--- Important: if we read here we must convert to Bohrs
          pot%cvxyz = pot%cvxyz*aatoau
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        end if
      else  !> if both failed we must return (the potential is not set up)
        return
      end if
      if (nat .ne. mol%nat) then !> can't do that! something is wrong
        if (allocated(pot%cvxyz)) deallocate (pot%cvxyz)
        pot%mtdtype = 0
        write (stdout,'(1x,a)') '*WARNING* static metadynamics setup failed! Mismatch of #atoms'
        !return
        error stop
      end if
      !> for safety, perturb all bias slightly (so the potential won't explode)
      do i = 1,nall
        call rmsdcv_perturb(nat,pot%cvxyz(:,:,i))
      end do
      pot%ncur = nall    !> will not change
      pot%maxsave = nall !> won't change either
      if (pot%ramp <= 0.0_wp) then          !> only if not set by the user
        pot%ramp = (tstep/5.0_wp)*0.015_wp  !> a default derived from the entropy mode at GFN2-xTB
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case default
      return

    end select

    call setup_mtd_whole(mol,pot,pr)

    !>--- printout
    if (pr) then
      call pot%info(stdout)
      write (stdout,'(1X,59("─"))')
    end if

    return
  end subroutine mtd_ini

!========================================================================================!
  subroutine setup_mtd_whole(mol,pot,pr)
!**********************************************
!* Set up PBC make-whole support for RMSD-MTD.
!* The selected graph is read once from the
!* freshly generated bondlengths file.
!**********************************************
    implicit none
    type(coord),intent(in) :: mol
    type(mtdpot),intent(inout) :: pot
    logical,intent(in) :: pr
    if (.not.pot%whole) return
    if (.not.(pot%mtdtype == cv_rmsd.or.pot%mtdtype == cv_rmsd_static)) return

    if (.not.allocated(mol%lat)) then
      write (stdout,'(1x,a)') 'ERROR: -whole requested, but no lattice is available.'
      error stop
    end if

    if (.not.whole_graph_ready) then
      write (stdout,'(1x,a)') 'ERROR: selected PBC make-whole graph was not prepared before MD.'
      error stop
    end if
    if (size(cached_whole_mask) /= mol%nat) then
      write (stdout,'(1x,a)') 'ERROR: selected PBC make-whole graph has the wrong system size.'
      error stop
    end if

    if (.not.allocated(pot%atinclude)) then
      allocate (pot%atinclude(mol%nat),source=cached_whole_mask)
    else if (size(pot%atinclude) /= mol%nat.or. &
    & any(pot%atinclude .neqv. cached_whole_mask)) then
      write (stdout,'(1x,a)') 'ERROR: RMSD-MTD selection differs from the prepared PBC graph.'
      error stop
    end if

    if (allocated(pot%whole_mask)) deallocate (pot%whole_mask)
    allocate (pot%whole_mask(mol%nat),source=cached_whole_mask)
    if (allocated(pot%whole_bonds)) deallocate (pot%whole_bonds)
    allocate (pot%whole_bonds(2,cached_whole_nbonds),source=cached_whole_bonds)
    pot%whole_nbonds = cached_whole_nbonds
    if (pr) write (stdout,'(1x,a,i0,a,i0,a)') 'Loaded ',count(pot%whole_mask,1), &
    & ' selected atoms and ',pot%whole_nbonds,' bonds for PBC make-whole.'

    if (pot%mtdtype == cv_rmsd_static .and. allocated(pot%cvxyz)) then
      block
        integer :: i
        real(wp),allocatable :: xyzwhole(:,:)
        allocate (xyzwhole(3,size(pot%cvxyz,2)),source=0.0_wp)
        do i = 1,size(pot%cvxyz,3)
          call pbc_make_whole(size(pot%cvxyz,2),pot%cvxyz(:,:,i),mol%lat, &
          & pot%whole_mask,pot%whole_bonds,pot%whole_nbonds,xyzwhole)
          pot%cvxyz(:,:,i) = xyzwhole(:,:)
        end do
        deallocate (xyzwhole)
      end block
    end if
  end subroutine setup_mtd_whole

!========================================================================================!
  subroutine prepare_whole_bond_graph(fname,nat,mask,nbonds,ncomponents,iostat)
!***********************************************************************
!* Parse and cache the selected graph once, before any MTD worker can  *
!* change directories. Each mtdpot receives its own immutable copy.    *
!***********************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    logical,intent(in) :: mask(nat)
    integer,intent(out) :: nbonds,ncomponents,iostat
    integer,allocatable :: bonds(:,:)

    whole_graph_ready = .false.
    cached_whole_nbonds = 0
    if (allocated(cached_whole_mask)) deallocate (cached_whole_mask)
    if (allocated(cached_whole_bonds)) deallocate (cached_whole_bonds)

    call read_whole_bond_graph(fname,nat,mask,bonds,nbonds,ncomponents,iostat)
    if (iostat /= 0) return
    allocate (cached_whole_mask(nat),source=mask)
    allocate (cached_whole_bonds(2,nbonds),source=bonds)
    cached_whole_nbonds = nbonds
    whole_graph_ready = .true.
    deallocate (bonds)
  end subroutine prepare_whole_bond_graph

!========================================================================================!
  subroutine read_whole_bond_graph(fname,nat,mask,bonds,nbonds,ncomponents,iostat)
!***********************************************************************
!* Read the lower-case distance records from a bondlengths-style file. *
!* Only pairs fully contained in the supplied RMSD-MTD mask are kept.  *
!***********************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    logical,intent(in) :: mask(nat)
    integer,allocatable,intent(inout) :: bonds(:,:)
    integer,intent(out) :: nbonds,ncomponents,iostat

    character(len=1024) :: line
    integer :: i,j,k,ich,io,nmax,itmp
    logical :: duplicate
    integer,allocatable :: tmp(:,:)

    nbonds = 0
    ncomponents = 0
    iostat = 0
    if (allocated(bonds)) deallocate (bonds)
    if (count(mask,1) < 2) then
      iostat = 4
      return
    end if

    nmax = count(mask,1)*(count(mask,1)-1)/2
    allocate (tmp(2,nmax),source=0)

    open (newunit=ich,file=fname,status='old',action='read',iostat=io)
    if (io /= 0) then
      iostat = 1
      deallocate (tmp)
      return
    end if
    do
      read (ich,'(a)',iostat=io) line
      if (io /= 0) exit
      line = adjustl(line)
      if (index(line,'distance:') /= 1) cycle
      read (line(10:),*,iostat=io) i,j
      if (io /= 0) cycle
      if (i < 1.or.i > nat.or.j < 1.or.j > nat.or.i == j) cycle
      if (.not.mask(i).or..not.mask(j)) cycle
      if (i > j) then
        itmp = i
        i = j
        j = itmp
      end if
      duplicate = .false.
      do k = 1,nbonds
        if (tmp(1,k) == i.and.tmp(2,k) == j) then
          duplicate = .true.
          exit
        end if
      end do
      if (duplicate) cycle
      nbonds = nbonds+1
      tmp(1,nbonds) = i
      tmp(2,nbonds) = j
    end do
    close (ich)

    if (nbonds < 1) then
      iostat = 2
      deallocate (tmp)
      return
    end if
    call count_whole_components(nat,mask,tmp(:,1:nbonds),nbonds,ncomponents)
    if (ncomponents /= 1) then
      iostat = 3
      deallocate (tmp)
      return
    end if
    allocate (bonds(2,nbonds),source=tmp(:,1:nbonds))
    deallocate (tmp)
  end subroutine read_whole_bond_graph

!========================================================================================!
  subroutine count_whole_components(nat,mask,bonds,nbonds,ncomponents)
    implicit none
    integer,intent(in) :: nat,nbonds
    logical,intent(in) :: mask(nat)
    integer,intent(in) :: bonds(2,nbonds)
    integer,intent(out) :: ncomponents

    logical,allocatable :: seen(:)
    integer,allocatable :: queue(:)
    integer :: seed,head,tail,i,j,b

    ncomponents = 0
    allocate (seen(nat),source=.false.)
    allocate (queue(nat),source=0)
    do seed = 1,nat
      if (.not.mask(seed).or.seen(seed)) cycle
      ncomponents = ncomponents+1
      head = 1
      tail = 1
      queue(tail) = seed
      seen(seed) = .true.
      do while (head <= tail)
        i = queue(head)
        head = head+1
        do b = 1,nbonds
          j = 0
          if (bonds(1,b) == i) then
            j = bonds(2,b)
          else if (bonds(2,b) == i) then
            j = bonds(1,b)
          end if
          if (j < 1.or.seen(j)) cycle
          seen(j) = .true.
          tail = tail+1
          queue(tail) = j
        end do
      end do
    end do
    deallocate (queue,seen)
  end subroutine count_whole_components

!========================================================================================!
  subroutine pbc_make_whole(nat,xyz,lat,mask,bonds,nbonds,whole)
!**********************************************
!* Reconstruct selected atoms by walking the
!* fixed bond graph with minimum-image bonded
!* displacements. Unselected atoms are copied.
!**********************************************
    implicit none
    integer,intent(in) :: nat,nbonds
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: lat(3,3)
    logical,intent(in) :: mask(nat)
    integer,intent(in) :: bonds(2,nbonds)
    real(wp),intent(out) :: whole(3,nat)

    logical,allocatable :: seen(:)
    integer,allocatable :: queue(:)
    integer :: head,tail,seed,i,j,b
    real(wp) :: invlat(3,3),dx(3),mic(3)

    whole(:,:) = xyz(:,:)
    if (nbonds < 1.or.count(mask,1) < 2) return

    call invert_lat3(lat,invlat)
    allocate (seen(nat),source=.false.)
    allocate (queue(nat),source=0)

    do seed = 1,nat
      if (.not.mask(seed).or.seen(seed)) cycle
      head = 1
      tail = 1
      queue(tail) = seed
      seen(seed) = .true.
      whole(:,seed) = xyz(:,seed)

      do while (head <= tail)
        i = queue(head)
        head = head+1
        do b = 1,nbonds
          j = 0
          if (bonds(1,b) == i) then
            j = bonds(2,b)
          else if (bonds(2,b) == i) then
            j = bonds(1,b)
          end if
          if (j < 1.or.seen(j)) cycle
          dx(:) = xyz(:,j)-xyz(:,i)
          call minimum_image_displacement(dx,lat,invlat,mic)
          whole(:,j) = whole(:,i)+mic(:)
          seen(j) = .true.
          tail = tail+1
          queue(tail) = j
        end do
      end do
    end do

    deallocate (queue,seen)
  end subroutine pbc_make_whole

!========================================================================================!
  subroutine minimum_image_displacement(dx,lat,invlat,mic)
    implicit none
    real(wp),intent(in) :: dx(3),lat(3,3),invlat(3,3)
    real(wp),intent(out) :: mic(3)
    real(wp) :: df(3)
    integer :: k
    df(:) = matmul(invlat,dx)
    do k = 1,3
      df(k) = df(k)-anint(df(k))
    end do
    mic(:) = matmul(lat,df)
  end subroutine minimum_image_displacement

!========================================================================================!
  subroutine invert_lat3(lat,invlat)
    implicit none
    real(wp),intent(in) :: lat(3,3)
    real(wp),intent(out) :: invlat(3,3)
    real(wp) :: det

    det = lat(1,1)*(lat(2,2)*lat(3,3)-lat(2,3)*lat(3,2)) &
    &   - lat(1,2)*(lat(2,1)*lat(3,3)-lat(2,3)*lat(3,1)) &
    &   + lat(1,3)*(lat(2,1)*lat(3,2)-lat(2,2)*lat(3,1))

    if (abs(det) <= epsilon(det)) then
      invlat = 0.0_wp
      return
    end if

    invlat(1,1) =  (lat(2,2)*lat(3,3)-lat(2,3)*lat(3,2))/det
    invlat(1,2) = -(lat(1,2)*lat(3,3)-lat(1,3)*lat(3,2))/det
    invlat(1,3) =  (lat(1,2)*lat(2,3)-lat(1,3)*lat(2,2))/det
    invlat(2,1) = -(lat(2,1)*lat(3,3)-lat(2,3)*lat(3,1))/det
    invlat(2,2) =  (lat(1,1)*lat(3,3)-lat(1,3)*lat(3,1))/det
    invlat(2,3) = -(lat(1,1)*lat(2,3)-lat(1,3)*lat(2,1))/det
    invlat(3,1) =  (lat(2,1)*lat(3,2)-lat(2,2)*lat(3,1))/det
    invlat(3,2) = -(lat(1,1)*lat(3,2)-lat(1,2)*lat(3,1))/det
    invlat(3,3) =  (lat(1,1)*lat(2,2)-lat(1,2)*lat(2,1))/det
  end subroutine invert_lat3

!========================================================================================!
  subroutine mtd_deallocate(self)
!**********************************************
!* subroutine mtd_deallocate
!* type internal procedure to deallocate data
!**********************************************
    class(mtdpot) :: self
    if (allocated(self%cvxyz)) deallocate (self%cvxyz)
    if (allocated(self%atinclude)) deallocate (self%atinclude)
    if (allocated(self%whole_mask)) deallocate (self%whole_mask)
    if (allocated(self%whole_bonds)) deallocate (self%whole_bonds)
    if (allocated(self%cv)) deallocate (self%cv)
    if (allocated(self%cvgrd)) deallocate (self%cvgrd)

    self%mtdtype = 0
    self%nmax = 0
    self%ncur = 0
    self%kpush = 1.0_wp
    self%alpha = 1.0_wp
    self%cvdump = 0 !xyz dump counter
    self%cvdump_fs = 0.0_wp !xyz dump frequency (in fs)
    self%cvdumpstep = 0 !xyz dump frequency (in MD steps)
    self%maxsave = 0
    self%whole = .false.
    self%whole_nbonds = 0
    self%damptype = 0
    self%ramp = -1.0_wp
    self%damp = 1.0_wp

    return
  end subroutine mtd_deallocate

!========================================================================================!
  subroutine mtd_info(self,iunit)
!**********************************************
!* subroutine mtd_info
!* print information about the MTD potential
!**********************************************
    class(mtdpot) :: self
    integer,intent(in) :: iunit

    !write (iunit,'(" --- metadynamics parameter ---")')
    select case (self%mtdtype)
    case (cv_std_mtd)
      write (iunit,'("  MTD/CV type",t25,":",1x,a)') 'standard'
    case (cv_rmsd)
      write (stdout,'("  MTD/CV type",t25,":",1x,a)') 'RMSD bias'
    case (cv_rmsd_static)
      write (iunit,'("  MTD/CV type",t25,":",1x,a)') 'RMSD bias (static)'
    end select
    write (iunit,'("  kpush /Eh",t25,":",f10.4)') self%kpush
    write (iunit,'("  alpha /Bohr⁻²",t28,":",f10.4)') self%alpha

    select case (self%mtdtype)
    case (cv_rmsd)
      write (iunit,'("  ramp rate",t25,":",f10.4,1x,"(",i0,")")') self%ramp,check_dump_steps_rmsd(self)
      write (iunit,'("  dump/fs",t25,":",f10.4,1x,"(",i0,")")') self%cvdump_fs,self%cvdumpstep
      write (iunit,'("  # CVs (max)",t25,":",i10 )') self%maxsave
    case (cv_rmsd_static)
      if (allocated(self%biasfile)) write (iunit,'(" reading from",t25,":",1x,a)') self%biasfile
      write (iunit,'("  ramp (adjust.)",t25,":",f10.4,1x,i0)') self%ramp,check_dump_steps_rmsd(self)
      write (iunit,'("  # CVs (loaded)",t25,":",i10 )') self%maxsave
    end select
    if (self%mtdtype == cv_rmsd.or.self%mtdtype == cv_rmsd_static) then
      if (allocated(self%atinclude)) then
        write (iunit,'("  # of atoms affected",t25,":",i10)') count(self%atinclude,1)
      end if
      if (self%whole) then
        write (iunit,'("  PBC make-whole",t25,":",1x,a,1x,"(",i0," bonds)")') &
        & 'on',self%whole_nbonds
      else
        write (iunit,'("  PBC make-whole",t25,":",1x,a)') 'off'
      end if
    end if

    return
  end subroutine mtd_info

!========================================================================================!
  subroutine cv_dump(mol,pot,cv,pr)
!*****************************************************
!* subroutine cv_dump
!* update the list of CVs at the current MD timestep
!*****************************************************
!$  use omp_lib
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(in) :: cv
    logical :: pr

    select case (pot%mtdtype)
    case (cv_std_mtd) !>--- CV update for "standard" MTD
      pot%ncur = pot%ncur+1
      pot%cv(pot%ncur) = cv

    case (cv_rmsd) !>--- structure mapping for RMSD MTD
      pot%cvdump = pot%cvdump+1  !> cvdump counts the MD step since the last CV was added
      if (pot%cvdump == pot%cvdumpstep) then !> the MTD tracks when it needs to be updated
        pot%cvdump = 0  !> reset if new CV is added
        pot%ncur = pot%ncur+1
        if (pot%whole .and. allocated(mol%lat) .and. allocated(pot%whole_bonds)) then
          call pbc_make_whole(mol%nat,mol%xyz,mol%lat,pot%whole_mask, &
          & pot%whole_bonds,pot%whole_nbonds,pot%cvxyz(:,:,pot%ncur))
        else
          pot%cvxyz(:,:,pot%ncur) = mol%xyz(:,:)
        end if
        if (pot%ncur == 1) then
          !>--- The first one should be sligthly distorted
          call rmsdcv_perturb(mol%nat,pot%cvxyz(:,:,pot%ncur))
        end if
        if (pr) then
          write (stdout,'(2x,"adding snapshot to metadynamics bias, now at ",i0," CVs")') pot%ncur
        end if
      end if

    case (cv_rmsd_static)
      pot%cvdump = pot%cvdump+1 !> the cvdump is equal to the MD step
      !> no further update necessary

    case default
      return
    end select

    return

  end subroutine cv_dump

!=========================================================================================!
  subroutine rmsdcv_perturb(nat,xyz)
!************************************************
!* Slightly perturb a given geometry for RMSD CV
!* to avoid singularities if the CV is exactly the
!* current structure
!*************************************************
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp) :: r(3)
    integer :: i,j
    real(wp),parameter :: tol = 1.0e-8_wp
    real(wp),parameter :: displace = 1.0e-6_wp
    do i = 1,nat
      do
        !> generate a random vector r in [-1,1]
        call random_number(r)
        r = (r-0.5_wp)*2.0_wp
        !> check that displacement is large enough
        if (norm2(r) >= 1e-8_wp) exit
      end do
      !> normalize
      r = r/norm2(r)
      !> displace
      xyz(:,i) = xyz(:,i)+displace*r
    end do
  end subroutine rmsdcv_perturb

!========================================================================================!
  subroutine calc_mtd(mol,pot,emtd,grdmtd)
!*******************************************************
!* subroutine calc_mtd
!* select how the MTD potential is calculated.
!* On Output:
!*             emtd - final MTD energy contribution
!*           mtdgrd - final MTD gradient contribution
!********************************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(out) :: emtd
    real(wp),intent(out) :: grdmtd(3,mol%nat)
    real(wp) :: dum
    emtd = 0.0_wp
    grdmtd = 0.0_wp

    select case (pot%mtdtype)
    case (cv_std_mtd)
      call calc_damp(pot,pot%damptype,0.0_wp)

    case (cv_rmsd)
      dum = float(pot%cvdump)
      call calc_damp(pot,cv_rmsd,dum)
      call calc_rmsd_mtd(mol,pot,emtd,grdmtd)

    case (cv_rmsd_static)
      dum = float(pot%cvdump)
      call calc_damp(pot,cv_rmsd_static,dum)
      call calc_rmsd_mtd(mol,pot,emtd,grdmtd)

    case default
      emtd = 0.0_wp
      grdmtd = 0.0_wp

    end select

    return
  end subroutine calc_mtd

!========================================================================================!
  subroutine calc_damp(pot,dt,x)
!**************************************************
!* subroutine calc_damp
!* calculate a MTD-type-specific damping factor.
!* If/how/where the damping factor is applied
!* depends on the MTD type
!**************************************************
    implicit none
    type(mtdpot) :: pot
    integer :: dt
    real(wp),intent(in) :: x
    real(wp),parameter :: tol = 0.9999_wp
    select case (dt) !>-- select damping parameter calculation
    case (cv_rmsd,cv_rmsd_static)
      pot%damp = (2.0_wp/(1.0_wp+ &
      &       exp(-pot%ramp*x))-1.0_wp)  !> x is pot%cvdump as float

    case (damp_heaviside) !> simple heaviside switch
      pot%damp = sign(0.5_wp,x)+0.5_wp

    case default
      pot%damp = 1.0_wp

    end select

    return
  end subroutine calc_damp

  function check_dump_steps_rmsd(pot) result(steps)
!**************************************************
!* Check the number of MD steps that are affected
!* by the damping
!**************************************************
    implicit none
    type(mtdpot) :: pot
    integer :: steps
    real(wp) :: dum
    real(wp),parameter :: tol = 0.9999_wp
    steps = 0
    do
      steps = steps+1
      dum = float(steps)
      call calc_damp(pot,cv_rmsd,dum)
      if (pot%damp > tol) then
        pot%damp = 0.0_wp
        exit
      end if
    end do
  end function check_dump_steps_rmsd

!========================================================================================!
  subroutine calc_damp2(pot,t,damp)
!*********************************************
!* subroutine calc_damp2
!* damping routine for snapshot-cv-specific
!* damping parameter
!*********************************************
    implicit none
    type(mtdpot) :: pot
    integer :: t !> snapshot
    real(wp) :: damp

    select case (pot%damptype) !>-- select damping parameter calculation
    case (damp_heaviside_cv) !> simple heaviside switch
      damp = sign(0.5_wp,pot%damping(t))+0.5_wp
    case default
      pot%damp = 1.0_wp
    end select

    return
  end subroutine calc_damp2

!========================================================================================!
  subroutine calc_rmsd_mtd(mol,pot,ebias,grdmtd)
!**************************************************************
!* subroutine calc_rmsd_mtd
!* calculate energy and gradient contribution from the RMSD
!* of the current structure (mol) to any structure in a list
!* of documented references.
!* Optionally, atoms for which the RMSD is to be calculated
!* can be specified.
!* Since RMSD calculation can be costly for many structures
!* there is some OMP parallelization going on.
!**************************************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(out) :: ebias
    real(wp),intent(out) :: grdmtd(3,mol%nat)

    real(wp),allocatable :: xyzref(:,:)
    real(wp),allocatable :: xyzcp(:,:)
    real(wp),allocatable :: xyzwhole(:,:)
    real(wp),allocatable :: grad(:,:)
    real(wp) :: U(3,3),x_center(3),y_center(3)
    real(wp) :: rmsdval,E,dEdr

    integer :: i,j,k,l
    logical :: usewhole

    ebias = 0.0_wp
    grdmtd = 0.0_wp

    if (pot%ncur < 1) return
    usewhole = pot%whole .and. allocated(mol%lat) .and. allocated(pot%whole_mask) &
    & .and. allocated(pot%whole_bonds) .and. pot%whole_nbonds > 0
    if (usewhole) then
      allocate (xyzwhole(3,mol%nat),source=0.0_wp)
      call pbc_make_whole(mol%nat,mol%xyz,mol%lat,pot%whole_mask, &
      & pot%whole_bonds,pot%whole_nbonds,xyzwhole)
    end if

    if (.not.allocated(pot%atinclude)) then !>-- include all atoms in RMSD
      allocate (xyzref(3,mol%nat),grad(3,mol%nat),source=0.0_wp)
      !$omp parallel default(none) &
      !$omp shared(pot,mol,xyzwhole,usewhole) &
      !$omp private(grad,xyzref,U,x_center,y_center,rmsdval,E,dEdr) &
      !$omp reduction(+:ebias,grdmtd)
      !$omp do schedule(dynamic)
      do i = 1,pot%ncur
        grad = 0.0_wp
        xyzref = pot%cvxyz(:,:,i)
        if (usewhole) then
          call rmsd(mol%nat,xyzwhole,xyzref,1,U,x_center,y_center,rmsdval, &
          &          .true.,grad)
        else
          call rmsd(mol%nat,mol%xyz,xyzref,1,U,x_center,y_center,rmsdval, &
          &          .true.,grad)
        end if
        E = pot%kpush*exp(-pot%alpha*rmsdval**2)
        if (i == pot%ncur.or.pot%mtdtype == cv_rmsd_static) then
          E = E*pot%damp
        end if
        ebias = ebias+E
        dEdr = -2.0_wp*pot%alpha*e*rmsdval
        grdmtd = grdmtd+dEdr*grad
      end do
      !$omp enddo
      !$omp end parallel
      deallocate (grad,xyzref)

    else !>--- use only selected atoms in RMSD
      k = count(pot%atinclude,1)
      if (k < 1) then
        if (allocated(xyzwhole)) deallocate (xyzwhole)
        return
      end if
      allocate (xyzcp(3,k),xyzref(3,k),grad(3,k),source=0.0_wp)
      !$omp parallel default(none) &
      !$omp shared(pot,mol,k,xyzwhole,usewhole) &
      !$omp private(grad,xyzref,U,x_center,y_center,rmsdval,E,dEdr) &
      !$omp private(xyzcp,j,l) &
      !$omp reduction(+:ebias,grdmtd)
      !$omp do schedule(dynamic)
      do i = 1,pot%ncur
        grad = 0.0_wp
        l = 0
        do j = 1,mol%nat
          if (pot%atinclude(j)) then
            l = l+1
            if (usewhole) then
              xyzcp(:,l) = xyzwhole(:,j)
            else
              xyzcp(:,l) = mol%xyz(:,j)
            end if
            xyzref(:,l) = pot%cvxyz(:,j,i)
          end if
        end do
        call rmsd(k,xyzcp,xyzref,1,U,x_center,y_center,rmsdval, &
        &          .true.,grad)
        E = pot%kpush*exp(-pot%alpha*rmsdval**2)
        if (i == pot%ncur.or.pot%mtdtype == cv_rmsd_static) then
          E = E*pot%damp
        end if
        ebias = ebias+E
        dEdr = -2.0_wp*pot%alpha*e*rmsdval
        l = 0
        do j = 1,mol%nat
          if (pot%atinclude(j)) then
            l = l+1
            grdmtd(:,j) = grdmtd(:,j)+dEdr*grad(:,l)
          end if
        end do
      end do
      !$omp enddo
      !$omp end parallel
      deallocate (grad,xyzref,xyzcp)
    end if

    if (allocated(xyzwhole)) deallocate (xyzwhole)

    return
  end subroutine calc_rmsd_mtd

!========================================================================================!
  subroutine calc_std_mtd(mol,pot,cvt,ebias,grdmtd)
!*******************************************************
!* subroutine calc_std_mtd
!* calculate energy and gradient contribution from the
!* standard MTD formulation (list of CVs)
!*******************************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp) :: cvt  !> value of the CV at the current timestep
    real(wp),intent(out) :: ebias
    real(wp),intent(out) :: grdmtd(3,mol%nat)

    real(wp) :: U(3,3),x_center(3),y_center(3)
    real(wp) :: rmsdval,E,dEdr,dcv,damp2

    integer :: i,j,k,l

    ebias = 0.0_wp
    grdmtd = 0.0_wp

    if (pot%ncur < 1) return

    !$omp parallel default(none) &
    !$omp shared(pot,mol,cvt) &
    !$omp private(U,x_center,y_center,dcv,damp2,E,dEdr) &
    !$omp reduction(+:ebias,grdmtd)
    !$omp do schedule(dynamic)
    do i = 1,pot%ncur
      dcv = cvt-pot%cv(i)
      E = pot%kpush*exp(-pot%alpha*dcv**2)  !> Gaussian shaped potential
      E = E*pot%damp
      call calc_damp2(pot,i,damp2)
      E = E*damp2
      ebias = ebias+E
      dEdr = -2.0_wp*pot%alpha*e*dcv
      grdmtd = grdmtd+dEdr*pot%cvgrd
    end do
    !$omp enddo
    !$omp end parallel

    return

  end subroutine calc_std_mtd

!========================================================================================!
end module metadynamics_module
