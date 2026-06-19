!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023-2026 Philipp Pracht
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
!================================================================================!

!====================================================!
! module gfnff_api
! An interface to GFN-FF standalone calculations
!====================================================!

module gfnff_api
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
  use molecule_parameters,only:aatoau
#ifdef WITH_GFNFF
  use gfnff_interface
#endif
  implicit none
  private

  !> link to subproject routines/datatypes
  public :: gfnff_data
  !> routines from this file
  public :: gfnff_api_setup
  public :: gfnff_sp
  public :: gfnff_printout
  public :: gfnff_getwbos
  public :: gfnff_dump_sasa

#ifndef WITH_GFNFF
  !> these are placeholders if no gfnff module is used!
  type :: gfnff_data
    integer :: id = 0
    character(len=:),allocatable :: parametrisation
    logical :: restart = .false.
    character(len=:),allocatable :: restartfile
    character(len=:),allocatable :: refgeo
  end type gfnff_data
#endif

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gfnff_api_setup(mol,chrg,ff_dat,io,pr,iunit,version)
!*************************************************************
!* Set up (initialize) a GFN-FF calculator from a coord mol. *
!* Lattice vectors are read from mol%lat when present, so    *
!* PBC calculations are automatically activated.             *
!*                                                           *
!* INPUT:                                                    *
!*   mol    - molecule (coords + optional lattice)           *
!*   chrg   - total molecular charge                         *
!*   pr     - optional verbosity flag (logical)              *
!*   iunit  - optional output unit                           *
!* OUTPUT:                                                   *
!*   ff_dat - initialized GFN-FF data object                 *
!*   io     - error status (0 = success)                     *
!*************************************************************
    implicit none
    type(coord),intent(in)      :: mol
    integer,intent(in)          :: chrg
    integer,intent(out)         :: io
    logical,intent(in),optional :: pr
    integer,intent(in),optional :: iunit
    integer,intent(in),optional :: version
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    type(coord) :: refmol
    !> LOCAL
    integer :: mylevel,myunit
    io = 0

    ! ── map legacy pr/iunit to integer printlevel/printunit ──────────────────
    mylevel = 0
    if (present(pr)) then
      if (pr) mylevel = 2
    end if
    if (present(iunit)) then
      myunit = iunit
    else
      myunit = stdout
    end if

#ifdef WITH_GFNFF
    if (allocated(ff_dat%refgeo)) then
      ! ── initialize from a separate reference structure ──────────────────────
      call refmol%open(ff_dat%refgeo)
      if (present(version)) then
        call gfnff_initialize_coord(refmol,chrg,ff_dat,io,mylevel,myunit,version=version)
      else
        call gfnff_initialize_coord(refmol,chrg,ff_dat,io,mylevel,myunit)
      end if
      call refmol%deallocate()
    else
      ! ── initialize from mol directly ────────────────────────────────────────
      if (present(version)) then
        call gfnff_initialize_coord(mol,chrg,ff_dat,io,mylevel,myunit,version=version)
      else
        call gfnff_initialize_coord(mol,chrg,ff_dat,io,mylevel,myunit)
      end if
    end if

#else /* WITH_GFNFF */
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_api_setup

!========================================================================================!

  subroutine gfnff_sp(mol,ff_dat,energy,gradient,iostatus,sigma)
!*************************************************************
!* GFN-FF single-point energy + gradient for mol.           *
!* When mol%lat is allocated the periodic singlepoint call  *
!* is used automatically (lattice passed to gfnff).         *
!*                                                           *
!* INPUT:                                                    *
!*   mol      - molecule (coords + optional lattice)        *
!*   ff_dat   - initialized GFN-FF data object              *
!* OUTPUT:                                                   *
!*   energy   - total energy (Hartree)                      *
!*   gradient - gradient (Eh/Bohr)                          *
!*   iostatus - error status (0 = success)                  *
!*   sigma    - optional stress tensor (Eh); zero non-PBC   *
!*************************************************************
    implicit none
    !> INPUT
    type(coord),intent(in) :: mol
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,mol%nat)
    integer,intent(out)  :: iostatus
    real(wp),intent(out),optional :: sigma(3,3)
    !> LOCAL
    real(wp),allocatable :: xyz_work(:,:)
    real(wp) :: sigma_loc(3,3)
    energy = 0.0_wp
    gradient = 0.0_wp
    iostatus = 0
    sigma_loc = 0.0_wp
#ifdef WITH_GFNFF
    if (allocated(mol%lat)) then
      call gfnff_pbc_make_whole(mol%nat,mol%at,mol%xyz,mol%lat,xyz_work)
      call gfnff_singlepoint(mol%nat,mol%at,xyz_work,ff_dat, &
      &   energy,gradient,lattice=mol%lat,sigma=sigma_loc,iostat=iostatus)
    else
      call gfnff_singlepoint(mol%nat,mol%at,mol%xyz,ff_dat, &
      &   energy,gradient,iostat=iostatus)
    end if
    if (present(sigma)) sigma = sigma_loc
#else
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_sp

!========================================================================================!
  subroutine gfnff_initialize_coord(mol,chrg,ff_dat,io,mylevel,myunit,version)
!*************************************************************
!* Common GFN-FF initialization path. Periodic inputs get a  *
!* temporary make-whole coordinate copy before topology setup *
!* so split molecules are not perceived as separate fragments.*
!*************************************************************
    implicit none
    type(coord),intent(in)      :: mol
    integer,intent(in)          :: chrg,mylevel,myunit
    integer,intent(out)         :: io
    integer,intent(in),optional :: version
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    real(wp),allocatable :: xyz_work(:,:)
    integer :: nbonds_whole

#ifdef WITH_GFNFF
    if (allocated(mol%lat)) then
      call gfnff_pbc_make_whole(mol%nat,mol%at,mol%xyz,mol%lat,xyz_work,nbonds_whole)
      if (mylevel >= 2) write (myunit,'(10x,"CREST PBC make-whole bonds:",1x,i0)') nbonds_whole
      if (present(version)) then
        call gfnff_initialize(mol%nat,mol%at,xyz_work,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io, &
        &   version=version,lattice=mol%lat,npbc=3)
      else
        call gfnff_initialize(mol%nat,mol%at,xyz_work,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io, &
        &   lattice=mol%lat,npbc=3)
      end if
    else
      if (present(version)) then
        call gfnff_initialize(mol%nat,mol%at,mol%xyz,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io, &
        &   version=version)
      else
        call gfnff_initialize(mol%nat,mol%at,mol%xyz,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io)
      end if
    end if
#else
    io = 1
#endif
  end subroutine gfnff_initialize_coord

!========================================================================================!
  subroutine gfnff_pbc_make_whole(nat,at,xyz,lat,whole,nbonds_found)
!*************************************************************
!* Build a simple MIC covalent graph and reconstruct each     *
!* connected component into a whole image. This affects only  *
!* the temporary coordinate array sent to GFN-FF.             *
!*************************************************************
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat),lat(3,3)
    real(wp),allocatable,intent(out) :: whole(:,:)
    integer,intent(out),optional :: nbonds_found
    integer,allocatable :: bonds(:,:)
    integer :: nbonds

    allocate (whole(3,nat),source=xyz)
    if (present(nbonds_found)) nbonds_found = 0
    if (nat < 2) return

    call gfnff_build_whole_bonds(nat,at,xyz,lat,bonds,nbonds)
    if (present(nbonds_found)) nbonds_found = nbonds
    if (nbonds > 0) call gfnff_apply_whole_bonds(nat,xyz,lat,bonds,nbonds,whole)
    if (allocated(bonds)) deallocate (bonds)
  end subroutine gfnff_pbc_make_whole

!========================================================================================!
  subroutine gfnff_build_whole_bonds(nat,at,xyz,lat,bonds,nbonds)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat),lat(3,3)
    integer,allocatable,intent(out) :: bonds(:,:)
    integer,intent(out) :: nbonds

    integer :: i,j,nmax
    integer,allocatable :: tmp(:,:)
    real(wp) :: invlat(3,3),dx(3),mic(3),rcut

    nbonds = 0
    nmax = nat*(nat-1)/2
    if (nmax < 1) return
    allocate (tmp(2,nmax),source=0)
    call gfnff_invert_lat3(lat,invlat)

    do i = 1,nat-1
      do j = i+1,nat
        dx(:) = xyz(:,j)-xyz(:,i)
        call gfnff_minimum_image(dx,lat,invlat,mic)
        rcut = 1.25_wp*(gfnff_whole_covrad(at(i))+gfnff_whole_covrad(at(j)))
        if (norm2(mic) <= rcut) then
          nbonds = nbonds+1
          tmp(1,nbonds) = i
          tmp(2,nbonds) = j
        end if
      end do
    end do

    if (nbonds > 0) allocate (bonds(2,nbonds),source=tmp(:,1:nbonds))
    deallocate (tmp)
  end subroutine gfnff_build_whole_bonds

!========================================================================================!
  subroutine gfnff_apply_whole_bonds(nat,xyz,lat,bonds,nbonds,whole)
    implicit none
    integer,intent(in) :: nat,nbonds
    real(wp),intent(in) :: xyz(3,nat),lat(3,3)
    integer,intent(in) :: bonds(2,nbonds)
    real(wp),intent(inout) :: whole(3,nat)

    logical,allocatable :: seen(:)
    integer,allocatable :: queue(:)
    integer :: seed,head,tail,i,j,b
    real(wp) :: invlat(3,3),dx(3),mic(3)

    if (nbonds < 1) return
    call gfnff_invert_lat3(lat,invlat)
    allocate (seen(nat),source=.false.)
    allocate (queue(nat),source=0)

    do seed = 1,nat
      if (seen(seed)) cycle
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
          call gfnff_minimum_image(dx,lat,invlat,mic)
          whole(:,j) = whole(:,i)+mic(:)
          seen(j) = .true.
          tail = tail+1
          queue(tail) = j
        end do
      end do
    end do

    deallocate (queue,seen)
  end subroutine gfnff_apply_whole_bonds

!========================================================================================!
  subroutine gfnff_minimum_image(dx,lat,invlat,mic)
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
  end subroutine gfnff_minimum_image

!========================================================================================!
  subroutine gfnff_invert_lat3(lat,invlat)
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
  end subroutine gfnff_invert_lat3

!========================================================================================!
  pure real(wp) function gfnff_whole_covrad(at)
    implicit none
    integer,intent(in) :: at
    real(wp) :: rad
    select case (at)
    case (1);  rad = 0.31_wp
    case (3);  rad = 1.28_wp
    case (4);  rad = 0.96_wp
    case (5);  rad = 0.84_wp
    case (6);  rad = 0.76_wp
    case (7);  rad = 0.71_wp
    case (8);  rad = 0.66_wp
    case (9);  rad = 0.57_wp
    case (11); rad = 1.66_wp
    case (12); rad = 1.41_wp
    case (13); rad = 1.21_wp
    case (14); rad = 1.11_wp
    case (15); rad = 1.07_wp
    case (16); rad = 1.05_wp
    case (17); rad = 1.02_wp
    case (19); rad = 2.03_wp
    case (20); rad = 1.76_wp
    case (22); rad = 1.60_wp
    case (23); rad = 1.53_wp
    case (24); rad = 1.39_wp
    case (25); rad = 1.39_wp
    case (26); rad = 1.32_wp
    case (27); rad = 1.26_wp
    case (28); rad = 1.24_wp
    case (29); rad = 1.32_wp
    case (30); rad = 1.22_wp
    case (35); rad = 1.20_wp
    case (40); rad = 1.75_wp
    case (53); rad = 1.39_wp
    case (72); rad = 1.75_wp
    case default
      rad = 1.00_wp
    end select
    gfnff_whole_covrad = rad*aatoau
  end function gfnff_whole_covrad

!========================================================================================!

  subroutine gfnff_printout(iunit,ff_dat)
    implicit none
    !> INPUT
    integer,intent(in)  :: iunit
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
#ifdef WITH_GFNFF
    call ff_dat%resultprint(printunit=iunit)
#else
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_printout

!========================================================================================!
  subroutine gfnff_getwbos(ff_dat,nat,wbo)
!********************************************************
!* obtain connectivity information from GFN-FF topology
!* This is obviously not a true WBO
!********************************************************
    implicit none
    type(gfnff_data),intent(in) :: ff_dat
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)

    wbo = 0.0_wp
#ifdef WITH_GFNFF
    call gfnff_get_fake_wbo(ff_dat,nat,wbo)
#endif
  end subroutine gfnff_getwbos

!========================================================================================!

  subroutine gfnff_dump_sasa(ff_dat,nat,atlist)
!***********************************************************
!* Dumps cumulative SASA for all the atoms .true. in atlist
!* into an file called fort.5454
!***********************************************************
    implicit none
    type(gfnff_data),intent(in) :: ff_dat
    logical,intent(in) :: atlist(nat)
    integer,intent(in) :: nat
    integer :: i
    real(wp) :: sumsasa
    if (allocated(ff_dat%solvation)) then
      if (allocated(ff_dat%solvation%sasa)) then
        sumsasa = 0.0_wp
        do i = 1,nat
          if (atlist(i)) sumsasa = sumsasa+ff_dat%solvation%sasa(i)
        end do
        write (5454,*) sumsasa
      end if
    end if
  end subroutine gfnff_dump_sasa
!========================================================================================!
!========================================================================================!
end module gfnff_api
