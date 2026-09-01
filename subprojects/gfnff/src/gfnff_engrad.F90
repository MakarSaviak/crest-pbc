!================================================================================!
! This file is part of gfnff.
!
! Copyright (C) 2023 Philipp Pracht
!
! gfnff is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! gfnff is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with gfnff. If not, see <https://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module gfnff_engrad_module
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use gfnff_ini2
  use gfnff_data_types,only:TGFFData,TGFFNeighbourList,new,TGFFTopology
  use gfnff_gbsa,only:TBorn
  use gfnff_gdisp0,only:d3_workspace
  use gfnff_param,only:sqrtZr4r2
  use gfnff_helpers
  use gfnff_math_wrapper
  implicit none
  private
  public :: gfnff_eg,gfnff_results,gfnff_workspace

  type :: gfnff_results
    real(wp) :: e_total = 0.0_wp
    real(wp) :: e_rep = 0.0_wp
    real(wp) :: e_es = 0.0_wp
    real(wp) :: e_disp = 0.0_wp
    real(wp) :: e_xb = 0.0_wp
    real(wp) :: g_born = 0.0_wp
    real(wp) :: g_sasa = 0.0_wp
    real(wp) :: g_hb = 0.0_wp
    real(wp) :: g_shift = 0.0_wp
    real(wp) :: dipole(3) = (/0.0_wp,0.0_wp,0.0_wp/)
    real(wp) :: g_solv = 0.0_wp
    real(wp) :: gnorm = 0.0_wp
    real(wp) :: e_bond = 0.0_wp
    real(wp) :: e_angl = 0.0_wp
    real(wp) :: e_tors = 0.0_wp
    real(wp) :: e_hb = 0.0_wp
    real(wp) :: e_batm = 0.0_wp
    real(wp) :: e_ext = 0.0_wp
  end type gfnff_results

  !> Persistent storage and exact caches for repeated frozen-host calls.
  !> No charge, CN, D3 coefficient, cutoff, or force-field approximation is used.
  type :: gfnff_workspace
    integer :: nat = 0
    integer :: nfrag = 0
    integer :: nbond = 0
    integer :: nangl = 0
    integer :: ntors = 0
    real(wp),allocatable :: grab0(:,:,:),rab0(:),eeqtmp(:,:)
    real(wp),allocatable :: cn(:),dcn(:,:,:),qtmp(:)
    real(wp),allocatable :: hb_cn(:),hb_dcn(:,:,:)
    real(wp),allocatable :: sqrab(:),srab(:),g5tmp(:,:),serial_grad(:,:)
    integer,allocatable :: d3list(:,:),d3count(:),d3offset(:),active_atoms(:)
    logical,allocatable :: frozen_pair(:),frozen_mask_ref(:)
    real(wp),allocatable :: frozen_sqrab(:),frozen_srab(:),frozen_reference(:,:)
    real(wp),allocatable :: raw_cn_static(:)
    type(d3_workspace) :: d3_scratch
    real(wp),allocatable :: eeq_a(:,:),eeq_x(:),eeq_lapack_work(:)
    integer,allocatable :: eeq_ipiv(:)
    real(wp),allocatable :: eeq_frozen_block(:,:)
    real(wp),allocatable :: eeq_host_inverse(:,:),eeq_block_k(:,:),eeq_block_y(:,:)
    real(wp),allocatable :: eeq_block_s(:,:),eeq_host_rhs(:),eeq_host_solution(:)
    real(wp),allocatable :: eeq_block_rhs(:),eeq_block_z(:)
    real(wp),allocatable :: eeq_host_fact_work(:),eeq_block_fact_work(:),eeq_host_inv_work(:)
    integer,allocatable :: eeq_host_ipiv(:),eeq_block_ipiv(:)
    integer,allocatable :: active_pair_i(:),active_pair_j(:),active_pair_idx(:)
    integer,allocatable :: dynamic_bond_idx(:),dynamic_angle_idx(:),dynamic_torsion_idx(:)
    integer :: n_active_pairs = 0
    integer :: frozen_prefix = 0
    logical :: prefix_frozen = .false.
    integer :: n_dynamic_bonds = 0
    integer :: n_dynamic_angles = 0
    integer :: n_dynamic_torsions = 0
    logical :: frozen_cache_valid = .false.
    logical :: static_cache_valid = .false.
    logical :: cn_cache_valid = .false.
    logical :: eeq_frozen_block_valid = .false.
    logical :: eeq_host_inverse_valid = .false.
    real(wp) :: static_repthr = -1.0_wp
    real(wp) :: static_dispthr = -1.0_wp
    real(wp) :: cn_cache_thr = -1.0_wp
    real(wp) :: cached_nb_rep_total = 0.0_wp
    real(wp) :: cached_bond_rep_total = 0.0_wp
    real(wp) :: cached_angle_total = 0.0_wp
    real(wp) :: cached_torsion_total = 0.0_wp
  contains
    procedure :: ensure => gfnff_workspace_ensure
    procedure :: prepare_frozen => gfnff_workspace_prepare_frozen
    procedure :: prepare_static => gfnff_workspace_prepare_static
    procedure :: prepare_eeq_frozen_block => gfnff_workspace_prepare_eeq_frozen_block
    procedure :: prepare_eeq_host_inverse => gfnff_workspace_prepare_eeq_host_inverse
    procedure :: release => gfnff_workspace_release
  end type gfnff_workspace

  real(wp),private,parameter :: pi = 3.1415926535897932385_wp
  real(wp),private,parameter :: sqrtpi = 1.77245385091_wp

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gfnff_workspace_ensure(self,n,nfrag,nbond,nangl,ntors)
    class(gfnff_workspace),intent(inout) :: self
    integer,intent(in) :: n,nfrag,nbond,nangl,ntors
    integer :: npair,m

    if (self%nat == n .and. self%nfrag == nfrag .and. self%nbond == nbond .and. self%nangl == nangl .and. &
   &    self%ntors == ntors .and. allocated(self%sqrab)) return

    call self%release()
    npair = n*(n+1)/2
    m = n+nfrag
    allocate(self%sqrab(npair),self%srab(npair),self%qtmp(n),self%g5tmp(3,n), &
   &         self%serial_grad(3,n), &
   &         self%eeqtmp(2,npair),self%d3list(2,npair),self%dcn(3,n,n), &
   &         self%cn(n),self%hb_dcn(3,n,n),self%hb_cn(n), &
   &         self%d3count(n),self%d3offset(n),self%active_atoms(n), &
   &         self%grab0(3,n,nbond),self%rab0(nbond), &
   &         self%frozen_pair(npair),self%frozen_mask_ref(n), &
   &         self%frozen_sqrab(npair),self%frozen_srab(npair), &
   &         self%frozen_reference(3,n),self%raw_cn_static(n), &
   &         self%active_pair_i(npair), &
   &         self%active_pair_j(npair),self%active_pair_idx(npair), &
   &         self%dynamic_bond_idx(max(1,nbond)), &
   &         self%dynamic_angle_idx(max(1,nangl)), &
   &         self%dynamic_torsion_idx(max(1,ntors)), &
   &         self%eeq_a(m,m),self%eeq_x(m), &
   &         self%eeq_ipiv(m),self%eeq_frozen_block(n,n))
    self%frozen_pair = .false.
    self%frozen_mask_ref = .false.
    self%frozen_cache_valid = .false.
    self%static_cache_valid = .false.
    self%cn_cache_valid = .false.
    self%eeq_frozen_block_valid = .false.
    self%eeq_host_inverse_valid = .false.
    self%n_active_pairs = 0
    self%frozen_prefix = 0
    self%prefix_frozen = .false.
    self%n_dynamic_bonds = 0
    self%n_dynamic_angles = 0
    self%n_dynamic_torsions = 0
    self%cached_nb_rep_total = 0.0_wp
    self%cached_bond_rep_total = 0.0_wp
    self%cached_angle_total = 0.0_wp
    self%cached_torsion_total = 0.0_wp
    self%nat = n
    self%nfrag = nfrag
    self%nbond = nbond
    self%nangl = nangl
    self%ntors = ntors
  end subroutine gfnff_workspace_ensure

  subroutine gfnff_workspace_prepare_frozen(self,n,xyz,frozen_mask)
    class(gfnff_workspace),intent(inout) :: self
    integer,intent(in) :: n
    real(wp),intent(in) :: xyz(3,n)
    logical,intent(in) :: frozen_mask(n)
    integer :: i,j,ij,k,m
    real(wp) :: frozen_delta,cache_scale,cache_tolerance

    if (count(frozen_mask) < 2) then
      self%frozen_pair = .false.
      self%frozen_cache_valid = .false.
      self%static_cache_valid = .false.
      self%cn_cache_valid = .false.
      self%eeq_frozen_block_valid = .false.
      self%eeq_host_inverse_valid = .false.
      self%n_active_pairs = 0
      self%frozen_prefix = 0
      self%prefix_frozen = .false.
      self%n_dynamic_bonds = 0
      self%n_dynamic_angles = 0
      self%n_dynamic_torsions = 0
      return
    end if

    if (self%frozen_cache_valid) then
      if (any(self%frozen_mask_ref .neqv. frozen_mask)) then
        self%frozen_cache_valid = .false.
        self%static_cache_valid = .false.
        self%cn_cache_valid = .false.
        self%eeq_frozen_block_valid = .false.
        self%eeq_host_inverse_valid = .false.
      else
        ! Preserve the tolerant frozen-coordinate validation while avoiding
        ! three full-size SPREAD temporaries on every energy/gradient call.
        frozen_delta = 0.0_wp
        cache_scale = 1.0_wp
        do i = 1,n
          if (.not.frozen_mask(i)) cycle
          do k = 1,3
            frozen_delta = max(frozen_delta, &
           &                   abs(xyz(k,i)-self%frozen_reference(k,i)))
            cache_scale = max(cache_scale,abs(xyz(k,i)), &
           &                  abs(self%frozen_reference(k,i)))
          end do
        end do
        cache_tolerance = 128.0_wp*epsilon(1.0_wp)*cache_scale
        if (frozen_delta > cache_tolerance) then
          self%frozen_cache_valid = .false.
          self%static_cache_valid = .false.
          self%cn_cache_valid = .false.
          self%eeq_frozen_block_valid = .false.
          self%eeq_host_inverse_valid = .false.
        end if
      end if
    end if
    if (self%frozen_cache_valid) return

    self%frozen_pair = .false.
    self%frozen_mask_ref = frozen_mask
    self%n_active_pairs = 0
    self%frozen_prefix = 0
    do i = 1,n
      if (.not.frozen_mask(i)) exit
      self%frozen_prefix = i
    end do
    if (self%frozen_prefix == n) then
      self%prefix_frozen = .true.
    else
      self%prefix_frozen = self%frozen_prefix >= 2 .and. &
     &                     all(.not.frozen_mask(self%frozen_prefix+1:n))
    end if
    do i = 1,n
      if (frozen_mask(i)) self%frozen_reference(:,i) = xyz(:,i)
      ij = i*(i-1)/2
      self%frozen_sqrab(ij+i) = 0.0_wp
      self%frozen_srab(ij+i) = 0.0_wp
      self%sqrab(ij+i) = 0.0_wp
      self%srab(ij+i) = 0.0_wp
      do j = 1,i-1
        k = ij+j
        if (frozen_mask(i).and.frozen_mask(j)) then
          self%frozen_pair(k) = .true.
          self%frozen_sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + &
         &                       (xyz(2,i)-xyz(2,j))**2 + &
         &                       (xyz(3,i)-xyz(3,j))**2
          self%frozen_srab(k) = sqrt(self%frozen_sqrab(k))
          self%sqrab(k) = self%frozen_sqrab(k)
          self%srab(k) = self%frozen_srab(k)
        else
          self%n_active_pairs = self%n_active_pairs+1
          m = self%n_active_pairs
          self%active_pair_i(m) = i
          self%active_pair_j(m) = j
          self%active_pair_idx(m) = k
        end if
      end do
    end do
    self%frozen_cache_valid = .true.
    self%static_cache_valid = .false.
    self%cn_cache_valid = .false.
    self%eeq_frozen_block_valid = .false.
    self%eeq_host_inverse_valid = .false.
  end subroutine gfnff_workspace_prepare_frozen

  subroutine gfnff_workspace_prepare_eeq_frozen_block(self,n,r,topo)
    class(gfnff_workspace),intent(inout) :: self
    integer,intent(in) :: n
    real(wp),intent(in) :: r(n*(n+1)/2)
    type(TGFFTopology),intent(in) :: topo
    integer :: i,j,k,ij,nf
    real(wp) :: gammij,tmp
    real(wp),parameter :: tsqrt2pi = 0.797884560802866_wp

    if (.not.self%frozen_cache_valid .or. .not.self%prefix_frozen) then
      self%eeq_frozen_block_valid = .false.
      self%eeq_host_inverse_valid = .false.
      return
    end if
    if (self%eeq_frozen_block_valid) return

    nf = self%frozen_prefix
    self%eeq_frozen_block(1:nf,1:nf) = 0.0_wp
    do i = 1,nf
      self%eeq_frozen_block(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i)
      k = i*(i-1)/2
      do j = 1,i-1
        ij = k+j
        gammij = 1.0_wp/sqrt(topo%alpeeq(i)+topo%alpeeq(j))
        tmp = erf(gammij*r(ij))
        self%eeqtmp(1,ij) = gammij
        self%eeqtmp(2,ij) = tmp
        self%eeq_frozen_block(j,i) = tmp/r(ij)
        self%eeq_frozen_block(i,j) = self%eeq_frozen_block(j,i)
      end do
    end do
    self%eeq_frozen_block_valid = .true.
  end subroutine gfnff_workspace_prepare_eeq_frozen_block

  subroutine gfnff_workspace_release_eeq_host_solver(self)
    class(gfnff_workspace),intent(inout) :: self

    if (allocated(self%eeq_host_inverse)) deallocate(self%eeq_host_inverse)
    if (allocated(self%eeq_block_k)) deallocate(self%eeq_block_k)
    if (allocated(self%eeq_block_y)) deallocate(self%eeq_block_y)
    if (allocated(self%eeq_block_s)) deallocate(self%eeq_block_s)
    if (allocated(self%eeq_host_rhs)) deallocate(self%eeq_host_rhs)
    if (allocated(self%eeq_host_solution)) deallocate(self%eeq_host_solution)
    if (allocated(self%eeq_block_rhs)) deallocate(self%eeq_block_rhs)
    if (allocated(self%eeq_block_z)) deallocate(self%eeq_block_z)
    if (allocated(self%eeq_host_fact_work)) deallocate(self%eeq_host_fact_work)
    if (allocated(self%eeq_block_fact_work)) deallocate(self%eeq_block_fact_work)
    if (allocated(self%eeq_host_inv_work)) deallocate(self%eeq_host_inv_work)
    if (allocated(self%eeq_host_ipiv)) deallocate(self%eeq_host_ipiv)
    if (allocated(self%eeq_block_ipiv)) deallocate(self%eeq_block_ipiv)
    self%eeq_host_inverse_valid = .false.
  end subroutine gfnff_workspace_release_eeq_host_solver

  subroutine gfnff_workspace_prepare_eeq_host_inverse(self,n,topo,io)
    class(gfnff_workspace),intent(inout) :: self
    integer,intent(in) :: n
    type(TGFFTopology),intent(in) :: topo
    integer,intent(out) :: io
    integer :: nf,na,p,i,j,io1,io2,stat_alloc
    logical :: layout_ok

    io = 0
    if (.not.self%eeq_frozen_block_valid .or. .not.self%prefix_frozen) return
    nf = self%frozen_prefix
    na = n-nf
    p = na+topo%nfrag
    if (nf < 1 .or. na < 1 .or. p >= nf) return
    if (self%eeq_host_inverse_valid) return

    layout_ok = allocated(self%eeq_host_inverse) .and. &
   &            allocated(self%eeq_block_k) .and. allocated(self%eeq_block_y) .and. &
   &            allocated(self%eeq_block_s) .and. allocated(self%eeq_host_rhs) .and. &
   &            allocated(self%eeq_host_solution) .and. &
   &            allocated(self%eeq_block_rhs) .and. allocated(self%eeq_block_z) .and. &
   &            allocated(self%eeq_host_ipiv) .and. allocated(self%eeq_block_ipiv) .and. &
   &            allocated(self%eeq_host_inv_work)
    if (layout_ok) then
      layout_ok = size(self%eeq_host_inverse,1) == nf .and. &
     &            size(self%eeq_host_inverse,2) == nf .and. &
     &            size(self%eeq_block_s,1) == p .and. size(self%eeq_block_s,2) == p
    end if
    if (.not.layout_ok) call gfnff_workspace_release_eeq_host_solver(self)
    if (.not.allocated(self%eeq_host_inverse)) then
      allocate(self%eeq_host_inverse(nf,nf),self%eeq_block_k(nf,p), &
     &         self%eeq_block_y(nf,p),self%eeq_block_s(p,p), &
     &         self%eeq_host_rhs(nf),self%eeq_host_solution(nf), &
     &         self%eeq_block_rhs(p),self%eeq_block_z(p), &
     &         self%eeq_host_ipiv(nf),self%eeq_block_ipiv(p), &
     &         self%eeq_host_inv_work(nf),stat=stat_alloc)
      if (stat_alloc /= 0) then
        call gfnff_workspace_release_eeq_host_solver(self)
        io = -1000
        return
      end if
    end if

    self%eeq_host_inverse = self%eeq_frozen_block(1:nf,1:nf)
    call sytrf_cached_wrap(self%eeq_host_inverse,self%eeq_host_ipiv, &
   &                       self%eeq_host_fact_work,io1)
    if (io1 /= 0) then
      io = io1
      return
    end if
    call lapack_sytri('U',nf,self%eeq_host_inverse,nf,self%eeq_host_ipiv, &
   &                  self%eeq_host_inv_work,io2)
    if (io2 /= 0) then
      io = io2
      return
    end if
    do i = 1,nf
      do j = 1,i-1
        self%eeq_host_inverse(i,j) = self%eeq_host_inverse(j,i)
      end do
    end do
    self%eeq_host_inverse_valid = .true.
  end subroutine gfnff_workspace_prepare_eeq_host_inverse

  subroutine gfnff_workspace_prepare_static(self,n,at,xyz,frozen_mask,repthr,dispthr,param,topo,sqrab,srab)
    class(gfnff_workspace),intent(inout) :: self
    integer,intent(in) :: n,at(n)
    real(wp),intent(in) :: xyz(3,n),repthr,dispthr,sqrab(:),srab(:)
    logical,intent(in) :: frozen_mask(n)
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer :: i,j,k,l,m,ij,iat,jat,ati,atj
    real(wp) :: r2,rab,t16,t19,t8,t26,alpha,repab,etmp
    real(wp) :: g3tmp(3,3),g4tmp(3,4)

    if (.not.self%frozen_cache_valid) then
      self%static_cache_valid = .false.
      return
    end if
    if (self%static_cache_valid.and.self%static_repthr == repthr .and. &
   &    self%static_dispthr == dispthr) return

    self%cached_nb_rep_total = 0.0_wp
    self%cached_bond_rep_total = 0.0_wp
    self%cached_angle_total = 0.0_wp
    self%cached_torsion_total = 0.0_wp
    self%n_dynamic_bonds = 0
    self%n_dynamic_angles = 0
    self%n_dynamic_torsions = 0

    do iat = 1,n
      m = iat*(iat-1)/2
      do jat = 1,iat-1
        if (.not.(frozen_mask(iat).and.frozen_mask(jat))) cycle
        ij = m+jat
        r2 = sqrab(ij)
        if (r2 .gt. repthr) cycle
        if (topo%bpair(ij) .eq. 1) cycle
        ati = at(iat)
        atj = at(jat)
        rab = srab(ij)
        t16 = r2**0.75_wp
        t19 = t16*t16
        t8 = t16*topo%alphanb(ij)
        t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
        self%cached_nb_rep_total = self%cached_nb_rep_total+t26/rab
      end do
    end do

    do i = 1,topo%nbond
      iat = topo%blist(1,i)
      jat = topo%blist(2,i)
      if (.not.(frozen_mask(iat).and.frozen_mask(jat))) then
        self%n_dynamic_bonds = self%n_dynamic_bonds+1
        self%dynamic_bond_idx(self%n_dynamic_bonds) = i
        cycle
      end if
      ij = iat*(iat-1)/2+jat
      r2 = sqrab(ij)
      rab = srab(ij)
      ati = at(iat)
      atj = at(jat)
      alpha = sqrt(param%repa(ati)*param%repa(atj))
      repab = param%repz(ati)*param%repz(atj)*param%repscalb
      t16 = r2**0.75_wp
      t19 = t16*t16
      t26 = exp(-alpha*t16)*repab
      self%cached_bond_rep_total = self%cached_bond_rep_total+t26/rab
    end do

    do m = 1,topo%nangl
      j = topo%alist(1,m)
      i = topo%alist(2,m)
      k = topo%alist(3,m)
      if (.not.(frozen_mask(i).and.frozen_mask(j).and.frozen_mask(k))) then
        self%n_dynamic_angles = self%n_dynamic_angles+1
        self%dynamic_angle_idx(self%n_dynamic_angles) = m
        cycle
      end if
      call egbend(m,j,i,k,n,at,xyz,etmp,g3tmp,param,topo)
      self%cached_angle_total = self%cached_angle_total+etmp
    end do

    do m = 1,topo%ntors
      i = topo%tlist(1,m)
      j = topo%tlist(2,m)
      k = topo%tlist(3,m)
      l = topo%tlist(4,m)
      if (.not.(frozen_mask(i).and.frozen_mask(j).and.frozen_mask(k).and.frozen_mask(l))) then
        self%n_dynamic_torsions = self%n_dynamic_torsions+1
        self%dynamic_torsion_idx(self%n_dynamic_torsions) = m
        cycle
      end if
      call egtors(m,i,j,k,l,n,at,xyz,etmp,g4tmp,param,topo)
      self%cached_torsion_total = self%cached_torsion_total+etmp
    end do

    self%static_repthr = repthr
    self%static_dispthr = dispthr
    self%static_cache_valid = .true.
  end subroutine gfnff_workspace_prepare_static

  subroutine gfnff_workspace_release(self)
    class(gfnff_workspace),intent(inout) :: self
    if (allocated(self%grab0)) deallocate(self%grab0)
    if (allocated(self%rab0)) deallocate(self%rab0)
    if (allocated(self%eeqtmp)) deallocate(self%eeqtmp)
    if (allocated(self%cn)) deallocate(self%cn)
    if (allocated(self%dcn)) deallocate(self%dcn)
    if (allocated(self%qtmp)) deallocate(self%qtmp)
    if (allocated(self%hb_cn)) deallocate(self%hb_cn)
    if (allocated(self%hb_dcn)) deallocate(self%hb_dcn)
    if (allocated(self%sqrab)) deallocate(self%sqrab)
    if (allocated(self%srab)) deallocate(self%srab)
    if (allocated(self%g5tmp)) deallocate(self%g5tmp)
    if (allocated(self%serial_grad)) deallocate(self%serial_grad)
    if (allocated(self%d3list)) deallocate(self%d3list)
    if (allocated(self%d3count)) deallocate(self%d3count)
    if (allocated(self%d3offset)) deallocate(self%d3offset)
    if (allocated(self%active_atoms)) deallocate(self%active_atoms)
    if (allocated(self%frozen_pair)) deallocate(self%frozen_pair)
    if (allocated(self%frozen_mask_ref)) deallocate(self%frozen_mask_ref)
    if (allocated(self%frozen_sqrab)) deallocate(self%frozen_sqrab)
    if (allocated(self%frozen_srab)) deallocate(self%frozen_srab)
    if (allocated(self%frozen_reference)) deallocate(self%frozen_reference)
    if (allocated(self%raw_cn_static)) deallocate(self%raw_cn_static)
    call self%d3_scratch%release()
    if (allocated(self%eeq_a)) deallocate(self%eeq_a)
    if (allocated(self%eeq_x)) deallocate(self%eeq_x)
    if (allocated(self%eeq_ipiv)) deallocate(self%eeq_ipiv)
    if (allocated(self%eeq_lapack_work)) deallocate(self%eeq_lapack_work)
    if (allocated(self%eeq_frozen_block)) deallocate(self%eeq_frozen_block)
    call gfnff_workspace_release_eeq_host_solver(self)
    if (allocated(self%active_pair_i)) deallocate(self%active_pair_i)
    if (allocated(self%active_pair_j)) deallocate(self%active_pair_j)
    if (allocated(self%active_pair_idx)) deallocate(self%active_pair_idx)
    if (allocated(self%dynamic_bond_idx)) deallocate(self%dynamic_bond_idx)
    if (allocated(self%dynamic_angle_idx)) deallocate(self%dynamic_angle_idx)
    if (allocated(self%dynamic_torsion_idx)) deallocate(self%dynamic_torsion_idx)
    self%nat = 0
    self%nfrag = 0
    self%nbond = 0
    self%nangl = 0
    self%ntors = 0
    self%frozen_cache_valid = .false.
    self%static_cache_valid = .false.
    self%cn_cache_valid = .false.
    self%eeq_frozen_block_valid = .false.
    self%eeq_host_inverse_valid = .false.
    self%static_repthr = -1.0_wp
    self%static_dispthr = -1.0_wp
    self%cn_cache_thr = -1.0_wp
    self%n_active_pairs = 0
    self%frozen_prefix = 0
    self%prefix_frozen = .false.
    self%n_dynamic_bonds = 0
    self%n_dynamic_angles = 0
    self%n_dynamic_torsions = 0
  end subroutine gfnff_workspace_release

!---------------------------------------------------
!> GFN-FF
!> energy and analytical gradient for given xyz and
!> charge ichrg
!> requires D3 ini (rcov,r2r4,copyc6) as well as
!> gfnff_ini call
!>
!> the total energy is
!> ees + edisp + erep + ebond + eangl + etors + ehb + exb + ebatm + eext
!>
!> uses EEQ charge and D3 routines
!> basic trigonometry for bending and torsion angles
!> taken slightly modified from QMDFF code
!> repulsion and rabguess from xtb GFN0 part
!>
!> requires setup of
!>     integer,allocatable :: blist(:,:)
!>     integer,allocatable :: alist(:,:)
!>     integer,allocatable :: tlist(:,:)
!>     integer,allocatable ::b3list(:,:)
!>     real(wp),allocatable:: vbond(:,:)
!>     real(wp),allocatable:: vangl(:,:)
!>     real(wp),allocatable:: vtors(:,:)
!>     chi,gam,alp,cnf
!>     repa,repz,alphanb
!> 
!---------------------------------------------------
  subroutine gfnff_eg(pr,n,ichrg,at,xyz,makeq,g,etot,res_gff, &
  &          param,topo,nlist,solvation,update,version,accuracy,io,work,frozen_mask)

    use gfnff_param,only:efield,gffVersion,gfnff_thresholds
    use gfnff_gdisp0,only:d3_gradient
    use gfnff_cn
    use gfnff_rab
    implicit none
    character(len=*),parameter :: source = 'gfnff_eg'
    type(gfnff_results),intent(out) :: res_gff

    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(TGFFNeighbourList),intent(inout) :: nlist
    type(gfnff_workspace),intent(inout),target :: work
    logical,intent(in),optional :: frozen_mask(n)

    type(TBorn),allocatable,intent(inout) :: solvation
    logical,intent(in) :: update
    integer,intent(in) :: version
    real(wp),intent(in) :: accuracy
    integer,intent(out) :: io
    integer,intent(in) :: n
    integer,intent(in) :: ichrg
    integer,intent(in) :: at(n)
    real(wp),intent(in) :: xyz(3,n)
    real(wp),intent(inout) :: g(3,n)
    real(wp),intent(inout) :: etot
    logical,intent(in) :: pr
    logical,intent(in) :: makeq

    real(wp) :: edisp,ees,ebond,eangl,etors,erep,ehb,exb,ebatm,eext
    real(wp) :: gsolv,gborn,ghb,gsasa,gshift

    integer  :: i,j,k,l,m,ij,nd3,nd3pos,nlocal,nthreads,tid,nactive,iact,nloop
    integer  :: ati,atj,iat,jat
    integer  :: hbA,hbB
    integer  :: lin
    logical  :: ex,require_update,force_hbond_update,serial_inner
    integer  :: nhb1,nhb2,nxb
    real(wp) ::  r2,rab,qq0,erff,dd,dum1,r3(3),t8,dum,t22,t39
    real(wp) ::  dx,dy,dz,yy,t4,t5,t6,alpha,t20
    real(wp) ::  repab,t16,t19,t26,t27,xa,ya,za,cosa,de,t28
    real(wp) ::  gammij,eesinf,etmp,phi
    real(wp) ::  rn,dr,g3tmp(3,3),g4tmp(3,4)
    real(wp) :: rij,drij(3,n),gactive(3)

    real(wp),pointer :: grab0(:,:,:),rab0(:),eeqtmp(:,:)
    real(wp),pointer :: cn(:),dcn(:,:,:),qtmp(:)
    real(wp),pointer :: hb_cn(:),hb_dcn(:,:,:)
    real(wp),pointer :: sqrab(:),srab(:),g5tmp(:,:),gserial(:,:)
    real(wp),allocatable :: ghb_thread(:,:,:),ehb_thread(:)
    integer,pointer :: d3list(:,:),d3count(:),d3offset(:),active_atoms(:)
    !type(tb_timer) :: timer
    real(wp) :: dispthr,cnthr,repthr,hbthr1,hbthr2,serial_energy

    call gfnff_thresholds(accuracy,dispthr,cnthr,repthr,hbthr1,hbthr2)

    io = 0 !> return status

    g = 0
    exb = 0
    ehb = 0
    erep = 0
    ees = 0
    edisp = 0
    ebond = 0
    eangl = 0
    etors = 0
    ebatm = 0
    eext = 0

    gsolv = 0.0d0
    gsasa = 0.0d0
    gborn = 0.0d0
    ghb = 0.0d0
    gshift = 0.0d0

    call work%ensure(n,topo%nfrag,topo%nbond,topo%nangl,topo%ntors)
    ! A CREST optimization worker sets its task-local OpenMP team size to one.
    ! Bypass the nested GFN-FF fork/reduction machinery only in that case;
    ! multi-threaded dynamics continues through the original OpenMP regions.
    serial_inner = .true.
!$  serial_inner = omp_get_max_threads() == 1
    if (present(frozen_mask)) then
      call work%prepare_frozen(n,xyz,frozen_mask)
      if (work%static_cache_valid) then
        if (work%static_repthr /= repthr .or. work%static_dispthr /= dispthr) then
          work%static_cache_valid = .false.
        end if
      end if
      if (work%cn_cache_valid.and.work%cn_cache_thr /= cnthr) then
        work%cn_cache_valid = .false.
      end if
    else
      work%frozen_cache_valid = .false.
      work%static_cache_valid = .false.
      work%cn_cache_valid = .false.
      work%frozen_pair = .false.
    end if
    sqrab => work%sqrab
    srab => work%srab
    qtmp => work%qtmp
    g5tmp => work%g5tmp
    gserial => work%serial_grad
    eeqtmp => work%eeqtmp
    d3list => work%d3list
    dcn => work%dcn
    cn => work%cn
    hb_dcn => work%hb_dcn
    hb_cn => work%hb_cn
    d3count => work%d3count
    d3offset => work%d3offset
    active_atoms => work%active_atoms
    grab0 => work%grab0
    rab0 => work%rab0

!      if (pr) call timer%new(10 + count([allocated(solvation)]),.false.)

!      if (pr) call timer%measure(1,'distance/D3 list')
    if (work%frozen_cache_valid) then
      if (work%prefix_frozen) then
!> Fast path for the CREST host--guest layout: all frozen atoms form a
!> contiguous prefix. Frozen rows are already present in the persistent packed
!> distance arrays; only rows belonging to active atoms are recomputed.
        if (.not.work%static_cache_valid) then
          d3count(1:work%frozen_prefix) = 0
          if (serial_inner) then
            do i = 1,work%frozen_prefix
              ij = i*(i-1)/2
              nlocal = 0
              do j = 1,i-1
                k = ij+j
                if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
              end do
              d3count(i) = nlocal
            end do
          else
            !$omp parallel do default(none) schedule(static) &
            !$omp shared(sqrab, dispthr, d3count, work) &
            !$omp private(i, j, k, ij, nlocal)
            do i = 1,work%frozen_prefix
              ij = i*(i-1)/2
              nlocal = 0
              do j = 1,i-1
                k = ij+j
                if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
              end do
              d3count(i) = nlocal
            end do
            !$omp end parallel do
          end if
        end if

        if (serial_inner) then
          do i = work%frozen_prefix+1,n
            ij = i*(i-1)/2
            nlocal = 0
            do j = 1,i-1
              k = ij+j
              sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + &
             &           (xyz(2,i)-xyz(2,j))**2 + &
             &           (xyz(3,i)-xyz(3,j))**2
              srab(k) = sqrt(sqrab(k))
              if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
            end do
            d3count(i) = nlocal
          end do
        else
          !$omp parallel do default(none) schedule(static) &
          !$omp shared(n, xyz, sqrab, srab, dispthr, d3count, work) &
          !$omp private(i, j, k, ij, nlocal)
          do i = work%frozen_prefix+1,n
            ij = i*(i-1)/2
            nlocal = 0
            do j = 1,i-1
              k = ij+j
              sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + &
             &           (xyz(2,i)-xyz(2,j))**2 + &
             &           (xyz(3,i)-xyz(3,j))**2
              srab(k) = sqrt(sqrab(k))
              if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
            end do
            d3count(i) = nlocal
          end do
          !$omp end parallel do
        end if
      else
!> General-mask fallback: update the packed list of every pair containing at
!> least one active atom, then rebuild threshold counts in original row order.
        d3count = 0
        if (serial_inner) then
          do m = 1,work%n_active_pairs
            i = work%active_pair_i(m)
            j = work%active_pair_j(m)
            k = work%active_pair_idx(m)
            sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + &
           &           (xyz(2,i)-xyz(2,j))**2 + &
           &           (xyz(3,i)-xyz(3,j))**2
            srab(k) = sqrt(sqrab(k))
          end do
        else
          !$omp parallel do default(none) schedule(static) &
          !$omp shared(xyz, sqrab, srab, work) &
          !$omp private(m, i, j, k)
          do m = 1,work%n_active_pairs
            i = work%active_pair_i(m)
            j = work%active_pair_j(m)
            k = work%active_pair_idx(m)
            sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + &
           &           (xyz(2,i)-xyz(2,j))**2 + &
           &           (xyz(3,i)-xyz(3,j))**2
            srab(k) = sqrt(sqrab(k))
          end do
          !$omp end parallel do
        end if

        if (serial_inner) then
          do i = 1,n
            ij = i*(i-1)/2
            nlocal = 0
            do j = 1,i-1
              k = ij+j
              if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
            end do
            d3count(i) = nlocal
          end do
        else
          !$omp parallel do default(none) schedule(static) &
          !$omp shared(n, sqrab, dispthr, d3count) &
          !$omp private(i, j, k, ij, nlocal)
          do i = 1,n
            ij = i*(i-1)/2
            nlocal = 0
            do j = 1,i-1
              k = ij+j
              if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
            end do
            d3count(i) = nlocal
          end do
          !$omp end parallel do
        end if
      end if
    else
      d3count = 0
      if (serial_inner) then
        do i = 1,n
          ij = i*(i-1)/2
          nlocal = 0
          do j = 1,i-1
            k = ij+j
            sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + &
           &           (xyz(2,i)-xyz(2,j))**2 + &
           &           (xyz(3,i)-xyz(3,j))**2
            srab(k) = sqrt(sqrab(k))
            if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
          end do
          d3count(i) = nlocal
          sqrab(ij+i) = 0.0_wp
          srab(ij+i) = 0.0_wp
        end do
      else
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(n, xyz, sqrab, srab, dispthr, d3count) &
        !$omp private(i, j, k, ij, nlocal)
        do i = 1,n
          ij = i*(i-1)/2
          nlocal = 0
          do j = 1,i-1
            k = ij+j
            sqrab(k) = (xyz(1,i)-xyz(1,j))**2 + &
           &           (xyz(2,i)-xyz(2,j))**2 + &
           &           (xyz(3,i)-xyz(3,j))**2
            srab(k) = sqrt(sqrab(k))
            if (sqrab(k) .lt. dispthr) nlocal = nlocal+1
          end do
          d3count(i) = nlocal
          sqrab(ij+i) = 0.0_wp
          srab(ij+i) = 0.0_wp
        end do
        !$omp end parallel do
      end if
    end if

!> Prefix offsets preserve the original deterministic (i,j) pair order.
    d3offset(1) = 0
    do i = 2,n
      d3offset(i) = d3offset(i-1)+d3count(i-1)
    end do
    nd3 = d3offset(n)+d3count(n)

    if (work%prefix_frozen.and.work%static_cache_valid) then
!> Frozen-prefix D3-list entries remain valid and already occupy the beginning
!> of d3list. Rebuild only rows containing an active atom.
      if (serial_inner) then
        do i = work%frozen_prefix+1,n
          ij = i*(i-1)/2
          nd3pos = d3offset(i)
          do j = 1,i-1
            k = ij+j
            if (sqrab(k) .lt. dispthr) then
              nd3pos = nd3pos+1
              d3list(1,nd3pos) = i
              d3list(2,nd3pos) = j
            end if
          end do
        end do
      else
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(n, sqrab, dispthr, d3offset, d3list, work) &
        !$omp private(i, j, k, ij, nd3pos)
        do i = work%frozen_prefix+1,n
          ij = i*(i-1)/2
          nd3pos = d3offset(i)
          do j = 1,i-1
            k = ij+j
            if (sqrab(k) .lt. dispthr) then
              nd3pos = nd3pos+1
              d3list(1,nd3pos) = i
              d3list(2,nd3pos) = j
            end if
          end do
        end do
        !$omp end parallel do
      end if
    else
      if (serial_inner) then
        do i = 1,n
          ij = i*(i-1)/2
          nd3pos = d3offset(i)
          do j = 1,i-1
            k = ij+j
            if (sqrab(k) .lt. dispthr) then
              nd3pos = nd3pos+1
              d3list(1,nd3pos) = i
              d3list(2,nd3pos) = j
            end if
          end do
        end do
      else
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(n, sqrab, dispthr, d3offset, d3list) &
        !$omp private(i, j, k, ij, nd3pos)
        do i = 1,n
          ij = i*(i-1)/2
          nd3pos = d3offset(i)
          do j = 1,i-1
            k = ij+j
            if (sqrab(k) .lt. dispthr) then
              nd3pos = nd3pos+1
              d3list(1,nd3pos) = i
              d3list(2,nd3pos) = j
            end if
          end do
        end do
        !$omp end parallel do
      end if
    end if
    if (present(frozen_mask)) then
      call work%prepare_static(n,at,xyz,frozen_mask,repthr,dispthr,param,topo,sqrab,srab)
    end if
!      if (pr) call timer%measure(1)

!!!!!!!!!!!!
! Setup HB
!!!!!!!!!!!!

!      if (pr) call timer%measure(10,'HB/XB (incl list setup)')
    !> Preserve an explicit caller request before the legacy count diagnostic
    !> below overwrites the public flag.  Ordinary GFN-FF/MTD calls never set
    !> this request, so their established RMSD/count behavior is unchanged.
    force_hbond_update = nlist%force_hbond_update
    if (allocated(nlist%q)) then
      nlist%initialized = size(nlist%q) == n
    end if
    call gfnff_hbset0(n,at,xyz,sqrab,topo,nhb1,nhb2,nxb,hbthr1,hbthr2)
    nlist%initialized = nlist%initialized.and.nhb1 <= nlist%nhb1 &
       & .and.nhb2 <= nlist%nhb2.and.nxb <= nlist%nxb
    require_update = .not.nlist%initialized
    nlist%force_hbond_update = nhb1 .ne. nlist%nhb1 &
                         & .or.nhb2 .ne. nlist%nhb2 &
                         & .or.nxb .ne. nlist%nxb   &
                         & .or. require_update
    if (.not.nlist%initialized) then
      call new(nlist,n,5*nhb1,5*nhb2,3*nxb)
      nlist%hbrefgeo(:,:) = xyz
    end if
    if (update.or.require_update.or.force_hbond_update) then
      call gfnff_hbset(n,at,xyz,sqrab,topo,nlist,hbthr1,hbthr2, &
         & force_update=force_hbond_update)
    end if
    !> The explicit request is one-shot.  Count differences remain governed
    !> by the release-v2 allocation/RMSD logic and do not leak into MTD calls.
    nlist%force_hbond_update = .false.
!      if (pr) call timer%measure(10)

!!!!!!!!!!!!!
! Setup
! GBSA
!!!!!!!!!!!!!
#ifdef WITH_GBSA
    if (allocated(solvation)) then
!      call timer%measure(11, "GBSA")
      call solvation%update(at,xyz)
!      call timer%measure(11)
    end if
#endif
!!!!!!!!!!!!!
! REP part
! non-bonded
!!!!!!!!!!!!!

!      if (pr) call timer%measure(2,'non bonded repulsion')
    if (work%static_cache_valid) then
      erep = erep+work%cached_nb_rep_total
      if (work%prefix_frozen) then
!> With a frozen prefix, every remaining pair belongs to an active outer atom.
!> Traverse the original packed rows directly to avoid indirect pair-list loads.
        if (serial_inner) then
          gserial = 0.0_wp
          serial_energy = 0.0_wp
          do iat = work%frozen_prefix+1,n
            m = iat*(iat-1)/2
            do jat = 1,iat-1
              ij = m+jat
              r2 = sqrab(ij)
              if (r2 .gt. repthr) cycle
              if (topo%bpair(ij) .eq. 1) cycle
              ati = at(iat)
              atj = at(jat)
              rab = srab(ij)
              t16 = r2**0.75_wp
              t19 = t16*t16
              t8 = t16*topo%alphanb(ij)
              t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
              serial_energy = serial_energy+t26/rab
              t27 = t26*(1.5_wp*t8+1.0_wp)/t19
              r3 = (xyz(:,iat)-xyz(:,jat))*t27
              gserial(:,iat) = gserial(:,iat)-r3
              if (.not.frozen_mask(jat)) gserial(:,jat) = gserial(:,jat)+r3
            end do
          end do
          erep = erep+serial_energy
          g = g+gserial
        else
          !$omp parallel do default(none) schedule(dynamic,8) reduction(+:erep, g) &
          !$omp shared(n, at, xyz, srab, sqrab, repthr, topo, param, frozen_mask, work) &
          !$omp private(iat, jat, m, ij, ati, atj, rab, r2, r3, t8, t16, t19, t26, t27)
          do iat = work%frozen_prefix+1,n
            m = iat*(iat-1)/2
            do jat = 1,iat-1
              ij = m+jat
              r2 = sqrab(ij)
              if (r2 .gt. repthr) cycle
              if (topo%bpair(ij) .eq. 1) cycle
              ati = at(iat)
              atj = at(jat)
              rab = srab(ij)
              t16 = r2**0.75_wp
              t19 = t16*t16
              t8 = t16*topo%alphanb(ij)
              t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
              erep = erep+t26/rab
              t27 = t26*(1.5_wp*t8+1.0_wp)/t19
              r3 = (xyz(:,iat)-xyz(:,jat))*t27
              g(:,iat) = g(:,iat)-r3
              if (.not.frozen_mask(jat)) g(:,jat) = g(:,jat)+r3
            end do
          end do
          !$omp end parallel do
        end if
      else
        if (serial_inner) then
          gserial = 0.0_wp
          serial_energy = 0.0_wp
          do m = 1,work%n_active_pairs
            iat = work%active_pair_i(m)
            jat = work%active_pair_j(m)
            ij = work%active_pair_idx(m)
            r2 = sqrab(ij)
            if (r2 .gt. repthr) cycle
            if (topo%bpair(ij) .eq. 1) cycle
            ati = at(iat)
            atj = at(jat)
            rab = srab(ij)
            t16 = r2**0.75_wp
            t19 = t16*t16
            t8 = t16*topo%alphanb(ij)
            t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
            serial_energy = serial_energy+t26/rab
            t27 = t26*(1.5_wp*t8+1.0_wp)/t19
            r3 = (xyz(:,iat)-xyz(:,jat))*t27
            if (.not.frozen_mask(iat)) gserial(:,iat) = gserial(:,iat)-r3
            if (.not.frozen_mask(jat)) gserial(:,jat) = gserial(:,jat)+r3
          end do
          erep = erep+serial_energy
          g = g+gserial
        else
          !$omp parallel do default(none) schedule(dynamic,32) reduction(+:erep, g) &
          !$omp shared(at, xyz, srab, sqrab, repthr, topo, param, frozen_mask, work) &
          !$omp private(m, iat, jat, ij, ati, atj, rab, r2, r3, t8, t16, t19, t26, t27)
          do m = 1,work%n_active_pairs
            iat = work%active_pair_i(m)
            jat = work%active_pair_j(m)
            ij = work%active_pair_idx(m)
            r2 = sqrab(ij)
            if (r2 .gt. repthr) cycle
            if (topo%bpair(ij) .eq. 1) cycle
            ati = at(iat)
            atj = at(jat)
            rab = srab(ij)
            t16 = r2**0.75_wp
            t19 = t16*t16
            t8 = t16*topo%alphanb(ij)
            t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
            erep = erep+t26/rab
            t27 = t26*(1.5_wp*t8+1.0_wp)/t19
            r3 = (xyz(:,iat)-xyz(:,jat))*t27
            if (.not.frozen_mask(iat)) g(:,iat) = g(:,iat)-r3
            if (.not.frozen_mask(jat)) g(:,jat) = g(:,jat)+r3
          end do
          !$omp end parallel do
        end if
      end if
    else
      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do iat = 1,n
          m = iat*(iat-1)/2
          do jat = 1,iat-1
            ij = m+jat
            r2 = sqrab(ij)
            if (r2 .gt. repthr) cycle
            if (topo%bpair(ij) .eq. 1) cycle
            ati = at(iat)
            atj = at(jat)
            rab = srab(ij)
            t16 = r2**0.75_wp
            t19 = t16*t16
            t8 = t16*topo%alphanb(ij)
            t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
            serial_energy = serial_energy+t26/rab
            if (present(frozen_mask)) then
              if (frozen_mask(iat).and.frozen_mask(jat)) cycle
            end if
            t27 = t26*(1.5_wp*t8+1.0_wp)/t19
            r3 = (xyz(:,iat)-xyz(:,jat))*t27
            if (present(frozen_mask)) then
              if (.not.frozen_mask(iat)) gserial(:,iat) = gserial(:,iat)-r3
              if (.not.frozen_mask(jat)) gserial(:,jat) = gserial(:,jat)+r3
            else
              gserial(:,iat) = gserial(:,iat)-r3
              gserial(:,jat) = gserial(:,jat)+r3
            end if
          end do
        end do
        erep = erep+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) schedule(dynamic,8) reduction(+:erep, g) &
        !$omp shared(n, at, xyz, srab, sqrab, repthr, topo, param, frozen_mask) &
        !$omp private(iat, jat, m, ij, ati, atj, rab, r2, r3, t8, t16, t19, t26, t27)
        do iat = 1,n
          m = iat*(iat-1)/2
          do jat = 1,iat-1
            ij = m+jat
            r2 = sqrab(ij)
            if (r2 .gt. repthr) cycle
            if (topo%bpair(ij) .eq. 1) cycle
            ati = at(iat)
            atj = at(jat)
            rab = srab(ij)
            t16 = r2**0.75_wp
            t19 = t16*t16
            t8 = t16*topo%alphanb(ij)
            t26 = exp(-t8)*param%repz(ati)*param%repz(atj)*param%repscaln
            erep = erep+t26/rab
            if (present(frozen_mask)) then
              if (frozen_mask(iat).and.frozen_mask(jat)) cycle
            end if
            t27 = t26*(1.5_wp*t8+1.0_wp)/t19
            r3 = (xyz(:,iat)-xyz(:,jat))*t27
            if (present(frozen_mask)) then
              if (.not.frozen_mask(iat)) g(:,iat) = g(:,iat)-r3
              if (.not.frozen_mask(jat)) g(:,jat) = g(:,jat)+r3
            else
              g(:,iat) = g(:,iat)-r3
              g(:,jat) = g(:,jat)+r3
            end if
          end do
        end do
        !$omp end parallel do
      end if
    end if
!      if (pr) call timer%measure(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! just a extremely crude mode for 2D-3D conversion
! i.e. an harmonic potential with estimated Re
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (version == gffVersion%harmonic2020) then
      ebond = 0
      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do i = 1,topo%nbond
          iat = topo%blist(1,i)
          jat = topo%blist(2,i)
          r3 = xyz(:,iat)-xyz(:,jat)
          rab = sqrt(sum(r3*r3))
          rn = 0.7*(param%rcov(at(iat))+param%rcov(at(jat)))
          r2 = rn-rab
          serial_energy = serial_energy+0.1d0*r2**2  ! fixfc = 0.1
          dum = 0.1d0*2.0d0*r2/rab
          gserial(:,jat) = gserial(:,jat)+dum*r3
          gserial(:,iat) = gserial(:,iat)-dum*r3
        end do
        ebond = ebond+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) reduction(+:ebond, g) &
        !$omp shared(topo, param, xyz, at) private(i, iat, jat, rab, r2, r3, rn, dum)
        do i = 1,topo%nbond
          iat = topo%blist(1,i)
          jat = topo%blist(2,i)
          r3 = xyz(:,iat)-xyz(:,jat)
          rab = sqrt(sum(r3*r3))
          rn = 0.7*(param%rcov(at(iat))+param%rcov(at(jat)))
          r2 = rn-rab
          ebond = ebond+0.1d0*r2**2  ! fixfc = 0.1
          dum = 0.1d0*2.0d0*r2/rab
          g(:,jat) = g(:,jat)+dum*r3
          g(:,iat) = g(:,iat)-dum*r3
        end do
        !$omp end parallel do
      end if
      etot = ebond+erep
      return
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! erf CN and gradient for disp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if (pr) call timer%measure(3,'dCN')
    if (present(frozen_mask).and.work%frozen_cache_valid) then
      call gfnff_dlogcoord(n,at,xyz,srab,cn,dcn,cnthr,param, &
     & frozen_mask,work%raw_cn_static,work%cn_cache_valid, &
     & work%frozen_prefix,work%prefix_frozen)
      work%cn_cache_thr = cnthr
    else
      call gfnff_dlogcoord(n,at,xyz,srab,cn,dcn,cnthr,param)
    end if
    if (sum(topo%nr_hb) .gt. 0) call dncoord_erf(n,at,xyz,param%rcov,hb_cn,hb_dcn,900.0d0,topo) ! HB erf CN
!      if (pr) call timer%measure(3)

!!!!!!
! EEQ
!!!!!!

!      if (pr) call timer%measure(4,'EEQ energy and q')
    ! The physical single point may reuse the exact validated frozen-host cache.
    call goed_gfnff(.true.,n,at,sqrab,srab,&
   &                dfloat(ichrg),eeqtmp,cn,nlist%q,ees,solvation,param,topo,work,io)  ! without dq/dr
!    if (pr) call timer%measure(4)

!!!!!!!!
!D3(BJ)
!!!!!!!!

!      if (pr) call timer%measure(5,'D3')
    if (nd3 .gt. 0) then
      if (present(frozen_mask)) then
        call d3_gradient(topo%dispm,n,at,xyz,nd3,d3list,topo%zetac6, &
           & param%d3r0,sqrtZr4r2,4.0d0,param%dispscale,cn,dcn,edisp,g, &
           & work%d3_scratch,frozen_mask)
      else
        call d3_gradient(topo%dispm,n,at,xyz,nd3,d3list,topo%zetac6, &
           & param%d3r0,sqrtZr4r2,4.0d0,param%dispscale,cn,dcn,edisp,g, &
           & work%d3_scratch)
      end if
    end if
!      if (pr) call timer%measure(5)

!!!!!!!!
! ES part
!!!!!!!!
!      if (pr) call timer%measure(6,'EEQ gradient')
    if (present(frozen_mask)) then
      nactive = count(.not.frozen_mask)
!> The atom-centred path is intended for large systems with a minority active
!> region.  Retain the original half-pair loop for small or mostly active cases.
      if (n < 256.or.nactive > n/2) then
        if (serial_inner) then
          gserial = 0.0_wp
          do i = 1,n
            k = i*(i-1)/2
            do j = 1,i-1
              if (frozen_mask(i).and.frozen_mask(j)) cycle
              ij = k+j
              r2 = sqrab(ij)
              rab = srab(ij)
              gammij = eeqtmp(1,ij)
              erff = eeqtmp(2,ij)
              dd = (2.0d0*gammij*exp(-gammij**2*r2) &
                 & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
              r3 = (xyz(:,i)-xyz(:,j))*dd
              if (.not.frozen_mask(i)) gserial(:,i) = gserial(:,i)+r3
              if (.not.frozen_mask(j)) gserial(:,j) = gserial(:,j)-r3
            end do
          end do
          g = g+gserial
        else
          !$omp parallel do default(none) reduction (+:g) &
          !$omp shared(nlist,n,sqrab,srab,eeqtmp,xyz,frozen_mask) &
          !$omp private(i,j,k,ij,r3,r2,rab,gammij,erff,dd)
          do i = 1,n
            k = i*(i-1)/2
            do j = 1,i-1
              if (frozen_mask(i).and.frozen_mask(j)) cycle
              ij = k+j
              r2 = sqrab(ij)
              rab = srab(ij)
              gammij = eeqtmp(1,ij)
              erff = eeqtmp(2,ij)
              dd = (2.0d0*gammij*exp(-gammij**2*r2) &
                 & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
              r3 = (xyz(:,i)-xyz(:,j))*dd
              if (.not.frozen_mask(i)) g(:,i) = g(:,i)+r3
              if (.not.frozen_mask(j)) g(:,j) = g(:,j)-r3
            end do
          end do
          !$omp end parallel do
        end if
      else
  !> Evaluate the direct EEQ pair gradient atom-centrically for active atoms.
  !> Each iteration owns one gradient vector, eliminating the full-array OpenMP
  !> reduction. Active-active pairs are visited once from each endpoint; their
  !> two endpoint forces are therefore accumulated without races.
        iact = 0
        do i = 1,n
          if (.not.frozen_mask(i)) then
            iact = iact+1
            active_atoms(iact) = i
          end if
        end do
        if (serial_inner) then
          do iact = 1,nactive
            i = active_atoms(iact)
            gactive = g(:,i)
            do j = 1,n
              if (j == i) cycle
              if (i > j) then
                ij = i*(i-1)/2+j
              else
                ij = j*(j-1)/2+i
              end if
              r2 = sqrab(ij)
              rab = srab(ij)
              gammij = eeqtmp(1,ij)
              erff = eeqtmp(2,ij)
              dd = (2.0d0*gammij*exp(-gammij**2*r2) &
                 & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
              r3 = (xyz(:,i)-xyz(:,j))*dd
              gactive = gactive+r3
            end do
            g(:,i) = gactive
          end do
        else
          !$omp parallel do default(none) schedule(dynamic,4) &
          !$omp shared(nlist,n,sqrab,srab,eeqtmp,xyz,active_atoms,nactive,g) &
          !$omp private(iact,i,j,k,ij,r3,gactive,r2,rab,gammij,erff,dd)
          do iact = 1,nactive
            i = active_atoms(iact)
            gactive = g(:,i)
            do j = 1,n
              if (j == i) cycle
              if (i > j) then
                ij = i*(i-1)/2+j
              else
                ij = j*(j-1)/2+i
              end if
              r2 = sqrab(ij)
              rab = srab(ij)
              gammij = eeqtmp(1,ij)
              erff = eeqtmp(2,ij)
              dd = (2.0d0*gammij*exp(-gammij**2*r2) &
                 & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
              r3 = (xyz(:,i)-xyz(:,j))*dd
              gactive = gactive+r3
            end do
            g(:,i) = gactive
          end do
          !$omp end parallel do
        end if
      end if
    else
      if (serial_inner) then
        gserial = 0.0_wp
        do i = 1,n
          k = i*(i-1)/2
          do j = 1,i-1
            ij = k+j
            r2 = sqrab(ij)
            rab = srab(ij)
            gammij = eeqtmp(1,ij)
            erff = eeqtmp(2,ij)
            dd = (2.0d0*gammij*exp(-gammij**2*r2) &
               & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
            r3 = (xyz(:,i)-xyz(:,j))*dd
            gserial(:,i) = gserial(:,i)+r3
            gserial(:,j) = gserial(:,j)-r3
          end do
        end do
        g = g+gserial
      else
        !$omp parallel do default(none) reduction (+:g) &
        !$omp shared(topo,nlist,n,sqrab,srab,eeqtmp,xyz,at) &
        !$omp private(i,j,k,ij,r3,r2,rab,gammij,erff,dd)
        do i = 1,n
          k = i*(i-1)/2
          do j = 1,i-1
            ij = k+j
            r2 = sqrab(ij)
            rab = srab(ij)
            gammij = eeqtmp(1,ij)
            erff = eeqtmp(2,ij)
            dd = (2.0d0*gammij*exp(-gammij**2*r2) &
               & /(sqrtpi*r2)-erff/(rab*r2))*nlist%q(i)*nlist%q(j)
            r3 = (xyz(:,i)-xyz(:,j))*dd
            g(:,i) = g(:,i)+r3
            g(:,j) = g(:,j)-r3
          end do
        end do
        !$omp end parallel do
      end if
    end if

#ifdef WITH_GBSA
    if (allocated(solvation)) then
!         call timer%measure(11, "GBSA")
      call solvation%addGradient(at,xyz,nlist%q,nlist%q,g)
      call solvation%getEnergyParts(nlist%q,nlist%q,gborn,ghb,gsasa, &
         & gshift)
      gsolv = gsasa+gborn+ghb+gshift
!         call timer%measure(11)
    else
      gborn = 0.0d0
      ghb = 0.0d0
    end if
#endif

    do i = 1,n
      qtmp(i) = nlist%q(i)*param%cnf(at(i))/(2.0d0*sqrt(cn(i))+1.d-16)
    end do

    call gemv(dcn,qtmp,g,alpha=-1.0_wp,beta=1.0_wp)
!      if (pr) call timer%measure(6)

!!!!!!!!!!!!!!!!!!
! SRB bonded part
!!!!!!!!!!!!!!!!!!

!      if (pr) call timer%measure(7,'bonds')
    if (topo%nbond .gt. 0) then
      rab0(:) = topo%vbond(1,:) ! shifts
      call gfnffdrab(n,at,xyz,cn,dcn,topo%nbond,topo%blist,rab0,grab0)

      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do i = 1,topo%nbond
          iat = topo%blist(1,i)
          jat = topo%blist(2,i)
          ati = at(iat)
          atj = at(jat)
          ij = iat*(iat-1)/2+jat
          rab = srab(ij)
          rij = rab0(i)
          drij = grab0(:,:,i)
          if (topo%nr_hb(i) .ge. 1) then
            call egbond_hb(i,iat,jat,rab,rij,drij,hb_cn,hb_dcn,n,at,xyz,serial_energy,gserial,param,topo)
          else
            call egbond(i,iat,jat,rab,rij,drij,n,at,xyz,serial_energy,gserial,topo)
          end if
        end do
        ebond = ebond+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) reduction(+:g, ebond) &
        !$omp shared(grab0, topo, param, rab0, srab, xyz, at, hb_cn, hb_dcn, n) &
        !$omp private(i, k, iat, jat, ij, rab, rij, drij, t8, dr, dum, yy, &
        !$omp& dx, dy, dz, t4, t5, t6, ati, atj)
        do i = 1,topo%nbond
          iat = topo%blist(1,i)
          jat = topo%blist(2,i)
          ati = at(iat)
          atj = at(jat)
          ij = iat*(iat-1)/2+jat
          rab = srab(ij)
          rij = rab0(i)
          drij = grab0(:,:,i)
          if (topo%nr_hb(i) .ge. 1) then
            call egbond_hb(i,iat,jat,rab,rij,drij,hb_cn,hb_dcn,n,at,xyz,ebond,g,param,topo)
          else
            call egbond(i,iat,jat,rab,rij,drij,n,at,xyz,ebond,g,topo)
          end if
        end do
        !$omp end parallel do
      end if


!!!!!!!!!!!!!!!!!!
! bonded REP
!!!!!!!!!!!!!!!!!!

      if (work%static_cache_valid) then
        erep = erep+work%cached_bond_rep_total
        nloop = work%n_dynamic_bonds
      else
        nloop = topo%nbond
      end if
      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do m = 1,nloop
          if (work%static_cache_valid) then
            i = work%dynamic_bond_idx(m)
          else
            i = m
          end if
          iat = topo%blist(1,i)
          jat = topo%blist(2,i)
          ij = iat*(iat-1)/2+jat
          xa = xyz(1,iat)
          ya = xyz(2,iat)
          za = xyz(3,iat)
          dx = xa-xyz(1,jat)
          dy = ya-xyz(2,jat)
          dz = za-xyz(3,jat)
          r2 = sqrab(ij)
          rab = srab(ij)
          ati = at(iat)
          atj = at(jat)
          alpha = sqrt(param%repa(ati)*param%repa(atj))
          repab = param%repz(ati)*param%repz(atj)*param%repscalb
          t16 = r2**0.75_wp
          t19 = t16*t16
          t26 = exp(-alpha*t16)*repab
          serial_energy = serial_energy+t26/rab
          t27 = t26*(1.5_wp*alpha*t16+1.0_wp)/t19
          if (present(frozen_mask)) then
            if (.not.frozen_mask(iat)) then
              gserial(1,iat) = gserial(1,iat)-dx*t27
              gserial(2,iat) = gserial(2,iat)-dy*t27
              gserial(3,iat) = gserial(3,iat)-dz*t27
            end if
            if (.not.frozen_mask(jat)) then
              gserial(1,jat) = gserial(1,jat)+dx*t27
              gserial(2,jat) = gserial(2,jat)+dy*t27
              gserial(3,jat) = gserial(3,jat)+dz*t27
            end if
          else
            gserial(1,iat) = gserial(1,iat)-dx*t27
            gserial(2,iat) = gserial(2,iat)-dy*t27
            gserial(3,iat) = gserial(3,iat)-dz*t27
            gserial(1,jat) = gserial(1,jat)+dx*t27
            gserial(2,jat) = gserial(2,jat)+dy*t27
            gserial(3,jat) = gserial(3,jat)+dz*t27
          end if
        end do
        erep = erep+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) reduction(+:erep, g) &
        !$omp shared(topo, param, at, sqrab, srab, xyz, frozen_mask, work, nloop) &
        !$omp private(m, i, iat, jat, ij, xa, ya, za, dx, dy, dz, r2, rab, ati, atj, &
        !$omp& alpha, repab, t16, t19, t26, t27)
        do m = 1,nloop
          if (work%static_cache_valid) then
            i = work%dynamic_bond_idx(m)
          else
            i = m
          end if
          iat = topo%blist(1,i)
          jat = topo%blist(2,i)
          ij = iat*(iat-1)/2+jat
          xa = xyz(1,iat)
          ya = xyz(2,iat)
          za = xyz(3,iat)
          dx = xa-xyz(1,jat)
          dy = ya-xyz(2,jat)
          dz = za-xyz(3,jat)
          r2 = sqrab(ij)
          rab = srab(ij)
          ati = at(iat)
          atj = at(jat)
          alpha = sqrt(param%repa(ati)*param%repa(atj))
          repab = param%repz(ati)*param%repz(atj)*param%repscalb
          t16 = r2**0.75_wp
          t19 = t16*t16
          t26 = exp(-alpha*t16)*repab
          erep = erep+t26/rab
          t27 = t26*(1.5_wp*alpha*t16+1.0_wp)/t19
          if (present(frozen_mask)) then
            if (.not.frozen_mask(iat)) then
              g(1,iat) = g(1,iat)-dx*t27
              g(2,iat) = g(2,iat)-dy*t27
              g(3,iat) = g(3,iat)-dz*t27
            end if
            if (.not.frozen_mask(jat)) then
              g(1,jat) = g(1,jat)+dx*t27
              g(2,jat) = g(2,jat)+dy*t27
              g(3,jat) = g(3,jat)+dz*t27
            end if
          else
            g(1,iat) = g(1,iat)-dx*t27
            g(2,iat) = g(2,iat)-dy*t27
            g(3,iat) = g(3,iat)-dz*t27
            g(1,jat) = g(1,jat)+dx*t27
            g(2,jat) = g(2,jat)+dy*t27
            g(3,jat) = g(3,jat)+dz*t27
          end if
        end do
        !$omp end parallel do
      end if
    end if
!      if (pr) call timer%measure(7)

!!!!!!!!!!!!!!!!!!
! bend
!!!!!!!!!!!!!!!!!!

!      if (pr) call timer%measure(8,'bend and torsion')
    if (topo%nangl .gt. 0) then
      if (work%static_cache_valid) then
        eangl = eangl+work%cached_angle_total
        nloop = work%n_dynamic_angles
      else
        nloop = topo%nangl
      end if
      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do l = 1,nloop
          if (work%static_cache_valid) then
            m = work%dynamic_angle_idx(l)
          else
            m = l
          end if
          j = topo%alist(1,m)
          i = topo%alist(2,m)
          k = topo%alist(3,m)
          call egbend(m,j,i,k,n,at,xyz,etmp,g3tmp,param,topo)
          gserial(1:3,j) = gserial(1:3,j)+g3tmp(1:3,1)
          gserial(1:3,i) = gserial(1:3,i)+g3tmp(1:3,2)
          gserial(1:3,k) = gserial(1:3,k)+g3tmp(1:3,3)
          serial_energy = serial_energy+etmp
        end do
        eangl = eangl+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) reduction (+:eangl, g) &
        !$omp shared(n, at, xyz, topo, param, work, nloop) &
        !$omp private(m, i, j, k, l, etmp, g3tmp)
        do l = 1,nloop
          if (work%static_cache_valid) then
            m = work%dynamic_angle_idx(l)
          else
            m = l
          end if
          j = topo%alist(1,m)
          i = topo%alist(2,m)
          k = topo%alist(3,m)
          call egbend(m,j,i,k,n,at,xyz,etmp,g3tmp,param,topo)
          g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)
          g(1:3,i) = g(1:3,i)+g3tmp(1:3,2)
          g(1:3,k) = g(1:3,k)+g3tmp(1:3,3)
          eangl = eangl+etmp
        end do
        !$omp end parallel do
      end if
    end if

!!!!!!!!!!!!!!!!!!
! torsion
!!!!!!!!!!!!!!!!!!

    if (topo%ntors .gt. 0) then
      if (work%static_cache_valid) then
        etors = etors+work%cached_torsion_total
        nloop = work%n_dynamic_torsions
      else
        nloop = topo%ntors
      end if
      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do nd3pos = 1,nloop
          if (work%static_cache_valid) then
            m = work%dynamic_torsion_idx(nd3pos)
          else
            m = nd3pos
          end if
          i = topo%tlist(1,m)
          j = topo%tlist(2,m)
          k = topo%tlist(3,m)
          l = topo%tlist(4,m)
          call egtors(m,i,j,k,l,n,at,xyz,etmp,g4tmp,param,topo)
          gserial(1:3,i) = gserial(1:3,i)+g4tmp(1:3,1)
          gserial(1:3,j) = gserial(1:3,j)+g4tmp(1:3,2)
          gserial(1:3,k) = gserial(1:3,k)+g4tmp(1:3,3)
          gserial(1:3,l) = gserial(1:3,l)+g4tmp(1:3,4)
          serial_energy = serial_energy+etmp
        end do
        etors = etors+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) reduction (+:etors, g) &
        !$omp shared(param, topo, n, at, xyz, work, nloop) &
        !$omp private(m, i, j, k, l, nd3pos, etmp, g4tmp)
        do nd3pos = 1,nloop
          if (work%static_cache_valid) then
            m = work%dynamic_torsion_idx(nd3pos)
          else
            m = nd3pos
          end if
          i = topo%tlist(1,m)
          j = topo%tlist(2,m)
          k = topo%tlist(3,m)
          l = topo%tlist(4,m)
          call egtors(m,i,j,k,l,n,at,xyz,etmp,g4tmp,param,topo)
          g(1:3,i) = g(1:3,i)+g4tmp(1:3,1)
          g(1:3,j) = g(1:3,j)+g4tmp(1:3,2)
          g(1:3,k) = g(1:3,k)+g4tmp(1:3,3)
          g(1:3,l) = g(1:3,l)+g4tmp(1:3,4)
          etors = etors+etmp
        end do
        !$omp end parallel do
      end if
    end if
!      if (pr) call timer%measure(8)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! triple bonded carbon torsion potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (allocated(topo%sTorsl)) then
      m = size(topo%sTorsl(1,:))
      if (m.ne.0) then
         do i=1, m
               call sTors_eg(m, n, xyz, topo, etmp, g5tmp)
               etors = etors + etmp
               g = g + g5tmp
         enddo
      endif
   endif

!!!!!!!!!!!!!!!!!!
! BONDED ATM
!!!!!!!!!!!!!!!!!!

!      if (pr) call timer%measure(9,'bonded ATM')
    if (topo%nbatm .gt. 0) then
      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do i = 1,topo%nbatm
          j = topo%b3list(1,i)
          k = topo%b3list(2,i)
          l = topo%b3list(3,i)
          call batmgfnff_eg(n,j,k,l,at,xyz,topo%qa,sqrab,srab,etmp,g3tmp,param)
          gserial(1:3,j) = gserial(1:3,j)+g3tmp(1:3,1)
          gserial(1:3,k) = gserial(1:3,k)+g3tmp(1:3,2)
          gserial(1:3,l) = gserial(1:3,l)+g3tmp(1:3,3)
          serial_energy = serial_energy+etmp
        end do
        ebatm = ebatm+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) reduction(+:ebatm, g) &
        !$omp shared(n, at, xyz, srab, sqrab, topo, param) &
        !$omp private(i, j, k, l, etmp, g3tmp)
        do i = 1,topo%nbatm
          j = topo%b3list(1,i)
          k = topo%b3list(2,i)
          l = topo%b3list(3,i)
          call batmgfnff_eg(n,j,k,l,at,xyz,topo%qa,sqrab,srab,etmp,g3tmp,param)
          g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)
          g(1:3,k) = g(1:3,k)+g3tmp(1:3,2)
          g(1:3,l) = g(1:3,l)+g3tmp(1:3,3)
          ebatm = ebatm+etmp
        end do
        !$omp end parallel do
      end if
    end if
!      if (pr) call timer%measure(9)

!!!!!!!!!!!!!!!!!!
! EHB
!!!!!!!!!!!!!!!!!!

!      if (pr) call timer%measure(10,'HB/XB (incl list setup)')

    if (nlist%nhb1 .gt. 0.or.nlist%nhb2 .gt. 0) then
!> Use one persistent gradient buffer per thread for both hydrogen-bond loops.
!> This removes repeated full-array OpenMP reductions and, for HB2, avoids
!> zeroing and adding a 3*n temporary gradient for every interaction.
      nthreads = 1
!$    nthreads = omp_get_max_threads()
      allocate(ghb_thread(3,n,nthreads),ehb_thread(nthreads))
      ghb_thread = 0.0_wp
      ehb_thread = 0.0_wp

      if (serial_inner) then
        tid = 1
        do i = 1,nlist%nhb1
          j = nlist%hblist1(1,i)
          k = nlist%hblist1(2,i)
          l = nlist%hblist1(3,i)
          call abhgfnff_eg1(n,j,k,l,at,xyz,topo%qa,sqrab,srab,etmp,g3tmp,param,topo)
          ghb_thread(1:3,j,tid) = ghb_thread(1:3,j,tid)+g3tmp(1:3,1)
          ghb_thread(1:3,k,tid) = ghb_thread(1:3,k,tid)+g3tmp(1:3,2)
          ghb_thread(1:3,l,tid) = ghb_thread(1:3,l,tid)+g3tmp(1:3,3)
          ehb_thread(tid) = ehb_thread(tid)+etmp
        end do

        do i = 1,nlist%nhb2
          j = nlist%hblist2(1,i)
          k = nlist%hblist2(2,i)
          l = nlist%hblist2(3,i)
          !Carbonyl case R-C=O...H_A
          if (at(k) .eq. 8.and.topo%nb(20,k) .eq. 1.and.at(topo%nb(1,k)) .eq. 6) then
            call abhgfnff_eg3(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
            !Nitro case R-N=O...H_A
          else if (at(k) .eq. 8.and.topo%nb(20,k) .eq. 1.and.at(topo%nb(1,k)) .eq. 7) then
            call abhgfnff_eg3(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
            !N hetero aromat
          else if (at(k) .eq. 7.and.topo%nb(20,k) .eq. 2) then
            call abhgfnff_eg2_rnr(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
          else
            !Default
            call abhgfnff_eg2new(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
          end if
          ehb_thread(tid) = ehb_thread(tid)+etmp
        end do
      else
        !$omp parallel default(none) &
        !$omp shared(topo, nlist, param, n, at, xyz, sqrab, srab, ghb_thread, ehb_thread) &
        !$omp private(tid, i, j, k, l, etmp, g3tmp)
        tid = 1
!$    tid = omp_get_thread_num()+1

        !$omp do schedule(static)
        do i = 1,nlist%nhb1
          j = nlist%hblist1(1,i)
          k = nlist%hblist1(2,i)
          l = nlist%hblist1(3,i)
          call abhgfnff_eg1(n,j,k,l,at,xyz,topo%qa,sqrab,srab,etmp,g3tmp,param,topo)
          ghb_thread(1:3,j,tid) = ghb_thread(1:3,j,tid)+g3tmp(1:3,1)
          ghb_thread(1:3,k,tid) = ghb_thread(1:3,k,tid)+g3tmp(1:3,2)
          ghb_thread(1:3,l,tid) = ghb_thread(1:3,l,tid)+g3tmp(1:3,3)
          ehb_thread(tid) = ehb_thread(tid)+etmp
        end do
        !$omp end do

        !$omp do schedule(static)
        do i = 1,nlist%nhb2
          j = nlist%hblist2(1,i)
          k = nlist%hblist2(2,i)
          l = nlist%hblist2(3,i)
          !Carbonyl case R-C=O...H_A
          if (at(k) .eq. 8.and.topo%nb(20,k) .eq. 1.and.at(topo%nb(1,k)) .eq. 6) then
            call abhgfnff_eg3(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
            !Nitro case R-N=O...H_A
          else if (at(k) .eq. 8.and.topo%nb(20,k) .eq. 1.and.at(topo%nb(1,k)) .eq. 7) then
            call abhgfnff_eg3(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
            !N hetero aromat
          else if (at(k) .eq. 7.and.topo%nb(20,k) .eq. 2) then
            call abhgfnff_eg2_rnr(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
          else
            !Default
            call abhgfnff_eg2new(n,j,k,l,at,xyz,topo%qa,sqrab,srab, &
               & etmp,ghb_thread(:,:,tid),param,topo)
          end if
          ehb_thread(tid) = ehb_thread(tid)+etmp
        end do
        !$omp end do
        !$omp end parallel
      end if

      do tid = 1,nthreads
        g = g+ghb_thread(:,:,tid)
        ehb = ehb+ehb_thread(tid)
      end do
      deallocate(ghb_thread,ehb_thread)
    end if

!!!!!!!!!!!!!!!!!!
! EXB
!!!!!!!!!!!!!!!!!!

    if (nlist%nxb .gt. 0) then
      if (serial_inner) then
        gserial = 0.0_wp
        serial_energy = 0.0_wp
        do i = 1,nlist%nxb
          j = nlist%hblist3(1,i)
          k = nlist%hblist3(2,i)
          l = nlist%hblist3(3,i)
          call rbxgfnff_eg(n,j,k,l,at,xyz,topo%qa,etmp,g3tmp,param)
          gserial(1:3,j) = gserial(1:3,j)+g3tmp(1:3,1)
          gserial(1:3,k) = gserial(1:3,k)+g3tmp(1:3,2)
          gserial(1:3,l) = gserial(1:3,l)+g3tmp(1:3,3)
          serial_energy = serial_energy+etmp
        end do
        exb = exb+serial_energy
        g = g+gserial
      else
        !$omp parallel do default(none) reduction(+:exb, g) &
        !$omp shared(topo, nlist, param, n, at, xyz) &
        !$omp private(i, j, k, l, etmp, g3tmp)
        do i = 1,nlist%nxb
          j = nlist%hblist3(1,i)
          k = nlist%hblist3(2,i)
          l = nlist%hblist3(3,i)
          call rbxgfnff_eg(n,j,k,l,at,xyz,topo%qa,etmp,g3tmp,param)
          g(1:3,j) = g(1:3,j)+g3tmp(1:3,1)
          g(1:3,k) = g(1:3,k)+g3tmp(1:3,2)
          g(1:3,l) = g(1:3,l)+g3tmp(1:3,3)
          exb = exb+etmp
        end do
        !$omp end parallel do
      end if
    end if
!     if (pr) call timer%measure(10)

!!!!!!!!!!!!!!!!!!
! external stuff
!!!!!!!!!!!!!!!!!!

    if (sum(abs(efield)) .gt. 1d-6) then
      do i = 1,n
        r3(:) = -nlist%q(i)*efield(:)
        g(:,i) = g(:,i)+r3(:)
        eext = eext+r3(1)*(xyz(1,i)-topo%xyze0(1,i))+&
 &                  r3(2)*(xyz(2,i)-topo%xyze0(2,i))+&
 &                  r3(3)*(xyz(3,i)-topo%xyze0(3,i))
      end do
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!
! total energy summation
!!!!!!!!!!!!!!!!!!!!!!!!!!
    etot = ees+edisp+erep+ebond &
   &       +eangl+etors+ehb+exb+ebatm+eext &
   &       +gsolv

!!!!!!!!!!!!!!!!!!
! printout
!!!!!!!!!!!!!!!!!!
    if (pr) then
!        call timer%write(6,'E+G')
      if (abs(sum(nlist%q)-ichrg) .gt. 1.d-1) then ! check EEQ only once
        write (stdout,*) nlist%q
        write (stdout,*) sum(nlist%q),ichrg
        write (stdout,'("EEQ charge constrain error ",a)') source
        io = 1
        return
      end if
      r3 = 0
      do i = 1,n
        r3(:) = r3(:)+nlist%q(i)*xyz(:,i)
      end do

!> just for fit De calc
      sqrab = 1.d+12
      srab = 1.d+6
      cn = 0
!> asymtotically for R=inf, Etot is the SIE contaminted EES
!> which is computed here to get the atomization energy De,n,at(n)
      ! This diagnostic matrix is defined by the supplied infinite-distance
      ! coordinates and must not reuse the physical-geometry frozen-host cache.
      call goed_gfnff(.false.,n,at,sqrab,srab,dfloat(ichrg),eeqtmp,cn,qtmp,eesinf,solvation,param,topo,work,io)
      ! The legacy fit-De diagnostic reuses the physical workspace arrays for
      ! its infinite-distance inputs and EEQ intermediates.  Treat every cache
      ! backed by those arrays as invalid before returning: the next ordinary
      ! call must reconstruct distances, CN data, EEQ pair data, the frozen
      ! block, and the host inverse from its physical geometry.  This avoids
      ! additional O(n**2) diagnostic copies while preventing fit-De state from
      ! being observed by a later cached physical calculation.
      work%frozen_cache_valid = .false.
      work%static_cache_valid = .false.
      work%cn_cache_valid = .false.
      work%eeq_frozen_block_valid = .false.
      work%eeq_host_inverse_valid = .false.
      de = -(etot-eesinf)
    end if
!> write resusts to res type
    res_gff%e_total = etot
    res_gff%gnorm = sqrt(sum(g**2))
    res_gff%e_bond = ebond
    res_gff%e_angl = eangl
    res_gff%e_tors = etors
    res_gff%e_es = ees
    res_gff%e_rep = erep
    res_gff%e_disp = edisp
    res_gff%e_hb = ehb
    res_gff%e_xb = exb
    res_gff%e_batm = ebatm
    res_gff%e_ext = eext
    res_gff%g_hb = ghb
    res_gff%g_born = gborn
    res_gff%g_solv = gsolv
    res_gff%g_shift = gshift
    res_gff%g_sasa = gsasa
    call gemv(xyz,nlist%q,res_gff%dipole)

  end subroutine gfnff_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbond(i,iat,jat,rab,rij,drij,n,at,xyz,e,g,topo)
    implicit none
    !Dummy
    type(TGFFTopology),intent(in) :: topo
    integer,intent(in)   :: i
    integer,intent(in)   :: n
    integer,intent(in)   :: iat
    integer,intent(in)   :: jat
    integer,intent(in)   :: at(n)
    real(wp),intent(in)    :: rab
    real(wp),intent(in)    :: rij
    real(wp),intent(in)    :: drij(3,n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(inout) :: e
    real(wp),intent(inout) :: g(3,n)
    !Stack
    integer j,k
    real(wp) :: dr,dum
    real(wp) :: dx,dy,dz
    real(wp) :: yy
    real(wp) :: t4,t5,t6,t8

    t8 = topo%vbond(2,i)
    dr = rab-rij
    dum = topo%vbond(3,i)*exp(-t8*dr**2)
    e = e+dum                      ! bond energy
    yy = 2.0d0*t8*dr*dum
    dx = xyz(1,iat)-xyz(1,jat)
    dy = xyz(2,iat)-xyz(2,jat)
    dz = xyz(3,iat)-xyz(3,jat)
    t4 = -yy*(dx/rab-drij(1,iat))
    t5 = -yy*(dy/rab-drij(2,iat))
    t6 = -yy*(dz/rab-drij(3,iat))
    g(1,iat) = g(1,iat)+t4-drij(1,iat)*yy ! to avoid if in loop below
    g(2,iat) = g(2,iat)+t5-drij(2,iat)*yy
    g(3,iat) = g(3,iat)+t6-drij(3,iat)*yy
    t4 = -yy*(-dx/rab-drij(1,jat))
    t5 = -yy*(-dy/rab-drij(2,jat))
    t6 = -yy*(-dz/rab-drij(3,jat))
    g(1,jat) = g(1,jat)+t4-drij(1,jat)*yy ! to avoid if in loop below
    g(2,jat) = g(2,jat)+t5-drij(2,jat)*yy
    g(3,jat) = g(3,jat)+t6-drij(3,jat)*yy
    do k = 1,n !3B gradient
      g(:,k) = g(:,k)+drij(:,k)*yy
    end do

  end subroutine egbond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbond_hb(i,iat,jat,rab,rij,drij,hb_cn,hb_dcn,n,at,xyz,e,g,param,topo)
    implicit none
    !Dummy
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer,intent(in)   :: i
    integer,intent(in)   :: n
    integer,intent(in)   :: iat
    integer,intent(in)   :: jat
    integer,intent(in)   :: at(n)
    real(wp),intent(in)    :: rab
    real(wp),intent(in)    :: rij
    real(wp),intent(in)    :: drij(3,n)
    real(wp),intent(in)    :: xyz(3,n)
    real(wp),intent(in)    :: hb_cn(n)
    real(wp),intent(in)    :: hb_dcn(3,n,n)
    real(wp),intent(inout) :: e
    real(wp),intent(inout) :: g(3,n)
    !Stack
    integer j,k
    integer jA,jH
    integer hbH,hbB,hbA
    real(wp) :: dr,dum
    real(wp) :: dx,dy,dz
    real(wp) :: yy,zz
    real(wp) :: t1,t4,t5,t6,t8

    if (at(iat) .eq. 1) then
      hbH = iat
      hbA = jat
    else if (at(jat) .eq. 1) then
      hbH = jat
      hbA = iat
    else
!      write (stdout,'(10x,"No H-atom found in this bond ",i0,1x,i0)') iat,jat
      return
    end if

    t1 = 1.0-param%vbond_scale
    t8 = (-t1*hb_cn(hbH)+1.0)*topo%vbond(2,i)
    dr = rab-rij
    dum = topo%vbond(3,i)*exp(-t8*dr**2)
    e = e+dum                      ! bond energy
    yy = 2.0d0*t8*dr*dum
    dx = xyz(1,iat)-xyz(1,jat)
    dy = xyz(2,iat)-xyz(2,jat)
    dz = xyz(3,iat)-xyz(3,jat)
    t4 = -yy*(dx/rab-drij(1,iat))
    t5 = -yy*(dy/rab-drij(2,iat))
    t6 = -yy*(dz/rab-drij(3,iat))
    g(1,iat) = g(1,iat)+t4-drij(1,iat)*yy ! to avoid if in loop below
    g(2,iat) = g(2,iat)+t5-drij(2,iat)*yy
    g(3,iat) = g(3,iat)+t6-drij(3,iat)*yy
    t4 = -yy*(-dx/rab-drij(1,jat))
    t5 = -yy*(-dy/rab-drij(2,jat))
    t6 = -yy*(-dz/rab-drij(3,jat))
    g(1,jat) = g(1,jat)+t4-drij(1,jat)*yy ! to avoid if in loop below
    g(2,jat) = g(2,jat)+t5-drij(2,jat)*yy
    g(3,jat) = g(3,jat)+t6-drij(3,jat)*yy
    do k = 1,n !3B gradient
      g(:,k) = g(:,k)+drij(:,k)*yy
    end do
    zz = dum*topo%vbond(2,i)*dr**2*t1
    do j = 1,topo%bond_hb_nr !CN gradient
      jH = topo%bond_hb_AH(2,j)
      jA = topo%bond_hb_AH(1,j)
      if (jH .eq. hbH.and.jA .eq. hbA) then
        g(:,hbH) = g(:,hbH)+hb_dcn(:,hbH,hbH)*zz
        do k = 1,topo%bond_hb_Bn(j)
          hbB = topo%bond_hb_B(k,j)
          g(:,hbB) = g(:,hbB)-hb_dcn(:,hbB,hbH)*zz
        end do
      end if
    end do

  end subroutine egbond_hb

  subroutine dncoord_erf(nat,at,xyz,rcov,cn,dcn,thr,topo)
    implicit none
    !Dummy
    type(TGFFTopology),intent(in) :: topo
    integer,intent(in)   :: nat
    integer,intent(in)   :: at(nat)
    real(wp),intent(in)  :: xyz(3,nat)
    real(wp),intent(in)  :: rcov(:)
    real(wp),intent(out) :: cn(nat)
    real(wp),intent(out) :: dcn(3,nat,nat)
    real(wp),intent(in),optional :: thr
    real(wp) :: cn_thr
    !Stack
    integer  :: i,j
    integer  :: lin,linAH
    integer  :: iat,jat
    integer  :: iA,jA,jH
    integer  :: ati,atj
    real(wp) :: r,r2,rij(3)
    real(wp) :: rcovij
    real(wp) :: dtmp,tmp
    real(wp),parameter :: hlfosqrtpi = 1.0_wp/1.77245385091_wp
    real(wp),parameter :: kn = 27.5_wp
    real(wp),parameter :: rcov_scal = 1.78

    cn = 0._wp
    dcn = 0._wp

    do i = 1,topo%bond_hb_nr
      iat = topo%bond_hb_AH(2,i)
      ati = at(iat)
      iA = topo%bond_hb_AH(1,i)
      do j = 1,topo%bond_hb_Bn(i)
        jat = topo%bond_hb_B(j,i)
        atj = at(jat)
        rij = xyz(:,jat)-xyz(:,iat)
        r2 = sum(rij**2)
        if (r2 .gt. thr) cycle
        r = sqrt(r2)
        rcovij = rcov_scal*(rcov(ati)+rcov(atj))
        tmp = 0.5_wp*(1.0_wp+erf(-kn*(r-rcovij)/rcovij))
        dtmp = -hlfosqrtpi*kn*exp(-kn**2*(r-rcovij)**2/rcovij**2)/rcovij
        cn(iat) = cn(iat)+tmp
        cn(jat) = cn(jat)+tmp
        dcn(:,jat,jat) = dtmp*rij/r+dcn(:,jat,jat)
        dcn(:,iat,jat) = dtmp*rij/r
        dcn(:,jat,iat) = -dtmp*rij/r
        dcn(:,iat,iat) = -dtmp*rij/r+dcn(:,iat,iat)
      end do
    end do

  end subroutine dncoord_erf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbend(m,j,i,k,n,at,xyz,e,g,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer m,n,at(n)
    integer i,j,k
    real(wp) :: xyz(3,n),g(3,3),e

    real(wp) ::  c0,kijk,va(3),vb(3),vc(3),cosa
    real(wp) ::  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) ::  term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) ::  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) ::  theta,deda(3),vp(3),et,dij,c1
    real(wp) ::  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ::  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) ::  rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) ::  dampjl,damp2jl,rn

    c0 = topo%vangl(1,m)
    kijk = topo%vangl(2,m)
    va(1:3) = xyz(1:3,i)
    vb(1:3) = xyz(1:3,j)
    vc(1:3) = xyz(1:3,k)
    call vsub(va,vb,vab,3)
    call vsub(vc,vb,vcb,3)
    rab2 = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rcb2 = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    call crprod(vcb,vab,vp)
    rp = vlen(vp)+1.d-14
    call impsc(vab,vcb,cosa)
    cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
    theta = dacos(cosa)

    call gfnffdampa(at(i),at(j),rab2,dampij,damp2ij,param)
    call gfnffdampa(at(k),at(j),rcb2,dampjk,damp2jk,param)
    damp = dampij*dampjk

    if (pi-c0 .lt. 1.d-6) then ! linear
      dt = theta-c0
      ea = kijk*dt**2
      deddt = 2.d0*kijk*dt
    else
      ea = kijk*(cosa-cos(c0))**2
      deddt = 2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
    end if

    e = ea*damp
    call crprod(vab,vp,deda)
    rmul1 = -deddt/(rab2*rp)
    deda = deda*rmul1
    call crprod(vcb,vp,dedc)
    rmul2 = deddt/(rcb2*rp)
    dedc = dedc*rmul2
    dedb = deda+dedc
    term1(1:3) = ea*damp2ij*dampjk*vab(1:3)
    term2(1:3) = ea*damp2jk*dampij*vcb(1:3)
    g(1:3,1) = -dedb(1:3)*damp-term1(1:3)-term2(1:3)
    g(1:3,2) = deda(1:3)*damp+term1(1:3)
    g(1:3,3) = dedc(1:3)*damp+term2(1:3)

  end subroutine egbend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbend_nci_mul(j,i,k,c0,fc,n,at,xyz,e,g)
    implicit none
    !Dummy
    integer n,at(n)
    integer i,j,k
    real(wp) ::  c0,fc
    real(wp) ::  xyz(3,n),g(3,3),e
    !Stack
    real(wp) ::  kijk,va(3),vb(3),vc(3),cosa
    real(wp) ::  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) ::  term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) ::  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) ::  theta,deda(3),vp(3),et,dij,c1
    real(wp) ::  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ::  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) ::  rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) ::  dampjl,damp2jl,rn

    kijk = fc/(cos(0.0d0)-cos(c0))**2
    va(1:3) = xyz(1:3,i)
    vb(1:3) = xyz(1:3,j)
    vc(1:3) = xyz(1:3,k)
    call vsub(va,vb,vab,3)
    call vsub(vc,vb,vcb,3)
    rab2 = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rcb2 = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    call crprod(vcb,vab,vp)
    rp = vlen(vp)+1.d-14
    call impsc(vab,vcb,cosa)
    cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
    theta = dacos(cosa)

    if (pi-c0 .lt. 1.d-6) then     ! linear
      dt = theta-c0
      ea = kijk*dt**2
      deddt = 2.d0*kijk*dt
    else
      ea = kijk*(cosa-cos(c0))**2  ! not linear
      deddt = 2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
    end if

    e = (1.0d0-ea)
    call crprod(vab,vp,deda)
    rmul1 = -deddt/(rab2*rp)
    deda = deda*rmul1
    call crprod(vcb,vp,dedc)
    rmul2 = deddt/(rcb2*rp)
    dedc = dedc*rmul2
    dedb = deda+dedc
    g(1:3,1) = dedb(1:3)
    g(1:3,2) = -deda(1:3)
    g(1:3,3) = -dedc(1:3)

  end subroutine egbend_nci_mul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egbend_nci(j,i,k,c0,kijk,n,at,xyz,e,g,param)
    implicit none
    !Dummy
    type(TGFFData),intent(in) :: param
    integer n,at(n)
    integer i,j,k
    real(wp) :: c0,kijk
    real(wp) :: xyz(3,n),g(3,3),e
    !Stack
    real(wp) ::  va(3),vb(3),vc(3),cosa
    real(wp) ::  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) ::  term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) ::  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) ::  theta,deda(3),vp(3),et,dij,c1
    real(wp) ::  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ::  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) ::  rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) ::  dampjl,damp2jl,rn

    va(1:3) = xyz(1:3,i)
    vb(1:3) = xyz(1:3,j)
    vc(1:3) = xyz(1:3,k)
    call vsub(va,vb,vab,3)
    call vsub(vc,vb,vcb,3)
    rab2 = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rcb2 = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    call crprod(vcb,vab,vp)
    rp = vlen(vp)+1.d-14
    call impsc(vab,vcb,cosa)
    cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
    theta = dacos(cosa)

    call gfnffdampa_nci(at(i),at(j),rab2,dampij,damp2ij,param)
    call gfnffdampa_nci(at(k),at(j),rcb2,dampjk,damp2jk,param)
    damp = dampij*dampjk

    if (pi-c0 .lt. 1.d-6) then ! linear
      dt = theta-c0
      ea = kijk*dt**2
      deddt = 2.d0*kijk*dt
    else
      ea = kijk*(cosa-cos(c0))**2
      deddt = 2.0d0*kijk*sin(theta)*(cos(c0)-cosa)
    end if

    e = ea*damp
    call crprod(vab,vp,deda)
    rmul1 = -deddt/(rab2*rp)
    deda = deda*rmul1
    call crprod(vcb,vp,dedc)
    rmul2 = deddt/(rcb2*rp)
    dedc = dedc*rmul2
    dedb = deda+dedc
    term1(1:3) = ea*damp2ij*dampjk*vab(1:3)
    term2(1:3) = ea*damp2jk*dampij*vcb(1:3)
    g(1:3,1) = -dedb(1:3)*damp-term1(1:3)-term2(1:3)
    g(1:3,2) = deda(1:3)*damp+term1(1:3)
    g(1:3,3) = dedc(1:3)*damp+term2(1:3)

  end subroutine egbend_nci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egtors(m,i,j,k,l,n,at,xyz,e,g,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer m,n,at(n)
    integer i,j,k,l
    real(wp) :: xyz(3,n),g(3,4),e

    real(wp) ::  c0,kijk,va(3),vb(3),vc(3),cosa
    real(wp) ::  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) ::  term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) ::  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) ::  theta,deda(3),vp(3),et,dij,c1
    real(wp) ::  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ::  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) ::  rij,rijk,phi0,rkl,rjk,dampkl,damp2kl
    real(wp) ::  dampjl,damp2jl,rn

    rn = dble(topo%tlist(5,m))
    phi0 = topo%vtors(1,m)
    if (topo%tlist(5,m) .gt. 0) then
      vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
      vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
      vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
      rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
      rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
      rkl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)
      call gfnffdampt(at(i),at(j),rij,dampij,damp2ij,param)
      call gfnffdampt(at(k),at(j),rjk,dampjk,damp2jk,param)
      call gfnffdampt(at(k),at(l),rkl,dampkl,damp2kl,param)
      damp = dampjk*dampij*dampkl
      phi = valijklff(n,xyz,i,j,k,l)
      call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
      dphi1 = phi-phi0
      c1 = rn*dphi1+pi
      x1cos = cos(c1)
      x1sin = sin(c1)
      et = (1.+x1cos)*topo%vtors(2,m)
      dij = -rn*x1sin*topo%vtors(2,m)*damp
      term1(1:3) = et*damp2ij*dampjk*dampkl*vab(1:3)
      term2(1:3) = et*damp2jk*dampij*dampkl*vcb(1:3)
      term3(1:3) = et*damp2kl*dampij*dampjk*vdc(1:3)
      g(1:3,1) = dij*dda(1:3)+term1
      g(1:3,2) = dij*ddb(1:3)-term1+term2
      g(1:3,3) = dij*ddc(1:3)+term3-term2
      g(1:3,4) = dij*ddd(1:3)-term3
      e = et*damp
    else
      vab(1:3) = xyz(1:3,j)-xyz(1:3,i)
      vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
      vdc(1:3) = xyz(1:3,j)-xyz(1:3,l)
      rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
      rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
      rjl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)
      call gfnffdampt(at(i),at(j),rij,dampij,damp2ij,param)
      call gfnffdampt(at(k),at(j),rjk,dampjk,damp2jk,param)
      call gfnffdampt(at(j),at(l),rjl,dampjl,damp2jl,param)
      damp = dampjk*dampij*dampjl
      phi = omega(n,xyz,i,j,k,l)
      call domegadr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
      if (topo%tlist(5,m) .eq. 0) then  ! phi0=0 case
        dphi1 = phi-phi0
        c1 = dphi1+pi
        x1cos = cos(c1)
        x1sin = sin(c1)
        et = (1.+x1cos)*topo%vtors(2,m)
        dij = -x1sin*topo%vtors(2,m)*damp
      else                     ! double min at phi0,-phi0
        et = topo%vtors(2,m)*(cos(phi)-cos(phi0))**2
        dij = 2.*topo%vtors(2,m)*sin(phi)*(cos(phi0)-cos(phi))*damp
      end if
      term1(1:3) = et*damp2ij*dampjk*dampjl*vab(1:3)
      term2(1:3) = et*damp2jk*dampij*dampjl*vcb(1:3)
      term3(1:3) = et*damp2jl*dampij*dampjk*vdc(1:3)
      g(1:3,1) = dij*dda(1:3)-term1
      g(1:3,2) = dij*ddb(1:3)+term1+term2+term3
      g(1:3,3) = dij*ddc(1:3)-term2
      g(1:3,4) = dij*ddd(1:3)-term3
      e = et*damp
    end if

  end subroutine egtors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !torsion without distance damping!!! damping is inherint in the HB term
  subroutine egtors_nci_mul(i,j,k,l,rn,phi0,tshift,n,at,xyz,e,g)
    implicit none
    !Dummy
    integer n,at(n)
    integer i,j,k,l
    integer rn
    real(wp) :: phi0,tshift
    real(wp) :: xyz(3,n),g(3,4),e
    !Stack
    real(wp) ::  c0,fc,kijk,va(3),vb(3),vc(3),cosa
    real(wp) ::  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) ::  term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) ::  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) ::  theta,deda(3),vp(3),et,dij,c1
    real(wp) ::  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ::  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) ::  rij,rijk,rkl,rjk,dampkl,damp2kl
    real(wp) ::  dampjl,damp2jl

    fc = (1.0d0-tshift)/2.0d0
    vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
    vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
    vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
    rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    rkl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)
    phi = valijklff(n,xyz,i,j,k,l)
    call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
    dphi1 = phi-phi0
    c1 = rn*dphi1+pi
    x1cos = cos(c1)
    x1sin = sin(c1)
    et = (1.+x1cos)*fc+tshift
    dij = -rn*x1sin*fc
    g(1:3,1) = dij*dda(1:3)
    g(1:3,2) = dij*ddb(1:3)
    g(1:3,3) = dij*ddc(1:3)
    g(1:3,4) = dij*ddd(1:3)
    e = et !*damp
  end subroutine egtors_nci_mul

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine egtors_nci(i,j,k,l,rn,phi0,fc,n,at,xyz,e,g,param)
    implicit none
    !Dummy
    type(TGFFData),intent(in) :: param
    integer n,at(n)
    integer i,j,k,l
    integer rn
    real(wp) :: phi0,fc
    real(wp) :: xyz(3,n),g(3,4),e
    !Stack
    real(wp) ::  c0,kijk,va(3),vb(3),vc(3),cosa
    real(wp) ::  dt,ea,dedb(3),dedc(3),rmul2,rmul1,deddt
    real(wp) ::  term1(3),term2(3),rab2,vab(3),vcb(3),rp
    real(wp) ::  rcb2,damp,dampij,damp2ij,dampjk,damp2jk
    real(wp) ::  theta,deda(3),vp(3),et,dij,c1
    real(wp) ::  term3(3),x1sin,x1cos,e1,dphi1,vdc(3)
    real(wp) ::  ddd(3),ddc(3),ddb(3),dda(3),rjl,phi
    real(wp) ::  rij,rijk,rkl,rjk,dampkl,damp2kl
    real(wp) ::  dampjl,damp2jl

    vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
    vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
    vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
    rij = vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
    rjk = vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
    rkl = vdc(1)*vdc(1)+vdc(2)*vdc(2)+vdc(3)*vdc(3)
    call gfnffdampt_nci(at(i),at(j),rij,dampij,damp2ij,param)
    call gfnffdampt_nci(at(k),at(j),rjk,dampjk,damp2jk,param)
    call gfnffdampt_nci(at(k),at(l),rkl,dampkl,damp2kl,param)
    damp = dampjk*dampij*dampkl
    phi = valijklff(n,xyz,i,j,k,l)
    call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
    dphi1 = phi-phi0
    c1 = rn*dphi1+pi
    x1cos = cos(c1)
    x1sin = sin(c1)
    et = (1.+x1cos)*fc
    dij = -rn*x1sin*fc*damp
    term1(1:3) = et*damp2ij*dampjk*dampkl*vab(1:3)
    term2(1:3) = et*damp2jk*dampij*dampkl*vcb(1:3)
    term3(1:3) = et*damp2kl*dampij*dampjk*vdc(1:3)
    g(1:3,1) = dij*dda(1:3)+term1
    g(1:3,2) = dij*ddb(1:3)-term1+term2
    g(1:3,3) = dij*ddc(1:3)+term3-term2
    g(1:3,4) = dij*ddd(1:3)-term3
    e = et*damp
  end subroutine egtors_nci

!cccccccccccccccccccccccccccccccccccccccccccccc
! damping of bend and torsion for long
! bond distances to allow proper dissociation
!cccccccccccccccccccccccccccccccccccccccccccccc

  subroutine gfnffdampa(ati,atj,r2,damp,ddamp,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) :: r2,damp,ddamp,rr,rcut
    rcut = param%atcuta*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)
  end subroutine gfnffdampa

  subroutine gfnffdampt(ati,atj,r2,damp,ddamp,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) :: r2,damp,ddamp,rr,rcut
    rcut = param%atcutt*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)
  end subroutine gfnffdampt

  subroutine gfnffdampa_nci(ati,atj,r2,damp,ddamp,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) :: r2,damp,ddamp,rr,rcut
    rcut = param%atcuta_nci*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)
  end subroutine gfnffdampa_nci

  subroutine gfnffdampt_nci(ati,atj,r2,damp,ddamp,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer ati,atj
    real(wp) :: r2,damp,ddamp,rr,rcut
    rcut = param%atcutt_nci*(param%rcov(ati)+param%rcov(atj))**2
    rr = (r2/rcut)**2
    damp = 1.0d0/(1.0d0+rr)
    ddamp = -2.d0*2*rr/(r2*(1.0d0+rr)**2)
  end subroutine gfnffdampt_nci

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ref.: S. Alireza Ghasemi, Albert Hofstetter, Santanu Saha, and Stefan Goedecker
!       PHYSICAL REVIEW B 92, 045131 (2015)
!       Interatomic potentials for ionic systems with density functional accuracy
!       based on charge densities obtained by a neural network
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine goed_gfnff(allow_frozen_cache,n,at,sqrab,r,chrg,eeqtmp,cn,q,es,gbsa,param,topo,work,io)
    implicit none
    character(len=*),parameter :: source = 'gfnff_eg_goed'
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    type(gfnff_workspace),intent(inout),target :: work
    logical,intent(in)  :: allow_frozen_cache
    integer,intent(in)  :: n          !> number of atoms
    integer,intent(in)  :: at(n)      !> ordinal numbers
    real(wp),intent(in)  :: sqrab(n*(n+1)/2)   !> squared dist
    real(wp),intent(in)  :: r(n*(n+1)/2)       !> dist
    real(wp),intent(in)  :: chrg       !> total charge on system
    real(wp),intent(in)  :: cn(n)      !> CN
    real(wp),intent(out) :: q(n)       !> output charges
    real(wp),intent(out) :: es         !> ES energy
    real(wp),intent(inout),target :: eeqtmp(2,n*(n+1)/2)   !> intermediates
    type(TBorn),allocatable,intent(in) :: gbsa    !> solvation object
    integer,intent(out) :: io   !> return status
    !> LOCAL
    integer :: m,i,j,k,ii,ij,nf
    integer :: io1,io2,io_block
    real(wp) :: gammij,tsqrt2pi,tmp
    real(wp),pointer :: A(:,:),x(:)
    integer,pointer :: ipiv(:)
!>  parameter
    parameter(tsqrt2pi=0.797884560802866_wp)
    logical :: exitRun,use_block_solver,lazy_block_assembly,serial_inner

    io = 0  !> return status
    serial_inner = .true.
!$  serial_inner = omp_get_max_threads() == 1

!> # atoms + fragment charge constraints
    m = n+topo%nfrag
    A => work%eeq_a
    x => work%eeq_x
    ipiv => work%eeq_ipiv

!> setup RHS
    do i = 1,n
      x(i) = topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))
    end do

    if (allow_frozen_cache) call work%prepare_eeq_frozen_block(n,r,topo)

    nf = 0
    if (allow_frozen_cache .and. work%eeq_frozen_block_valid) nf = work%frozen_prefix
    lazy_block_assembly = allow_frozen_cache .and. work%eeq_frozen_block_valid .and. &
   &                      nf >= 1 .and. nf < n .and. &
   &                      n-nf+topo%nfrag < nf
#ifdef WITH_GBSA
    lazy_block_assembly = lazy_block_assembly .and. .not.allocated(gbsa)
#endif

    if (lazy_block_assembly) then
      ! The block solver reads only K and D. Active atomic entries below are
      ! fully overwritten; only the fragment rows/columns require clearing.
      if (topo%nfrag > 0) then
        A(1:n,n+1:m) = 0.0_wp
        A(n+1:m,1:m) = 0.0_wp
      end if
    else
      ! Retain the complete legacy assembly for every non-block route.
      A = 0.0_wp
    end if

!> setup A matrix. For a contiguous frozen prefix, retain the exact
!> frozen-frozen block and rebuild only rows that contain active atoms.
    if (allow_frozen_cache .and. work%eeq_frozen_block_valid) then
      nf = work%frozen_prefix
      if (.not.lazy_block_assembly) then
        A(1:nf,1:nf) = work%eeq_frozen_block(1:nf,1:nf)
      end if
      if (serial_inner) then
        do i = nf+1,n
          A(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i)
          k = i*(i-1)/2
          do j = 1,i-1
            ij = k+j
            gammij = 1.0_wp/sqrt(topo%alpeeq(i)+topo%alpeeq(j))
            tmp = erf(gammij*r(ij))
            eeqtmp(1,ij) = gammij
            eeqtmp(2,ij) = tmp
            A(j,i) = tmp/r(ij)
            A(i,j) = A(j,i)
          end do
        end do
      else
        !$omp parallel default(none) &
        !$omp shared(topo,n,nf,r,eeqtmp,A) &
        !$omp private(i,j,k,ij,gammij,tmp)
        !$omp do schedule(dynamic)
        do i = nf+1,n
          A(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i)
          k = i*(i-1)/2
          do j = 1,i-1
            ij = k+j
            gammij = 1.0_wp/sqrt(topo%alpeeq(i)+topo%alpeeq(j))
            tmp = erf(gammij*r(ij))
            eeqtmp(1,ij) = gammij
            eeqtmp(2,ij) = tmp
            A(j,i) = tmp/r(ij)
            A(i,j) = A(j,i)
          end do
        end do
        !$omp enddo
        !$omp end parallel
      end if
    else
      if (serial_inner) then
        do i = 1,n
          A(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i)
          k = i*(i-1)/2
          do j = 1,i-1
            ij = k+j
            gammij = 1.0_wp/sqrt(topo%alpeeq(i)+topo%alpeeq(j))
            tmp = erf(gammij*r(ij))
            eeqtmp(1,ij) = gammij
            eeqtmp(2,ij) = tmp
            A(j,i) = tmp/r(ij)
            A(i,j) = A(j,i)
          end do
        end do
      else
        !$omp parallel default(none) &
        !$omp shared(topo,n,r,eeqtmp,A) &
        !$omp private(i,j,k,ij,gammij,tmp)
        !$omp do schedule(dynamic)
        do i = 1,n
          A(i,i) = tsqrt2pi/sqrt(topo%alpeeq(i))+topo%gameeq(i)
          k = i*(i-1)/2
          do j = 1,i-1
            ij = k+j
            gammij = 1.0_wp/sqrt(topo%alpeeq(i)+topo%alpeeq(j))
            tmp = erf(gammij*r(ij))
            eeqtmp(1,ij) = gammij
            eeqtmp(2,ij) = tmp
            A(j,i) = tmp/r(ij)
            A(i,j) = A(j,i)
          end do
        end do
        !$omp enddo
        !$omp end parallel
      end if
    end if

!> fragment charge constraints
    do i = 1,topo%nfrag
      x(n+i) = topo%qfrag(i)
      do j = 1,n
        if (topo%fraglist(j) .eq. i) then
          A(n+i,j) = 1.0_wp
          A(j,n+i) = 1.0_wp
        end if
      end do
    end do

#ifdef WITH_GBSA
    if (allocated(gbsa)) then
      A(:n,:n) = A(:n,:n)+gbsa%bornMat(:,:)
    end if
#endif

    io_block = 0
    use_block_solver = allow_frozen_cache
#ifdef WITH_GBSA
    use_block_solver = use_block_solver .and. .not.allocated(gbsa)
#endif
    if (use_block_solver) then
      call work%prepare_eeq_host_inverse(n,topo,io_block)
    end if
    use_block_solver = use_block_solver .and. io_block == 0 &
   &                   .and. work%eeq_host_inverse_valid
    if (use_block_solver) then
      call solve_eeq_host_block(n,topo%nfrag,work%frozen_prefix,topo%fraglist, &
     &                         topo%qfrag,A,x,work,q,io1,io2)
    else
      io1 = io_block
      io2 = 0
    end if
    if (.not.use_block_solver .or. io1 /= 0 .or. io2 /= 0) then
      if (lazy_block_assembly) then
        ! A lazy successful block solve never reads H. A full-system
        ! fallback does, so restore the pristine cached host block first.
        A(1:nf,1:nf) = work%eeq_frozen_block(1:nf,1:nf)
      end if
      call sytrf_cached_wrap(A,ipiv,work%eeq_lapack_work,io1)
      call sytrs_wrap(A,x,ipiv,io2)
      q(1:n) = x(1:n)
    end if

    exitRun = (io1 /= 0).or.(io2 /= 0)
    if (exitRun) then
      write (stdout,'("Solving linear equations failed ",a)') source
      io = max(abs(io1),abs(io2))
      return
    end if

    if (n .eq. 1) q(1) = chrg

!> energy
    es = 0.0_wp
    do i = 1,n
      ii = i*(i-1)/2
      do j = 1,i-1
        ij = ii+j
        tmp = eeqtmp(2,ij)
        es = es+q(i)*q(j)*tmp/r(ij)
      end do
      es = es-q(i)*(topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))) &
     &        +q(i)*q(i)*0.5_wp*(topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i)))
    end do

  end subroutine goed_gfnff

  subroutine solve_eeq_host_block(n,nfrag,nf,fraglist,qfrag,A,x,work,q,io1,io2)
    integer,intent(in) :: n,nfrag,nf,fraglist(n)
    real(wp),intent(in) :: qfrag(nfrag),A(n+nfrag,n+nfrag),x(n+nfrag)
    type(gfnff_workspace),intent(inout) :: work
    real(wp),intent(out) :: q(n)
    integer,intent(out) :: io1,io2
    integer :: na,p,i,j,k,ia

    na = n-nf
    p = na+nfrag
    io1 = 0
    io2 = 0
    if (.not.work%eeq_host_inverse_valid) then
      io1 = -1
      return
    end if
    if (size(work%eeq_block_s,1) /= p) then
      io1 = -2
      return
    end if

    ! K couples the frozen host block to active atomic charges and all
    ! fragment-charge multipliers. The fragment columns are fixed 0/1 entries.
    do j = 1,na
      do i = 1,nf
        work%eeq_block_k(i,j) = A(i,nf+j)
      end do
    end do
    if (nfrag > 0) work%eeq_block_k(:,na+1:p) = 0.0_wp
    do i = 1,nf
      k = fraglist(i)
      if (k >= 1 .and. k <= nfrag) work%eeq_block_k(i,na+k) = 1.0_wp
    end do

    ! Y = H^{-1} K and y0 = H^{-1} b_H, using the cached exact inverse of
    ! the invariant frozen-host EEQ block H.
    ! GEMM's default beta=0 fully overwrites Y.
    call gemm(work%eeq_host_inverse,work%eeq_block_k,work%eeq_block_y)
    work%eeq_host_rhs = x(1:nf)
    work%eeq_host_solution = 0.0_wp
    call gemv(work%eeq_host_inverse,work%eeq_host_rhs,work%eeq_host_solution)

    ! Dynamic Schur complement S = D - K^T H^{-1} K.
    do j = 1,na
      do i = 1,na
        work%eeq_block_s(i,j) = A(nf+i,nf+j)
      end do
    end do
    if (nfrag > 0) then
      ! The active block is fully overwritten; initialize only constraint
      ! rows and columns before inserting their exact 0/1 entries.
      work%eeq_block_s(1:na,na+1:p) = 0.0_wp
      work%eeq_block_s(na+1:p,1:p) = 0.0_wp
    end if
    do ia = 1,na
      k = fraglist(nf+ia)
      if (k >= 1 .and. k <= nfrag) then
        work%eeq_block_s(ia,na+k) = 1.0_wp
        work%eeq_block_s(na+k,ia) = 1.0_wp
      end if
    end do
    call gemm(work%eeq_block_k,work%eeq_block_y,work%eeq_block_s, &
   &          transa='t',alpha=-1.0_wp,beta=1.0_wp)

    work%eeq_block_rhs(1:na) = x(nf+1:n)
    if (nfrag > 0) work%eeq_block_rhs(na+1:p) = qfrag
    call gemv(work%eeq_block_k,work%eeq_host_solution,work%eeq_block_rhs, &
   &          alpha=-1.0_wp,beta=1.0_wp,trans='t')
    work%eeq_block_z = work%eeq_block_rhs

    call sytrf_cached_wrap(work%eeq_block_s,work%eeq_block_ipiv, &
   &                       work%eeq_block_fact_work,io1)
    if (io1 /= 0) return
    call sytrs_wrap(work%eeq_block_s,work%eeq_block_z,work%eeq_block_ipiv,io2)
    if (io2 /= 0) return

    q(nf+1:n) = work%eeq_block_z(1:na)
    q(1:nf) = work%eeq_host_solution
    call gemv(work%eeq_block_y,work%eeq_block_z,q(1:nf), &
   &          alpha=-1.0_wp,beta=1.0_wp)
  end subroutine solve_eeq_host_block

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! HB energy and analytical gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!> subroutine for case 1: A...H...B
  subroutine abhgfnff_eg1(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer A,B,H,n,at(n)
    real(wp) :: xyz(3,n),energy,gdr(3,3)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp) :: srab(n*(n+1)/2)    ! dist

    real(wp) :: outl,dampl,damps,rdamp,damp,dd24a,dd24b
    real(wp) :: ratio1,ratio2,ratio3
    real(wp) :: xm,ym,zm
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) :: drah(3),drbh(3),drab(3),drm(3)
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3)
    real(wp) :: ga(3),gb(3),gh(3)
    real(wp) :: gi,denom,ratio,tmp,qhoutl,radab,rahprbh
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo
    real(wp) :: bas,aci
    real(wp) :: eabh
    real(wp) :: aterm,rterm,dterm,sterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: gqa,gqb,gqh
    real(wp) :: caa,cbb
    real(wp) :: shortcut

    integer i,j,ij,lina
    lina(i,j) = min(i,j)+max(i,j)*(max(i,j)-1)/2

    gdr = 0
    energy = 0

    call hbonds(A,B,ca,cb,param,topo)

!     A-B distance
    ij = lina(A,B)
    rab2 = sqrab(ij)
    rab = srab(ij)

!     A-H distance
    ij = lina(A,H)
    rah2 = sqrab(ij)
    rah = srab(ij)

!     B-H distance
    ij = lina(B,H)
    rbh2 = sqrab(ij)
    rbh = srab(ij)

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!     out-of-line damp
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!     long damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!     short damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    rdamp = damp/rab2/rab

!     hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!     hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!     hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

!     donor-acceptor term
    rah4 = rah2*rah2
    rbh4 = rbh2*rbh2
    denom = 1.d0/(rah4+rbh4)

    caa = qa*ca(1)
    cbb = qb*cb(1)
    qhoutl = qh*outl

    bas = (caa*rah4+cbb*rbh4)*denom
    aci = (cb(2)*rah4+ca(2)*rbh4)*denom

!     energy
    rterm = -aci*rdamp*qhoutl
    energy = bas*rterm

!     gradient
    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)
    drab(1:3) = xyz(1:3,A)-xyz(1:3,B)

    aterm = -aci*bas*rdamp*qh
    sterm = -rdamp*bas*qhoutl
    dterm = -aci*bas*qhoutl

    tmp = denom*denom*4.0d0
    dd24a = rah2*rbh4*tmp
    dd24b = rbh2*rah4*tmp

!     donor-acceptor part: bas
    gi = (caa-cbb)*dd24a*rterm
    ga(1:3) = gi*drah(1:3)
    gi = (cbb-caa)*dd24b*rterm
    gb(1:3) = gi*drbh(1:3)
    gh(1:3) = -ga(1:3)-gb(1:3)

!     donor-acceptor part: aci
    gi = (cb(2)-ca(2))*dd24a
    dga(1:3) = gi*drah(1:3)*sterm
    ga(1:3) = ga(1:3)+dga(1:3)

    gi = (ca(2)-cb(2))*dd24b
    dgb(1:3) = gi*drbh(1:3)*sterm
    gb(1:3) = gb(1:3)+dgb(1:3)

    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!     damping part rab
    gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1+ratio1))+(2.d0*param%hbalp*ratio3/(1+ratio3))-3.d0)/rab2
    dg(1:3) = gi*drab(1:3)*dterm
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rab
    gi = aterm*2.d0*ratio2*expo*rahprbh/(1+ratio2)**2/(rahprbh-rab)/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rah,rbh
    tmp = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    dga(1:3) = drah(1:3)*tmp/rah
    ga(1:3) = ga(1:3)+dga(1:3)
    dgb(1:3) = drbh(1:3)*tmp/rbh
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!     move gradients into place
    gdr(1:3,1) = ga(1:3)
    gdr(1:3,2) = gb(1:3)
    gdr(1:3,3) = gh(1:3)

  end subroutine abhgfnff_eg1

!subroutine for case 2: A-H...B including orientation of neighbors at B
  subroutine abhgfnff_eg2new(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer A,B,H,n,at(n)
    real(wp) :: xyz(3,n),energy,gdr(3,n)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp) :: srab(n*(n+1)/2)    ! dist

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ddamp,rabdamp,rbhdamp
    real(wp) :: ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
    real(wp) :: xm,ym,zm
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) :: ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
    real(wp) :: drah(3),drbh(3),drab(3),drm(3)
    real(wp) :: dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) :: ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
    real(wp) :: denom,ratio,qhoutl,radab
    real(wp) :: gi,gi_nb(topo%nb(20,B))
    real(wp) :: tmp1,tmp2(topo%nb(20,B))
    real(wp) :: rahprbh,ranbprbnb(topo%nb(20,B))
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
    real(wp) :: eabh
    real(wp) :: aterm,dterm,nbterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: gqa,gqb,gqh
    real(wp) :: shortcut
    real(wp) :: const
    real(wp) :: outl_nb(topo%nb(20,B)),outl_nb_tot
    real(wp) :: hbnbcut_save
    logical mask_nb(topo%nb(20,B))

!     proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer i,j,ij,lina,nbb
    lina(i,j) = min(i,j)+max(i,j)*(max(i,j)-1)/2

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    energy = 0

    call hbonds(A,B,ca,cb,param,topo)

    nbb = topo%nb(20,B)
!     Neighbours of B
    do i = 1,nbb
!        compute distances
      dranb(1:3,i) = xyz(1:3,A)-xyz(1:3,topo%nb(i,B))
      drbnb(1:3,i) = xyz(1:3,B)-xyz(1:3,topo%nb(i,B))
!        A-nb(B) distance
      ranb2(i) = sum(dranb(1:3,i)**2)
      ranb(i) = sqrt(ranb2(i))
!        B-nb(B) distance
      rbnb2(i) = sum(drbnb(1:3,i)**2)
      rbnb(i) = sqrt(rbnb2(i))
    end do

!     A-B distance
    ij = lina(A,B)
    rab2 = sqrab(ij)
    rab = srab(ij)

!     A-H distance
    ij = lina(A,H)
    rah2 = sqrab(ij)
    rah = srab(ij)

!     B-H distance
    ij = lina(B,H)
    rbh2 = sqrab(ij)
    rbh = srab(ij)

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!     out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!     out-of-line damp: A...nb(B)-B
    if (at(B) .eq. 7.and.topo%nb(20,B) .eq. 1) then
      hbnbcut_save = 2.0
    else
      hbnbcut_save = param%hbnbcut
    end if
    do i = 1,nbb
      ranbprbnb(i) = ranb(i)+rbnb(i)+1.d-12
      expo_nb(i) = (hbnbcut_save/radab)*(ranbprbnb(i)/rab-1.d0)
      ratio2_nb(i) = exp(-expo_nb(i))**(1.0)
      outl_nb(i) = (2.d0/(1.d0+ratio2_nb(i)))-1.0d0
    end do
    outl_nb_tot = product(outl_nb)

!     long damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!     short damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

!     hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!     hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!     hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot

!     constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!     energy
    energy = -rdamp*qhoutl*const

!     gradient
    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)
    drab(1:3) = xyz(1:3,A)-xyz(1:3,B)

    aterm = -rdamp*qh*outl_nb_tot*const
    nbterm = -rdamp*qh*outl*const
    dterm = -qhoutl*const

!------------------------------------------------------------------------------
!     damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!     angular A-H...B term
!------------------------------------------------------------------------------
!     out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!     out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!------------------------------------------------------------------------------
!     angular A...nb(B)-B term
!------------------------------------------------------------------------------
!     out of line term: rab
    mask_nb = .true.
    do i = 1,nbb
      mask_nb(i) = .false.
      tmp2(i) = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
               & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
      gi_nb(i) = -tmp2(i)*ranbprbnb(i)/rab2
      dg(1:3) = gi_nb(i)*drab(1:3)
      ga(1:3) = ga(1:3)+dg(1:3)
      gb(1:3) = gb(1:3)-dg(1:3)
      mask_nb = .true.
    end do

!     out of line term: ranb,rbnb
    do i = 1,nbb
      gi_nb(i) = tmp2(i)/ranb(i)
      dga(1:3) = gi_nb(i)*dranb(1:3,i)
      ga(1:3) = ga(1:3)+dga(1:3)
      gi_nb(i) = tmp2(i)/rbnb(i)
      dgb(1:3) = gi_nb(i)*drbnb(1:3,i)
      gb(1:3) = gb(1:3)+dgb(1:3)
      dgnb(1:3) = -dga(1:3)-dgb(1:3)
      gnb(1:3,i) = dgnb(1:3)
    end do

!------------------------------------------------------------------------------
    if (nbb .lt. 1) then
      gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
      gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
      gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
      return
    end if

!------------------------------------------------------------------------------
!     move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    do i = 1,nbb
      gdr(1:3,topo%nb(i,B)) = gdr(1:3,topo%nb(i,B))+gnb(1:3,i)
    end do

  end subroutine abhgfnff_eg2new

!> subroutine for case 2: A-H...B including LP position
  subroutine abhgfnff_eg2_rnr(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer A,B,H,n,at(n)
    real(wp) :: xyz(3,n),energy,gdr(3,n)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   !> squared dist
    real(wp) :: srab(n*(n+1)/2)    !> dist

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ddamp,rabdamp,rbhdamp
    real(wp) :: ratio1,ratio2,ratio2_lp,ratio2_nb(topo%nb(20,B)),ratio3
    real(wp) :: xm,ym,zm
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) :: ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
    real(wp) :: drah(3),drbh(3),drab(3),drm(3),dralp(3),drblp(3)
    real(wp) :: dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) :: ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B)),gnb_lp(3),glp(3)
    real(wp) :: denom,ratio,qhoutl,radab
    real(wp) :: gi,gi_nb(topo%nb(20,B))
    real(wp) :: tmp1,tmp2(topo%nb(20,B)),tmp3
    real(wp) :: rahprbh,ranbprbnb(topo%nb(20,B))
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_lp,expo_nb(topo%nb(20,B))
    real(wp) :: eabh
    real(wp) :: aterm,dterm,nbterm,lpterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: gqa,gqb,gqh
    real(wp) :: shortcut
    real(wp) :: const
    real(wp) :: outl_nb(topo%nb(20,B)),outl_nb_tot,outl_lp
    real(wp) :: vector(3),vnorm
    real(wp) :: gii(3,3)
    real(wp) :: unit_vec(3)
    real(wp) :: drnb(3,topo%nb(20,B))
    real(wp) :: lp(3)   !> lonepair position
    real(wp) :: lp_dist !> distance parameter between B and lonepair
    real(wp) :: ralp,ralp2,rblp,rblp2,ralpprblp
    logical mask_nb(topo%nb(20,B))

!> proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab
!> lone-pair out-of-line damping
    real(wp) :: hblpcut

    integer i,j,ij,lina,nbb
    lina(i,j) = min(i,j)+max(i,j)*(max(i,j)-1)/2

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    energy = 0
    vector = 0
    lp_dist = 0.50-0.018*param%repz(at(B))
    hblpcut = 56

    call hbonds(A,B,ca,cb,param,topo)

    nbb = topo%nb(20,B)
!> Neighbours of B
    do i = 1,nbb
!> compute distances
      dranb(1:3,i) = xyz(1:3,A)-xyz(1:3,topo%nb(i,B))
      drbnb(1:3,i) = xyz(1:3,B)-xyz(1:3,topo%nb(i,B))
!> A-nb(B) distance
      ranb2(i) = sum(dranb(1:3,i)**2)
      ranb(i) = sqrt(ranb2(i))
!> B-nb(B) distance
      rbnb2(i) = sum(drbnb(1:3,i)**2)
      rbnb(i) = sqrt(rbnb2(i))
    end do

!> Neighbours of B
    do i = 1,nbb
      drnb(1:3,i) = xyz(1:3,topo%nb(i,B))-xyz(1:3,B)
      vector = vector+drnb(1:3,i)
    end do

    vnorm = norm2(vector)
!> lonepair coordinates
    if (vnorm .gt. 1.d-10) then
      lp = xyz(1:3,B)-lp_dist*(vector/vnorm)
    else
      lp = xyz(1:3,B)
      nbb = 0
    end if

!> A-B distance
    ij = lina(A,B)
    rab2 = sqrab(ij)
    rab = srab(ij)

!> A-H distance
    ij = lina(A,H)
    rah2 = sqrab(ij)
    rah = srab(ij)

!> B-H distance
    ij = lina(B,H)
    rbh2 = sqrab(ij)
    rbh = srab(ij)

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!> out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!> out-of-line damp: A...LP-B
    rblp2 = sum((xyz(1:3,B)-lp(1:3))**2)
    rblp = sqrt(rblp2)
    ralp2 = sum((xyz(1:3,A)-lp(1:3))**2)
    ralp = sqrt(ralp2)
    ralpprblp = ralp+rblp+1.d-12
    expo_lp = (hblpcut/radab)*(ralpprblp/rab-1.d0)
    ratio2_lp = exp(expo_lp)
    outl_lp = 2.d0/(1.d0+ratio2_lp)

!> out-of-line damp: A...nb(B)-B
    do i = 1,nbb
      ranbprbnb(i) = ranb(i)+rbnb(i)+1.d-12
      expo_nb(i) = (param%hbnbcut/radab)*(ranbprbnb(i)/rab-1.d0)
      ratio2_nb(i) = exp(-expo_nb(i))**(1.0)
      outl_nb(i) = (2.d0/(1.d0+ratio2_nb(i)))-1.0d0
    end do
    outl_nb_tot = product(outl_nb)

!> long range damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!> short range damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

!> hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!> hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!> hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot*outl_lp

!> constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!> energy
    energy = -rdamp*qhoutl*const

!> gradient
    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)
    drab(1:3) = xyz(1:3,A)-xyz(1:3,B)
    dralp(1:3) = xyz(1:3,A)-lp(1:3)
    drblp(1:3) = xyz(1:3,B)-lp(1:3)

    aterm = -rdamp*qh*outl_nb_tot*outl_lp*const
    nbterm = -rdamp*qh*outl*outl_lp*const
    lpterm = -rdamp*qh*outl*outl_nb_tot*const
    dterm = -qhoutl*const

!------------------------------------------------------------------------------
!> damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!> damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!> angular A-H...B term
!------------------------------------------------------------------------------
!> out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!> out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!------------------------------------------------------------------------------
!> angular A...LP-B term
!------------------------------------------------------------------------------
!> out of line term: rab
    tmp3 = -2.d0*lpterm*ratio2_lp*expo_lp/(1+ratio2_lp)**2/(ralpprblp-rab)
    gi = -tmp3*ralpprblp/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!> out of line term: ralp,rblp
    gi = tmp3/ralp
    dga(1:3) = gi*dralp(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp3/(rblp+1.0d-12)
    dgb(1:3) = gi*drblp(1:3)
    gb(1:3) = gb(1:3)-dga(1:3)
    glp(1:3) = -dga(1:3)!-dgb(1:3)

!> neighbor part: LP
    unit_vec = 0
    do i = 1,3
      unit_vec(i) = -1
      gii(1:3,i) = -lp_dist*dble(nbb)*(unit_vec/vnorm+(vector*vector(i)/sum(vector**2)**(1.5d0)))
      unit_vec = 0
    end do
    gnb_lp = matmul(gii,glp)

!------------------------------------------------------------------------------
!> angular A...nb(B)-B term
!------------------------------------------------------------------------------
!> out of line term: rab
    mask_nb = .true.
    do i = 1,nbb
      mask_nb(i) = .false.
      tmp2(i) = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
               & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
      gi_nb(i) = -tmp2(i)*ranbprbnb(i)/rab2
      dg(1:3) = gi_nb(i)*drab(1:3)
      ga(1:3) = ga(1:3)+dg(1:3)
      gb(1:3) = gb(1:3)-dg(1:3)
      mask_nb = .true.
    end do

!> out of line term: ranb,rbnb
    do i = 1,nbb
      gi_nb(i) = tmp2(i)/ranb(i)
      dga(1:3) = gi_nb(i)*dranb(1:3,i)
      ga(1:3) = ga(1:3)+dga(1:3)
      gi_nb(i) = tmp2(i)/rbnb(i)
      dgb(1:3) = gi_nb(i)*drbnb(1:3,i)
      gb(1:3) = gb(1:3)+dgb(1:3)
      dgnb(1:3) = -dga(1:3)-dgb(1:3)
      gnb(1:3,i) = dgnb(1:3)
    end do

!------------------------------------------------------------------------------
    if (nbb .lt. 1) then
      gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
      gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
      gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
      return
    end if

!------------------------------------------------------------------------------
!> move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)+gnb_lp(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    do i = 1,nbb
      gdr(1:3,topo%nb(i,B)) = gdr(1:3,topo%nb(i,B))+gnb(1:3,i)-gnb_lp(1:3)/dble(nbb)
    end do

  end subroutine abhgfnff_eg2_rnr

!> subroutine for case 3: A-H...B, B is 0=C including two in plane LPs at B
!> this is the multiplicative version of incorporationg etors and ebend
!> equal to abhgfnff_eg2_new multiplied by etors and eangl
  subroutine abhgfnff_eg3(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer A,B,H,n,at(n)
    real(wp) :: xyz(3,n),energy,gdr(3,n)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp) :: srab(n*(n+1)/2)    ! dist

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ddamp,rabdamp,rbhdamp
    real(wp) :: ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
    real(wp) :: xm,ym,zm
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) :: ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
    real(wp) :: drah(3),drbh(3),drab(3),drm(3)
    real(wp) :: dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) :: ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
    real(wp) :: phi,phi0,r0,t0,fc,tshift,bshift
    real(wp) :: eangl,etors,gangl(3,n),gtors(3,n)
    real(wp) :: etmp(20),g3tmp(3,3),g4tmp(3,4,20)
    real(wp) :: ratio,qhoutl,radab
    real(wp) :: gi,gi_nb(topo%nb(20,B))
    real(wp) :: tmp1,tmp2(topo%nb(20,B))
    real(wp) :: rahprbh,ranbprbnb(topo%nb(20,B))
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
    real(wp) :: eabh
    real(wp) :: aterm,dterm,nbterm,bterm,tterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: gqa,gqb,gqh
    real(wp) :: shortcut
    real(wp) :: tlist(5,topo%nb(20,topo%nb(1,B)))
    real(wp) :: vtors(2,topo%nb(20,topo%nb(1,B)))
    real(wp) :: const
    real(wp) :: outl_nb(topo%nb(20,B)),outl_nb_tot
    logical mask_nb(topo%nb(20,B)),t_mask(20)

!> proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer C,D
    integer i,j,ii,jj,kk,ll,ij,lina
    integer nbb,nbc
    integer ntors,rn

    lina(i,j) = min(i,j)+max(i,j)*(max(i,j)-1)/2

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    energy = 0
    etors = 0
    gtors = 0
    eangl = 0
    gangl = 0

    call hbonds(A,B,ca,cb,param,topo)

!> Determine all neighbors for torsion term
!>   A
!>    \         tors:
!>     H        ll
!>      :        \
!>       O        jj
!>       ||       |
!>       C        kk
!>      / \       \
!>     R1  R2      ii
!------------------------------------------
    nbb = topo%nb(20,B)
    C = topo%nb(nbb,B)
    nbc = topo%nb(20,C)
    ntors = nbc-nbb

    nbb = topo%nb(20,B)
!> Neighbours of B
    do i = 1,nbb
!> compute distances
      dranb(1:3,i) = xyz(1:3,A)-xyz(1:3,topo%nb(i,B))
      drbnb(1:3,i) = xyz(1:3,B)-xyz(1:3,topo%nb(i,B))
!> A-nb(B) distance
      ranb2(i) = sum(dranb(1:3,i)**2)
      ranb(i) = sqrt(ranb2(i))
!> B-nb(B) distance
      rbnb2(i) = sum(drbnb(1:3,i)**2)
      rbnb(i) = sqrt(rbnb2(i))
    end do

!> A-B distance
    ij = lina(A,B)
    rab2 = sqrab(ij)
    rab = srab(ij)

!> A-H distance
    ij = lina(A,H)
    rah2 = sqrab(ij)
    rah = srab(ij)

!> B-H distance
    ij = lina(B,H)
    rbh2 = sqrab(ij)
    rbh = srab(ij)

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!> out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!> out-of-line damp: A...nb(B)-B
    do i = 1,nbb
      ranbprbnb(i) = ranb(i)+rbnb(i)+1.d-12
      expo_nb(i) = (param%hbnbcut/radab)*(ranbprbnb(i)/rab-1.d0)
      ratio2_nb(i) = exp(-expo_nb(i))
      outl_nb(i) = (2.d0/(1.d0+ratio2_nb(i)))-1.0d0
    end do
    outl_nb_tot = product(outl_nb)

!> long range damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!> short range damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**6
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

!> Set up torsion paramter
    j = 0
    do i = 1,nbc
      if (topo%nb(i,C) == B) cycle
      j = j+1
      tlist(1,j) = topo%nb(i,C)
      tlist(2,j) = B
      tlist(3,j) = C
      tlist(4,j) = H
      tlist(5,j) = 2
      t0 = 180
      phi0 = t0*pi/180.
      vtors(1,j) = phi0-(pi/2.0)
      vtors(2,j) = param%tors_hb
    end do

!> Calculate etors
    do i = 1,ntors
      ii = tlist(1,i)
      jj = tlist(2,i)
      kk = tlist(3,i)
      ll = tlist(4,i)
      rn = tlist(5,i)
      phi0 = vtors(1,i)
      tshift = vtors(2,i)
      phi = valijklff(n,xyz,ii,jj,kk,ll)
      call egtors_nci_mul(ii,jj,kk,ll,rn,phi0,tshift,n,at,xyz,etmp(i),g4tmp(:,:,i))
    end do
    etors = product(etmp(1:ntors))

!> Calculate gtors
    t_mask = .true.
    do i = 1,ntors
      t_mask(i) = .false.
      ii = tlist(1,i)
      jj = tlist(2,i)
      kk = tlist(3,i)
      ll = tlist(4,i)
      gtors(1:3,ii) = gtors(1:3,ii)+g4tmp(1:3,1,i)*product(etmp(1:ntors),t_mask(1:ntors))
      gtors(1:3,jj) = gtors(1:3,jj)+g4tmp(1:3,2,i)*product(etmp(1:ntors),t_mask(1:ntors))
      gtors(1:3,kk) = gtors(1:3,kk)+g4tmp(1:3,3,i)*product(etmp(1:ntors),t_mask(1:ntors))
      gtors(1:3,ll) = gtors(1:3,ll)+g4tmp(1:3,4,i)*product(etmp(1:ntors),t_mask(1:ntors))
      t_mask = .true.
    end do

!> Calculate eangl + gangl
    r0 = 120
    phi0 = r0*pi/180.
    bshift = param%bend_hb
    fc = 1.0d0-bshift
    call bangl(xyz,kk,jj,ll,phi)
    call egbend_nci_mul(jj,kk,ll,phi0,fc,n,at,xyz,eangl,g3tmp)
    gangl(1:3,jj) = gangl(1:3,jj)+g3tmp(1:3,1)
    gangl(1:3,kk) = gangl(1:3,kk)+g3tmp(1:3,2)
    gangl(1:3,ll) = gangl(1:3,ll)+g3tmp(1:3,3)

!> hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!> hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!> hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

    qhoutl = qh*outl*outl_nb_tot

!> constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_coh

!> energy
    energy = -rdamp*qhoutl*eangl*etors*const

!> gradient
    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)
    drab(1:3) = xyz(1:3,A)-xyz(1:3,B)

    aterm = -rdamp*qh*outl_nb_tot*eangl*etors*const
    nbterm = -rdamp*qh*outl*eangl*etors*const
    dterm = -qhoutl*eangl*etors*const
    tterm = -rdamp*qhoutl*eangl*const
    bterm = -rdamp*qhoutl*etors*const

!------------------------------------------------------------------------------
!> damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!> damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!> angular A-H...B term
!------------------------------------------------------------------------------
!> out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!> out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!------------------------------------------------------------------------------
!> angular A...nb(B)-B term
!------------------------------------------------------------------------------
!> out of line term: rab
    mask_nb = .true.
    do i = 1,nbb
      mask_nb(i) = .false.
      tmp2(i) = 2.d0*nbterm*product(outl_nb,mask_nb)*ratio2_nb(i)*expo_nb(i)/&
               & (1+ratio2_nb(i))**2/(ranbprbnb(i)-rab)
      gi_nb(i) = -tmp2(i)*ranbprbnb(i)/rab2
      dg(1:3) = gi_nb(i)*drab(1:3)
      ga(1:3) = ga(1:3)+dg(1:3)
      gb(1:3) = gb(1:3)-dg(1:3)
      mask_nb = .true.
    end do

!> out of line term: ranb,rbnb
    do i = 1,nbb
      gi_nb(i) = tmp2(i)/ranb(i)
      dga(1:3) = gi_nb(i)*dranb(1:3,i)
      ga(1:3) = ga(1:3)+dga(1:3)
      gi_nb(i) = tmp2(i)/rbnb(i)
      dgb(1:3) = gi_nb(i)*drbnb(1:3,i)
      gb(1:3) = gb(1:3)+dgb(1:3)
      dgnb(1:3) = -dga(1:3)-dgb(1:3)
      gnb(1:3,i) = dgnb(1:3)
    end do

!------------------------------------------------------------------------------
    if (nbb .lt. 1) then
      gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
      gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
      gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
      return
    end if

!------------------------------------------------------------------------------
!> torsion term H...B=C<R1,R2
!------------------------------------------------------------------------------
    do i = 1,ntors
      ii = tlist(1,i)
      gdr(1:3,ii) = gdr(1:3,ii)+gtors(1:3,ii)*tterm
    end do
    gdr(1:3,jj) = gdr(1:3,jj)+gtors(1:3,jj)*tterm
    gdr(1:3,kk) = gdr(1:3,kk)+gtors(1:3,kk)*tterm
    gdr(1:3,ll) = gdr(1:3,ll)+gtors(1:3,ll)*tterm

!------------------------------------------------------------------------------
!> angle term H...B=C
!------------------------------------------------------------------------------
    gdr(1:3,jj) = gdr(1:3,jj)+gangl(1:3,jj)*bterm
    gdr(1:3,kk) = gdr(1:3,kk)+gangl(1:3,kk)*bterm
    gdr(1:3,ll) = gdr(1:3,ll)+gangl(1:3,ll)*bterm

!------------------------------------------------------------------------------
!> move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)
    do i = 1,nbb
      gdr(1:3,topo%nb(i,B)) = gdr(1:3,topo%nb(i,B))+gnb(1:3,i)
    end do

  end subroutine abhgfnff_eg3

!> subroutine for case 3: A-H...B, B is 0=C including two in plane LPs at B
!> this is the multiplicative version of incorporationg etors and ebend without neighbor LP
  subroutine abhgfnff_eg3_mul(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer A,B,H,C,D,n,at(n)
    real(wp) :: xyz(3,n),energy,gdr(3,n)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp) :: srab(n*(n+1)/2)    ! dist

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ddamp,rabdamp,rbhdamp
    real(wp) :: ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
    real(wp) :: xm,ym,zm
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) :: ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
    real(wp) :: drah(3),drbh(3),drab(3),drm(3)
    real(wp) :: dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) :: ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
    real(wp) :: phi,phi0,r0,fc,tshift,bshift
    real(wp) :: eangl,etors,gangl(3,n),gtors(3,n)
    real(wp) :: etmp,g3tmp(3,3),g4tmp(3,4)
    real(wp) :: denom,ratio,qhoutl,radab
    real(wp) :: gi,gi_nb(topo%nb(20,B))
    real(wp) :: tmp1,tmp2(topo%nb(20,B))
    real(wp) :: rahprbh,ranbprbnb(topo%nb(20,B))
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
    real(wp) :: eabh
    real(wp) :: aterm,dterm,bterm,tterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: gqa,gqb,gqh
    real(wp) :: shortcut
    real(wp) :: const
    real(wp) :: tlist(5,topo%nb(20,topo%nb(1,B)))
    real(wp) :: vtors(2,topo%nb(20,topo%nb(1,B)))
    logical mask_nb(topo%nb(20,B))

!> proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer i,j,ii,jj,kk,ll,ij,lina
    integer nbb,nbc
    integer ntors,rn

    lina(i,j) = min(i,j)+max(i,j)*(max(i,j)-1)/2

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0
    etors = 0
    gtors = 0
    eangl = 0
    gangl = 0

    call hbonds(A,B,ca,cb,param,topo)

!> Determine all neighbors for torsion term
!>   A
!>    \         tors:
!>     H        ll
!>      :        \
!>       O        jj
!>       ||       |
!>       C        kk
!>      / \       \
!>     R1  R2      ii
!------------------------------------------
    nbb = topo%nb(20,B)
    C = topo%nb(nbb,B)
    nbc = topo%nb(20,C)
    ntors = nbc-nbb

!> A-B distance
    ij = lina(A,B)
    rab2 = sqrab(ij)
    rab = srab(ij)

!> A-H distance
    ij = lina(A,H)
    rah2 = sqrab(ij)
    rah = srab(ij)

!> B-H distance
    ij = lina(B,H)
    rbh2 = sqrab(ij)
    rbh = srab(ij)

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!> out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!> long range damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!> short range damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

!> Set up torsion paramter
    j = 0
    do i = 1,nbc
      if (topo%nb(i,C) == B) cycle
      j = j+1
      tlist(1,j) = topo%nb(i,C)
      tlist(2,j) = B
      tlist(3,j) = C
      tlist(4,j) = H
      tlist(5,j) = 2
      vtors(1,j) = pi/2.0
      vtors(2,j) = 0.70
    end do

!> Calculate etors
    do i = 1,ntors
      ii = tlist(1,i)
      jj = tlist(2,i)
      kk = tlist(3,i)
      ll = tlist(4,i)
      rn = tlist(5,i)
      phi0 = vtors(1,i)
      tshift = vtors(2,i)
      phi = valijklff(n,xyz,ii,jj,kk,ll)
      call egtors_nci_mul(ii,jj,kk,ll,rn,phi0,tshift,n,at,xyz,etmp,g4tmp)
      gtors(1:3,ii) = gtors(1:3,ii)+g4tmp(1:3,1)
      gtors(1:3,jj) = gtors(1:3,jj)+g4tmp(1:3,2)
      gtors(1:3,kk) = gtors(1:3,kk)+g4tmp(1:3,3)
      gtors(1:3,ll) = gtors(1:3,ll)+g4tmp(1:3,4)
      etors = etors+etmp
    end do
    etors = etors/ntors

    r0 = 120
    phi0 = r0*pi/180.
    bshift = 0.1
    fc = 1.0d0-bshift
    call bangl(xyz,kk,jj,ll,phi)
    call egbend_nci_mul(jj,kk,ll,phi0,fc,n,at,xyz,etmp,g3tmp)
    gangl(1:3,jj) = gangl(1:3,jj)+g3tmp(1:3,1)
    gangl(1:3,kk) = gangl(1:3,kk)+g3tmp(1:3,2)
    gangl(1:3,ll) = gangl(1:3,ll)+g3tmp(1:3,3)
    eangl = eangl+etmp

!> hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!> hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!> hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

!> max distance to neighbors excluded, would lead to linear C=O-H
    qhoutl = qh*outl

!> constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!> energy
    energy = -rdamp*qhoutl*const*eangl*etors

!> gradient
    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)
    drab(1:3) = xyz(1:3,A)-xyz(1:3,B)

    aterm = -rdamp*qh*etors*eangl*const
    dterm = -qhoutl*etors*eangl*const
    tterm = -rdamp*qhoutl*eangl*const/ntors
    bterm = -rdamp*qhoutl*etors*const

!------------------------------------------------------------------------------
!> damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!>  damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!> angular A-H...B term
!------------------------------------------------------------------------------
!> out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!> out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!------------------------------------------------------------------------------
!> torsion term H...B=C<R1,R2
!------------------------------------------------------------------------------
    do i = 1,ntors
      ii = tlist(1,i)
      gdr(1:3,ii) = gdr(1:3,ii)+gtors(1:3,ii)*tterm
    end do
    gdr(1:3,jj) = gdr(1:3,jj)+gtors(1:3,jj)*tterm
    gdr(1:3,kk) = gdr(1:3,kk)+gtors(1:3,kk)*tterm
    gdr(1:3,ll) = gdr(1:3,ll)+gtors(1:3,ll)*tterm

!------------------------------------------------------------------------------
!> angle term H...B=C
!------------------------------------------------------------------------------
    gdr(1:3,jj) = gdr(1:3,jj)+gangl(1:3,jj)*bterm
    gdr(1:3,kk) = gdr(1:3,kk)+gangl(1:3,kk)*bterm
    gdr(1:3,ll) = gdr(1:3,ll)+gangl(1:3,ll)*bterm

!------------------------------------------------------------------------------
!> move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)

  end subroutine abhgfnff_eg3_mul

!> subroutine for case 3: A-H...B, B is 0=C including two in plane LPs at B
!> this is the additive version of incorporationg etors and ebend
!> This subroutine is currently unused
  subroutine abhgfnff_eg3_add(n,A,B,H,at,xyz,q,sqrab,srab,energy,gdr,param,topo)
    implicit none
    type(TGFFData),intent(in) :: param
    type(TGFFTopology),intent(in) :: topo
    integer A,B,H,C,D,n,at(n)
    real(wp) :: xyz(3,n),energy,gdr(3,n)
    real(wp) :: q(n)
    real(wp) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp) :: srab(n*(n+1)/2)    ! dist

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ddamp,rabdamp,rbhdamp
    real(wp) :: ratio1,ratio2,ratio2_nb(topo%nb(20,B)),ratio3
    real(wp) :: xm,ym,zm
    real(wp) :: rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
    real(wp) :: ranb(topo%nb(20,B)),ranb2(topo%nb(20,B)),rbnb(topo%nb(20,B)),rbnb2(topo%nb(20,B))
    real(wp) :: drah(3),drbh(3),drab(3),drm(3)
    real(wp) :: dranb(3,topo%nb(20,B)),drbnb(3,topo%nb(20,B))
    real(wp) :: dg(3),dga(3),dgb(3),dgh(3),dgnb(3)
    real(wp) :: ga(3),gb(3),gh(3),gnb(3,topo%nb(20,B))
    real(wp) :: phi,phi0,r0,fc
    real(wp) :: eangl,etors
    real(wp) :: etmp,g3tmp(3,3),g4tmp(3,4)
    real(wp) :: denom,ratio,qhoutl,radab
    real(wp) :: gi,gi_nb(topo%nb(20,B))
    real(wp) :: tmp1,tmp2(topo%nb(20,B))
    real(wp) :: rahprbh,ranbprbnb(topo%nb(20,B))
    real(wp) :: ex1a,ex2a,ex1b,ex2b,ex1h,ex2h,expo,expo_nb(topo%nb(20,B))
    real(wp) :: eabh
    real(wp) :: aterm,dterm,nbterm
    real(wp) :: qa,qb,qh
    real(wp) :: ca(2),cb(2)
    real(wp) :: gqa,gqb,gqh
    real(wp) :: shortcut
    real(wp) :: const
    real(wp) :: outl_nb(topo%nb(20,B)),outl_nb_tot
    real(wp) :: tlist(5,topo%nb(20,topo%nb(1,B)))
    real(wp) :: vtors(2,topo%nb(20,topo%nb(1,B)))
    logical mask_nb(topo%nb(20,B))

!> proportion between Rbh und Rab distance dependencies
    real(wp) :: p_bh
    real(wp) :: p_ab

    integer i,j,ii,jj,kk,ll,ij,lina
    integer nbb,nbc
    integer ntors,rn

    lina(i,j) = min(i,j)+max(i,j)*(max(i,j)-1)/2

    p_bh = 1.d0+param%hbabmix
    p_ab = -param%hbabmix

    gdr = 0
    energy = 0
    etors = 0
    eangl = 0

    call hbonds(A,B,ca,cb,param,topo)

!> Determine all neighbors for torsion term
!>   A
!>    \         tors:
!>     H        ll
!>      :        \
!>       O        jj
!>       ||       |
!>       C        kk
!>      / \       \
!>     R1  R2      ii
!------------------------------------------
    nbb = topo%nb(20,B)
    C = topo%nb(nbb,B)
    nbc = topo%nb(20,C)
    ntors = nbc-nbb

!> A-B distance
    ij = lina(A,B)
    rab2 = sqrab(ij)
    rab = srab(ij)

!> A-H distance
    ij = lina(A,H)
    rah2 = sqrab(ij)
    rah = srab(ij)

!> B-H distance
    ij = lina(B,H)
    rbh2 = sqrab(ij)
    rbh = srab(ij)

    rahprbh = rah+rbh+1.d-12
    radab = param%rad(at(A))+param%rad(at(B))

!> out-of-line damp: A-H...B
    expo = (param%hbacut/radab)*(rahprbh/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!> long range damping
    ratio1 = (rab2/param%hblongcut)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!> short range damping
    shortcut = param%hbscut*radab
    ratio3 = (shortcut/rab2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    ddamp = (-2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3/(1.d0+ratio3))
    rbhdamp = damp*((p_bh/rbh2/rbh))
    rabdamp = damp*((p_ab/rab2/rab))
    rdamp = rbhdamp+rabdamp

!> Set up torsion paramter
    j = 0
    do i = 1,nbc
      if (topo%nb(i,C) == B) cycle
      j = j+1
      tlist(1,j) = topo%nb(i,C)
      tlist(2,j) = B
      tlist(3,j) = C
      tlist(4,j) = H
      tlist(5,j) = 2
      vtors(1,j) = pi
      vtors(2,j) = 0.30
    end do

!> Calculate etors
    do i = 1,ntors
      ii = tlist(1,i)
      jj = tlist(2,i)
      kk = tlist(3,i)
      ll = tlist(4,i)
      rn = tlist(5,i)
      phi0 = vtors(1,i)
      fc = vtors(2,i)
      phi = valijklff(n,xyz,ii,jj,kk,ll)
      call egtors_nci(ii,jj,kk,ll,rn,phi0,fc,n,at,xyz,etmp,g4tmp,param)
      gdr(1:3,ii) = gdr(1:3,ii)+g4tmp(1:3,1)
      gdr(1:3,jj) = gdr(1:3,jj)+g4tmp(1:3,2)
      gdr(1:3,kk) = gdr(1:3,kk)+g4tmp(1:3,3)
      gdr(1:3,ll) = gdr(1:3,ll)+g4tmp(1:3,4)
      etors = etors+etmp
    end do

!> Calculate eangl + gangl
    write (stdout,*) 'angle atoms          phi0   phi      FC'
    r0 = 120
    phi0 = r0*pi/180.
    fc = 0.20
    call bangl(xyz,kk,jj,ll,phi)
    write (stdout,'(3i5,2x,3f8.3)') &
    &   jj,kk,ll,phi0*180./pi,phi*180./pi,fc
    call egbend_nci(jj,kk,ll,phi0,fc,n,at,xyz,etmp,g3tmp,param)
    gdr(1:3,jj) = gdr(1:3,jj)+g3tmp(1:3,1)
    gdr(1:3,kk) = gdr(1:3,kk)+g3tmp(1:3,2)
    gdr(1:3,ll) = gdr(1:3,ll)+g3tmp(1:3,3)
    eangl = eangl+etmp

!> hydrogen charge scaled term
    ex1h = exp(param%hbst*q(H))
    ex2h = ex1h+param%hbsf
    qh = ex1h/ex2h

!> hydrogen charge scaled term
    ex1a = exp(-param%hbst*q(A))
    ex2a = ex1a+param%hbsf
    qa = ex1a/ex2a

!> hydrogen charge scaled term
    ex1b = exp(-param%hbst*q(B))
    ex2b = ex1b+param%hbsf
    qb = ex1b/ex2b

!> max distance to neighbors excluded, would lead to linear C=O-H
    qhoutl = qh*outl

!> constant values, no gradient
    const = ca(2)*qa*cb(1)*qb*param%xhaci_globabh

!> energy
    energy = -rdamp*qhoutl*const+etors+eangl

!> gradient
    drah(1:3) = xyz(1:3,A)-xyz(1:3,H)
    drbh(1:3) = xyz(1:3,B)-xyz(1:3,H)
    drab(1:3) = xyz(1:3,A)-xyz(1:3,B)

    aterm = -rdamp*qh*const
    dterm = -qhoutl*const

!------------------------------------------------------------------------------
!> damping part: rab
    gi = ((rabdamp+rbhdamp)*ddamp-3.d0*rabdamp)/rab2
    gi = gi*dterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = dg(1:3)
    gb(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!> damping part: rbh
    gi = -3.d0*rbhdamp/rbh2
    gi = gi*dterm
    dg(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dg(1:3)
    gh(1:3) = -dg(1:3)

!------------------------------------------------------------------------------
!> angular A-H...B term
!------------------------------------------------------------------------------
!> out of line term: rab
    tmp1 = -2.d0*aterm*ratio2*expo/(1+ratio2)**2/(rahprbh-rab)
    gi = -tmp1*rahprbh/rab2
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = ga(1:3)+dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!> out of line term: rah,rbh
    gi = tmp1/rah
    dga(1:3) = gi*drah(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = tmp1/rbh
    dgb(1:3) = gi*drbh(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgh(1:3) = -dga(1:3)-dgb(1:3)
    gh(1:3) = gh(1:3)+dgh(1:3)

!------------------------------------------------------------------------------
!> move gradients into place
    gdr(1:3,A) = gdr(1:3,A)+ga(1:3)
    gdr(1:3,B) = gdr(1:3,B)+gb(1:3)
    gdr(1:3,H) = gdr(1:3,H)+gh(1:3)

  end subroutine abhgfnff_eg3_add

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! XB energy and analytical gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine rbxgfnff_eg(n,A,B,X,at,xyz,q,energy,gdr,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer               :: A,B,X,n,at(n)
    real(wp)                :: xyz(3,n)
    real(wp),intent(inout)  :: energy,gdr(3,3)
    real(wp)                :: q(n)

    real(wp) :: outl,dampl,damps,rdamp,damp
    real(wp) :: ratio1,ratio2,ratio3
    real(wp) :: rab,rax,rbx,rab2,rax2,rbx2,rax4,rbx4
    real(wp) :: drax(3),drbx(3),drab(3),drm(3)
    real(wp) :: dg(3),dga(3),dgb(3),dgx(3)
    real(wp) :: gi,ga(3),gb(3),gx(3)
    real(wp) :: ex1_a,ex2_a,ex1_b,ex2_b,ex1_x,ex2_x,expo
    real(wp) :: aterm,dterm
    real(wp) :: qa,qb,qx
    real(wp) :: cx,cb
    real(wp) :: gqa,gqb,gqx
    real(wp) :: shortcut,const

    integer i,j

    gdr = 0
    energy = 0

    cb = 1.!param%xhbas(at(B))
    cx = param%xbaci(at(X))

!> compute distances
    drax(1:3) = xyz(1:3,A)-xyz(1:3,X)
    drbx(1:3) = xyz(1:3,B)-xyz(1:3,X)
    drab(1:3) = xyz(1:3,A)-xyz(1:3,B)

!> A-B distance
    rab2 = sum(drab**2)
    rab = sqrt(rab2)

!> A-X distance
    rax2 = sum(drax**2)
    rax = sqrt(rax2)+1.d-12

!> B-X distance
    rbx2 = sum(drbx**2)
    rbx = sqrt(rbx2)+1.d-12

!> out-of-line damp
    expo = param%xbacut*((rax+rbx)/rab-1.d0)
    if (expo .gt. 15.0d0) return ! avoid overflow
    ratio2 = exp(expo)
    outl = 2.d0/(1.d0+ratio2)

!> long damping
    ratio1 = (rbx2/param%hblongcut_xb)**param%hbalp
    dampl = 1.d0/(1.d0+ratio1)

!> short damping
    shortcut = param%xbscut*(param%rad(at(A))+param%rad(at(B)))
    ratio3 = (shortcut/rbx2)**param%hbalp
    damps = 1.d0/(1.d0+ratio3)

    damp = damps*dampl
    rdamp = damp/rbx2/rbx ! **2

!> halogen charge scaled term
    ex1_x = exp(param%xbst*q(X))
    ex2_x = ex1_x+param%xbsf
    qx = ex1_x/ex2_x

!> donor charge scaled term
    ex1_b = exp(-param%xbst*q(B))
    ex2_b = ex1_b+param%xbsf
    qb = ex1_b/ex2_b

!> constant values, no gradient
    const = cb*qb*cx*qx

!> r^3 only sligxtly better than r^4
    aterm = -rdamp*const
    dterm = -outl*const
    energy = -rdamp*outl*const

!> damping part rab
    gi = rdamp*(-(2.d0*param%hbalp*ratio1/(1.d0+ratio1))+(2.d0*param%hbalp*ratio3&
   &     /(1.d0+ratio3))-3.d0)/rbx2   ! 4,5,6 instead of 3.
    gi = gi*dterm
    dg(1:3) = gi*drbx(1:3)
    gb(1:3) = dg(1:3)
    gx(1:3) = -dg(1:3)

!> out of line term: rab
    gi = 2.d0*ratio2*expo*(rax+rbx)/(1.d0+ratio2)**2/(rax+rbx-rab)/rab2
    gi = gi*aterm
    dg(1:3) = gi*drab(1:3)
    ga(1:3) = +dg(1:3)
    gb(1:3) = gb(1:3)-dg(1:3)

!> out of line term: rax,rbx
    gi = -2.d0*ratio2*expo/(1.d0+ratio2)**2/(rax+rbx-rab)/rax
    gi = gi*aterm
    dga(1:3) = gi*drax(1:3)
    ga(1:3) = ga(1:3)+dga(1:3)
    gi = -2.d0*ratio2*expo/(1.d0+ratio2)**2/(rax+rbx-rab)/rbx
    gi = gi*aterm
    dgb(1:3) = gi*drbx(1:3)
    gb(1:3) = gb(1:3)+dgb(1:3)
    dgx(1:3) = -dga(1:3)-dgb(1:3)
    gx(1:3) = gx(1:3)+dgx(1:3)

!> move gradients into place
    gdr(1:3,1) = ga(1:3)
    gdr(1:3,2) = gb(1:3)
    gdr(1:3,3) = gx(1:3)

    return
  end subroutine rbxgfnff_eg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! taken from D3 ATM code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine batmgfnff_eg(n,iat,jat,kat,at,xyz,q,sqrab,srab,energy,g,param)
    implicit none
    type(TGFFData),intent(in) :: param
    integer,intent(in) :: iat,jat,kat,n,at(n)
    real(wp),intent(in) :: xyz(3,n),q(n)
    real(wp),intent(out) :: energy,g(3,3)
    real(wp),intent(in) :: sqrab(n*(n+1)/2)   ! squared dist
    real(wp),intent(in) :: srab(n*(n+1)/2)   ! dist

    real(wp) :: r2ij,r2jk,r2ik,c9,mijk,imjk,ijmk,rijk3,ang,angr9,rav3
    real(wp) :: rij(3),rik(3),rjk(3),drij,drik,drjk,dang,ff,fi,fj,fk,fqq
    parameter(fqq=3.0d0)
    integer linij,linik,linjk,lina,i,j
    lina(i,j) = min(i,j)+max(i,j)*(max(i,j)-1)/2

    fi = (1.d0-fqq*q(iat))
    fi = min(max(fi,-4.0d0),4.0d0)
    fj = (1.d0-fqq*q(jat))
    fj = min(max(fj,-4.0d0),4.0d0)
    fk = (1.d0-fqq*q(kat))
    fk = min(max(fk,-4.0d0),4.0d0)
!> charge term
    ff = fi*fj*fk
!> strength of interaction
    c9 = ff*param%zb3atm(at(iat))*param%zb3atm(at(jat))*param%zb3atm(at(kat))
    linij = lina(iat,jat)
    linik = lina(iat,kat)
    linjk = lina(jat,kat)
    r2ij = sqrab(linij)
    r2jk = sqrab(linjk)
    r2ik = sqrab(linik)
    mijk = -r2ij+r2jk+r2ik
    imjk = r2ij-r2jk+r2ik
    ijmk = r2ij+r2jk-r2ik
    rijk3 = r2ij*r2jk*r2ik
    rav3 = rijk3**1.5 ! R^9
    ang = 0.375d0*ijmk*imjk*mijk/rijk3
    angr9 = (ang+1.0d0)/rav3
    energy = c9*angr9 ! energy

!> derivatives of each part w.r.t. r_ij,jk,ik
    dang = -0.375d0*(r2ij**3+r2ij**2*(r2jk+r2ik) &
&             +r2ij*(3.0d0*r2jk**2+2.0*r2jk*r2ik+3.0*r2ik**2) &
&             -5.0*(r2jk-r2ik)**2*(r2jk+r2ik)) &
&             /(srab(linij)*rijk3*rav3)
    drij = -dang*c9
    dang = -0.375d0*(r2jk**3+r2jk**2*(r2ik+r2ij) &
&             +r2jk*(3.0d0*r2ik**2+2.0*r2ik*r2ij+3.0*r2ij**2) &
&             -5.0*(r2ik-r2ij)**2*(r2ik+r2ij)) &
&             /(srab(linjk)*rijk3*rav3)
    drjk = -dang*c9
    dang = -0.375d0*(r2ik**3+r2ik**2*(r2jk+r2ij) &
&             +r2ik*(3.0d0*r2jk**2+2.0*r2jk*r2ij+3.0*r2ij**2) &
&             -5.0*(r2jk-r2ij)**2*(r2jk+r2ij)) &
&             /(srab(linik)*rijk3*rav3)
    drik = -dang*c9

    rij = xyz(:,jat)-xyz(:,iat)
    rik = xyz(:,kat)-xyz(:,iat)
    rjk = xyz(:,kat)-xyz(:,jat)
    g(:,1) = drij*rij/srab(linij)
    g(:,1) = g(:,1)+drik*rik/srab(linik)
    g(:,2) = drjk*rjk/srab(linjk)
    g(:,2) = g(:,2)-drij*rij/srab(linij)
    g(:,3) = -drik*rik/srab(linik)
    g(:,3) = g(:,3)-drjk*rjk/srab(linjk)

  end subroutine batmgfnff_eg

!========================================================================================!
!> torsion term for rotation around triple bonded carbon
subroutine sTors_eg(m, n, xyz, topo, energy, dg)
   integer, intent(in) :: m
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3,n)
   type(TGFFTopology), intent(in) :: topo
   real(wp), intent(out) :: energy
   real(wp), intent(out) :: dg(3,n)
   integer :: c1,c2,c3,c4
   integer :: i

   !> torsion angle between C1-C4
   real(wp) :: phi
   real(wp) :: erefhalf
   real(wp) :: dp1(3),dp2(3),dp3(3),dp4(3)

   energy = 0.0_wp
   dg(:,:) = 0.0_wp

   if ( .not. any(topo%sTorsl(:,m) .eq. 0)) then

      c1 = topo%sTorsl(1,m)
      c2 = topo%sTorsl(2,m)
      c3 = topo%sTorsl(5,m)
      c4 = topo%sTorsl(6,m)

      ! dihedral angle in radians!
      phi=valijklff(n,xyz,c1,c2,c3,c4)
      call dphidr(n,xyz,c1,c2,c3,c4,phi,dp1,dp2,dp3,dp4)

      ! reference energy for torsion of 90° !
      ! calculated with DLPNO-CCSD(T) CBS on diphenylacetylene !
      erefhalf = 3.75_wp*1.0e-4_wp  ! approx 1.97 kJ/mol !
      energy = -erefhalf*cos(2.0_wp*phi) + erefhalf
      do i=1, 3
         dg(i, c1) = dg(i, c1) + erefhalf*2.0_wp*sin(2.0_wp*phi)*dp1(i)
         dg(i, c2) = dg(i, c2) + erefhalf*2.0_wp*sin(2.0_wp*phi)*dp2(i)
         dg(i, c3) = dg(i, c3) + erefhalf*2.0_wp*sin(2.0_wp*phi)*dp3(i)
         dg(i, c4) = dg(i, c4) + erefhalf*2.0_wp*sin(2.0_wp*phi)*dp4(i)
      enddo
   endif

end subroutine sTors_eg



!========================================================================================!
end module gfnff_engrad_module
