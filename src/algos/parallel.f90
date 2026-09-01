!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2023  Philipp Pracht
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

!> A collection of routines to set up OMP-parallel runs of MDs and optimizations.

!========================================================================================!
!========================================================================================!
!> Interfaces to handle optional arguments
!========================================================================================!
!========================================================================================!
module parallel_interface
!*******************************************************
!* module to load an interface to the parallel routines
!* mandatory to handle any optional input arguments
!*******************************************************
  implicit none
  interface
    subroutine crest_sploop(env,nat,nall,at,xyz,eread)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
      use optimize_module
      use iomod,only:makedir,directory_exist,remove
      implicit none
      type(systemdata),intent(inout) :: env
      real(wp),intent(inout) :: xyz(3,nat,nall)
      integer,intent(in)  :: at(nat)
      real(wp),intent(inout) :: eread(nall)
      integer,intent(in) :: nat,nall
    end subroutine crest_sploop
  end interface

  interface
    subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump,customcalc, &
    & poststage_buffer,poststage_written)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
      use crest_poststage_ensemble,only:poststage_ensemble
      use optimize_module
      use iomod,only:makedir,directory_exist,remove
      use crest_restartlog,only:trackrestart,restart_write_dummy
      implicit none
      type(systemdata),target,intent(inout) :: env
      real(wp),intent(inout) :: xyz(3,nat,nall)
      integer,intent(in)  :: at(nat)
      real(wp),intent(inout) :: eread(nall)
      integer,intent(in) :: nat,nall
      logical,intent(in) :: dump
      type(calcdata),intent(in),target,optional :: customcalc
      type(poststage_ensemble),intent(inout),optional :: poststage_buffer
      logical,intent(out),optional :: poststage_written
    end subroutine crest_oloop
  end interface

  interface
    subroutine crest_search_multimd2(env,mols,mddats,nsim,poststage_out)
      use crest_data,only:systemdata
      use strucrd,only:coord
      use dynamics_module,only:mddata
      use crest_poststage_ensemble,only:poststage_ensemble
      implicit none
      type(systemdata),intent(inout) :: env
      integer,intent(in) :: nsim
      type(coord),intent(in) :: mols(nsim)
      type(mddata),intent(inout) :: mddats(nsim)
      type(poststage_ensemble),intent(inout),optional :: poststage_out
    end subroutine crest_search_multimd2
  end interface
end module parallel_interface

!========================================================================================!
!========================================================================================!
!> Routines for concurrent singlepoint evaluations
!========================================================================================!
!========================================================================================!
subroutine crest_sploop(env,nat,nall,at,xyz,eread)
!***************************************************************
!* subroutine crest_sploop
!* This subroutine performs concurrent singlepoint evaluations
!* for the given ensemble. Input eread is overwritten
!***************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
  use optimize_module
  use iomod,only:makedir,directory_exist,remove
  implicit none
  type(systemdata),intent(inout) :: env
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(in)  :: at(nat)
  real(wp),intent(inout) :: eread(nall)
  integer,intent(in) :: nat,nall

  type(coord),allocatable :: mols(:)
  integer :: i,j,k,l,io,ich,ich2,c,z,job_id,zcopy
  logical :: pr,wr,ex
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:),grads(:,:,:)
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent,runtime

  type(timer) :: profiler
  integer :: T,Tn  !> threads and threads per core
  logical :: nested

!>--- check if we have any calculation settings allocated
  if (env%calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
  call new_ompautoset(env,'auto_nested',nall,T,Tn)
  nested = env%omp_allow_nested


!>--- prepare objects for parallelization
  T = env%threads
  allocate (calculations(T),source=env%calc)
  allocate (mols(T))
  do i = 1,T
    do j = 1,env%calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if(allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
      if(allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
    allocate (mols(i)%at(nat),mols(i)%xyz(3,nat))
  end do

!>--- printout directions and timer initialization
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log
  call profiler%init(1)
  call profiler%start(1)

!>--- first progress printout (initializes progress variables)
  call crest_oloop_pr_progress(env,nall,0)

!>--- shared variables
  allocate (grads(3,nat,T),source=0.0_wp)
  c = 0  !> counter of successfull optimizations
  k = 0  !> counter of total optimization (fail+success)
  z = 0  !> counter to perform optimization in right order (1...nall)
  eread(:) = 0.0_wp
  grads(:,:,:) = 0.0_wp
!>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nat,nall,at,xyz,eread,grads,c,k,z,pr,wr) &
  !$omp shared(ich,ich2,mols, nested,Tn)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(i,j,job,energy,io,thread_id,zcopy)
    call initsignal()

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !>--- modify calculation spaces
    !$omp critical
    z = z+1
    zcopy = z
    mols(job)%nat = nat
    mols(job)%at(:) = at(:)
    mols(job)%xyz(:,:) = xyz(:,:,z)
    !$omp end critical

    !>-- engery+gradient call
    call engrad(mols(job),calculations(job),energy,grads(:,:,job),io)

    !$omp critical
    if (io == 0) then
      !>--- successful optimization (io==0)
      c = c+1
      eread(zcopy) = energy
    else
      eread(zcopy) = 0.0_wp
    end if
    k = k+1
    !>--- print progress
    call crest_oloop_pr_progress(env,nall,k)
    !$omp end critical
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

!>--- finalize progress printout
  call crest_oloop_pr_progress(env,nall,-1)

!>--- stop timer
  call profiler%stop(1)

!>--- prepare some summary printout
  percent = float(c)/float(nall)*100.0_wp
  write (atmp,'(f5.1,a)') percent,'% success)'
  write (stdout,'(">",1x,i0,a,i0,a,a)') c,' of ',nall,' structures successfully evaluated (', &
  &     trim(adjustl(atmp))
  write (atmp,'(">",1x,a,i0,a)') 'Total runtime for ',nall,' singlepoint calculations:'
  call profiler%write_timing(stdout,1,trim(atmp),.true.)
  runtime = profiler%get(1)
  write (atmp,'(f16.3,a)') runtime/real(nall,wp),' sec'
  write (stdout,'(a,a,a)') '> Corresponding to approximately ',trim(adjustl(atmp)), &
  &                       ' per processed structure'

  deallocate (grads)
  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  return
end subroutine crest_sploop

!========================================================================================!
!========================================================================================!
!> Routines for concurrent geometry optimization
!========================================================================================!
!========================================================================================!
subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump,customcalc, &
& poststage_buffer,poststage_written)
!*******************************************************************************
!* subroutine crest_oloop
!* This subroutine performs concurrent geometry optimizations
!* for the given ensemble. Inputs xyz and eread are overwritten
!* env        - contains parallelization and other program settings
!* dump       - decides on whether to dump an ensemble file
!*              WARNING: the ensemble file will NOT be in the same order
!*              as the input xyz array. However, the overwritten xyz will be! 
!* customcalc - customized (optional) calculation level data
!* poststage_buffer - when present with dump=.true., write the exact canonical
!*              XYZ produced by an empty crest_refine pass and retain its
!*              reparsed Angstrom data in task-completion order
!* poststage_written - true only after that canonical file and buffer complete
!*
!* IMPORTANT: xyz should be in Bohr(!) for this routine
!******************************************************************************
  use crest_parameters,only:wp,stdout,sep,bohr
  use iso_fortran_env,only:int64
  use iso_c_binding,only:c_int
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
  use crest_poststage_ensemble,only:poststage_ensemble,poststage_comment_length, &
  & round_fixed_decimal_10,write_canonical_ensemble_fast
  use optimize_module
  use iomod,only:makedir,directory_exist,remove
  use crest_restartlog,only:trackrestart,restart_write_dummy
  implicit none
  type(systemdata),target,intent(inout) :: env
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(in)  :: at(nat)
  real(wp),intent(inout) :: eread(nall)
  integer,intent(in) :: nat,nall
  logical,intent(in) :: dump
  type(calcdata),intent(in),target,optional :: customcalc
  type(poststage_ensemble),intent(inout),optional :: poststage_buffer
  logical,intent(out),optional :: poststage_written

  type(coord),allocatable :: mols(:)
  type(coord),allocatable :: molsnew(:)
  type(coord) :: canonical_mol
  integer :: i,j,k,l,io,ich,ich2,c,z,job_id,zcopy,ninitialized,ngfnff
  integer :: ndump,canonical_status,canonical_iostat
  integer :: canonical_workers,canonical_team
  integer(int64) :: canonical_fallbacks,canonical_frame_fallbacks
  integer(int64) :: canonical_artifact_bytes
  logical :: pr,wr,ex,canonical_output
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grads(:,:,:)
  real(wp),allocatable :: canonical_xyz(:,:,:),canonical_eread(:)
  integer,allocatable :: canonical_frame_status(:)
  character(len=poststage_comment_length),allocatable :: canonical_comments(:)
  real(wp) :: poststage_clock
  real(wp) :: canonical_write_seconds
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  character(len=1024) :: canonical_message
  real(wp) :: percent,runtime
  type(calcdata),pointer :: mycalc
  type(calcdata),target :: local_template
  type(timer) :: profiler
  integer :: T,Tn  !> threads and threads per core
  logical :: nested
  logical :: optimizer_affinity,affinity_tasks_started
  integer :: affinity_rank,affinity_team_size
  integer :: affinity_bind_errors,affinity_restore_errors
  integer :: affinity_first_bind_error,affinity_first_restore_error
  integer(c_int) :: affinity_cstat,affinity_finish_status,affinity_gate_enabled
  integer(c_int) :: affinity_standalone_allowed

  interface
    integer(c_int) function c_optimizer_affinity_resolve(outer_threads,inner_threads, &
    & standalone_allowed,enabled_out) &
    & bind(C,name='crest_optimizer_affinity_resolve')
      import :: c_int
      integer(c_int),value :: outer_threads,inner_threads,standalone_allowed
      integer(c_int),intent(out) :: enabled_out
    end function c_optimizer_affinity_resolve

    integer(c_int) function c_optimizer_affinity_prepare(expected_threads) &
    & bind(C,name='crest_optimizer_affinity_prepare')
      import :: c_int
      integer(c_int),value :: expected_threads
    end function c_optimizer_affinity_prepare

    integer(c_int) function c_optimizer_affinity_bind(rank,team_size) &
    & bind(C,name='crest_optimizer_affinity_bind')
      import :: c_int
      integer(c_int),value :: rank,team_size
    end function c_optimizer_affinity_bind

    integer(c_int) function c_optimizer_affinity_restore(rank,team_size) &
    & bind(C,name='crest_optimizer_affinity_restore')
      import :: c_int
      integer(c_int),value :: rank,team_size
    end function c_optimizer_affinity_restore

    integer(c_int) function c_optimizer_affinity_finish() &
    & bind(C,name='crest_optimizer_affinity_finish')
      import :: c_int
    end function c_optimizer_affinity_finish
  end interface

  canonical_output = present(poststage_buffer).and.dump
  if (present(poststage_written)) poststage_written = .false.
  if (present(poststage_buffer)) call poststage_buffer%clear()

!>--- decide wether to skip this call
  if (trackrestart(env)) then
    call restart_write_dummy(ensemblefile)
    return
  end if

!>--- check which calc to use
  if(present(customcalc))then
    local_template = customcalc
  else
    local_template = env%calc
  endif
  !> Always seed a call-local template. Neither env%calc nor an INTENT(IN)
  !> custom calculator is mutated by canonical topology preparation.
  mycalc => local_template

!>--- check if we have any calculation settings allocated
  if (mycalc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

!>--- establish the outer worker count before preparing the shared topology
  call new_ompautoset(env,'auto_nested',nall,T,Tn)
  nested = env%omp_allow_nested
  affinity_gate_enabled = 0_c_int
  affinity_standalone_allowed = 0_c_int
  if (env%crestver == crest_screen) affinity_standalone_allowed = 1_c_int
  affinity_cstat = c_optimizer_affinity_resolve(int(T,c_int),int(Tn,c_int), &
  & affinity_standalone_allowed,affinity_gate_enabled)
  if (affinity_cstat /= 0_c_int) then
    write (stdout,'(1x,a,i0)') 'Optimizer affinity opt-in/state is invalid; errno=', &
    & int(affinity_cstat)
    env%iostatus_meta = status_config
    return
  end if
  optimizer_affinity = affinity_gate_enabled == 1_c_int

!>--- initialize from ensemble frame 1 before cloning the worker calculators
!  The historical worker-local setup ran inside the outer OpenMP region and
!  therefore generated every topology with one active thread.  Preserve those
!  exact parameters while doing the work only once: temporarily serialize the
!  canonical setup, then restore the selected outer worker count before copies
!  or optimization work are created.
  canonical_mol%nat = nat
  allocate (canonical_mol%at(nat),canonical_mol%xyz(3,nat))
  canonical_mol%at = at
  canonical_mol%xyz = xyz(:,:,1)
  call omp_set_num_threads(1)
  call prepare_gfnff_topology(canonical_mol,mycalc,io,ninitialized,ngfnff)
  call omp_set_num_threads(T)
  call canonical_mol%deallocate()
  if (io /= 0) then
    write (stdout,'(1x,a,i0)') 'Canonical GFN-FF topology initialization failed; iostat=',io
    env%iostatus_meta = status_failed
    return
  end if
  if (ngfnff > 0) then
    if (ninitialized > 0) then
      write (stdout,'(1x,a)') 'Canonical GFN-FF topology ready before worker cloning (initialized).'
    else
      write (stdout,'(1x,a)') 'Canonical GFN-FF topology ready before worker cloning (retained).'
    end if
  end if

!>--- arm placement only after explicit screen opt-in or a validated process batch
  if (optimizer_affinity) then
    if (omp_get_proc_bind() /= omp_proc_bind_false .or. omp_get_num_places() /= 0) then
      write (stdout,'(1x,a)') 'Optimizer affinity requires an unbound process parent.'
      env%iostatus_meta = status_config
      return
    end if
    affinity_cstat = c_optimizer_affinity_prepare(int(T,c_int))
    if (affinity_cstat /= 0_c_int) then
      write (stdout,'(1x,a,i0)') 'Optimizer affinity preparation failed; errno=', &
      & int(affinity_cstat)
      env%iostatus_meta = status_config
      return
    end if
  end if

!>--- prepare objects for parallelization
  allocate (calculations(T),source=mycalc)
  allocate (mols(T),molsnew(T))
  do i = 1,T
    do j = 1,mycalc%ncalculations
      !>--- directories and io preparation
      ex = directory_exist(mycalc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(mycalc%calcs(j)%calcspace))
      end if
      if(calculations(i)%calcs(j)%id == jobtype%tblite)then
         calculations(i)%optnewinit=.true.
      endif
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = mycalc%calcs(j)%calcspace//trim(atmp)
      if(allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
      if(allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
    allocate (mols(i)%at(nat),mols(i)%xyz(3,nat))
    allocate (molsnew(i)%at(nat),molsnew(i)%xyz(3,nat))
  end do

!>--- printout directions and timer initialization
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log
  if (dump) then
    if (canonical_output) then
      ich = -1
      open (newunit=ich2,file=ensembleelog,iostat=io)
      if (io /= 0) then
        write (stdout,'(1x,a,i0)') &
        & 'Canonical optimizer energy-log open failed; iostat=',io
        env%iostatus_meta = status_failed
        call remove(ensemblefile)
        call remove(ensembleelog)
        deallocate (calculations)
        if (allocated(mols)) deallocate (mols)
        if (allocated(molsnew)) deallocate (molsnew)
        if (optimizer_affinity) then
          affinity_finish_status = c_optimizer_affinity_finish()
          if (affinity_finish_status /= 0_c_int) write (stdout,'(1x,a,i0)') &
          & 'Optimizer affinity finish after canonical setup failure returned errno=', &
          & int(affinity_finish_status)
        end if
        return
      end if
    else
      open (newunit=ich,file=ensemblefile)
      open (newunit=ich2,file=ensembleelog)
    end if
  end if
  if (canonical_output) then
    allocate (canonical_xyz(3,nat,nall),canonical_eread(nall), &
    & canonical_comments(nall),canonical_frame_status(nall),stat=io)
    if (io /= 0) then
      write (stdout,'(1x,a,i0)') &
      & 'Canonical optimizer buffer allocation failed; stat=',io
      env%iostatus_meta = status_failed
      close (ich2,iostat=canonical_status)
      call remove(ensemblefile)
      call remove(ensembleelog)
      if (allocated(canonical_xyz)) deallocate (canonical_xyz)
      if (allocated(canonical_eread)) deallocate (canonical_eread)
      if (allocated(canonical_comments)) deallocate (canonical_comments)
      if (allocated(canonical_frame_status)) deallocate (canonical_frame_status)
      deallocate (calculations)
      if (allocated(mols)) deallocate (mols)
      if (allocated(molsnew)) deallocate (molsnew)
      if (optimizer_affinity) then
        affinity_finish_status = c_optimizer_affinity_finish()
        if (affinity_finish_status /= 0_c_int) write (stdout,'(1x,a,i0)') &
        & 'Optimizer affinity finish after canonical setup failure returned errno=', &
        & int(affinity_finish_status)
      end if
      return
    end if
    canonical_frame_status = 0
  end if
  call profiler%init(1)
  call profiler%start(1)

!>--- first progress printout (initializes progress variables)
  call crest_oloop_pr_progress(env,nall,0)

!>--- shared variables
  allocate (grads(3,nat,T),source=0.0_wp)
  c = 0  !> counter of successfull optimizations
  ndump = 0 !> successful structures written, in task-completion order
  k = 0  !> counter of total optimization (fail+success)
  z = 0  !> counter to perform optimization in right order (1...nall)
  eread(:) = 0.0_wp
  grads(:,:,:) = 0.0_wp
  affinity_bind_errors = 0
  affinity_restore_errors = 0
  affinity_first_bind_error = 0
  affinity_first_restore_error = 0
  affinity_finish_status = 0_c_int
  affinity_tasks_started = .false.
  canonical_iostat = 0
!>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nat,nall,at,xyz,eread,grads,c,k,z,pr,wr,dump) &
  !$omp shared(ich,ich2,mols,molsnew,nested,Tn,optimizer_affinity) &
  !$omp shared(canonical_output,canonical_xyz,canonical_eread,canonical_comments) &
  !$omp shared(ndump,canonical_iostat) &
  !$omp shared(affinity_bind_errors,affinity_restore_errors) &
  !$omp shared(affinity_first_bind_error,affinity_first_restore_error,affinity_tasks_started) &
  !$omp private(affinity_rank,affinity_team_size,affinity_cstat)
  affinity_rank = omp_get_thread_num()
  affinity_team_size = omp_get_num_threads()
  if (optimizer_affinity) then
    affinity_cstat = c_optimizer_affinity_bind(int(affinity_rank,c_int), &
    & int(affinity_team_size,c_int))
    if (affinity_cstat /= 0_c_int) then
      !$omp critical(crest_optimizer_affinity_status)
      affinity_bind_errors = affinity_bind_errors+1
      if (affinity_first_bind_error == 0) affinity_first_bind_error = int(affinity_cstat)
      !$omp end critical(crest_optimizer_affinity_status)
    end if
    !$omp barrier
  end if
  !$omp single
  if (.not.optimizer_affinity .or. affinity_bind_errors == 0) then
  affinity_tasks_started = .true.
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) &
    !$omp private(j,job,energy,io,atmp,gnorm,thread_id,zcopy)
    call initsignal()

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !>--- modify calculation spaces
    !$omp critical
    z = z+1
    zcopy = z
    mols(job)%nat = nat
    mols(job)%at(:) = at(:)
    mols(job)%xyz(:,:) = xyz(:,:,z)

    molsnew(job)%nat = nat
    molsnew(job)%at(:) = at(:)
    molsnew(job)%xyz(:,:) = xyz(:,:,z)
    !$omp end critical

    !>-- geometry optimization
    !> Each worker calculator persists across ensemble members.  Rebuild only
    !> its geometry-dependent GFN-FF HB/XB list at this input boundary; retain
    !> canonical topology, exact EEQ, frozen-host caches, and all workspaces.
    call request_gfnff_hbond_update(calculations(job))
    call optimize_geometry(mols(job),molsnew(job),calculations(job),energy,grads(:,:,job),pr,wr,io)

    !$omp critical
    if (io == 0) then
      !>--- successful optimization (io==0)
      c = c+1
      if (dump) then
        gnorm = norm2(grads(:,:,job))
        write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
        molsnew(job)%comment = trim(atmp)
        if (canonical_output) then
          ndump = ndump+1
          canonical_xyz(:,:,ndump) = molsnew(job)%xyz
          canonical_comments(ndump) = trim(atmp)
          !> Preserve coord%append's multiply/divide side effect on the array
          !> returned to the optimization caller, without serial formatting.
          molsnew(job)%xyz = molsnew(job)%xyz*bohr
          molsnew(job)%xyz = molsnew(job)%xyz/bohr
        else
          call molsnew(job)%append(ich)
        end if
        call calc_eprint(calculations(job),energy,calculations(job)%etmp,gnorm,ich2)
      end if
      eread(zcopy) = energy
      xyz(:,:,zcopy) = molsnew(job)%xyz(:,:)
    else if(io==calculations(job)%maxcycle .and. calculations(job)%anopt) then
      !>--- allow partial optimization?
      c = c+1
      eread(zcopy) = energy
      xyz(:,:,zcopy) = molsnew(job)%xyz(:,:)
    else
      eread(zcopy) = 1.0_wp
    end if
    k = k+1
    !>--- print progress
    call crest_oloop_pr_progress(env,nall,k)
    !$omp end critical
    !$omp end task
  end do
  !$omp taskwait
  end if
  !$omp end single
  if (optimizer_affinity) then
    affinity_cstat = c_optimizer_affinity_restore(int(affinity_rank,c_int), &
    & int(affinity_team_size,c_int))
    if (affinity_cstat /= 0_c_int) then
      !$omp critical(crest_optimizer_affinity_status)
      affinity_restore_errors = affinity_restore_errors+1
      if (affinity_first_restore_error == 0) affinity_first_restore_error = int(affinity_cstat)
      !$omp end critical(crest_optimizer_affinity_status)
    end if
  end if
  !$omp end parallel

  if (optimizer_affinity) then
    affinity_finish_status = c_optimizer_affinity_finish()
    if (affinity_bind_errors /= 0 .or. affinity_restore_errors /= 0 .or. &
    & affinity_finish_status /= 0_c_int) then
      if (affinity_bind_errors /= 0) then
        write (stdout,'(1x,a,i0,a,i0)') 'Optimizer affinity bind failed on ', &
        & affinity_bind_errors,' workers; first errno=',affinity_first_bind_error
      end if
      if (affinity_restore_errors /= 0) then
        write (stdout,'(1x,a,i0,a,i0)') 'Optimizer affinity restore failed on ', &
        & affinity_restore_errors,' workers; first errno=',affinity_first_restore_error
      end if
      if (affinity_finish_status /= 0_c_int) write (stdout,'(1x,a,i0)') &
      & 'Optimizer affinity finish failed; errno=',int(affinity_finish_status)
      if (affinity_tasks_started .and. affinity_bind_errors /= 0) then
        write (stdout,'(1x,a)') 'Optimizer affinity invariant violated: tasks were started.'
        env%iostatus_meta = status_failed
      else if (affinity_restore_errors /= 0 .or. affinity_finish_status /= 0_c_int) then
        env%iostatus_meta = status_failed
      else
        env%iostatus_meta = status_config
      end if
      if (dump) then
        if (canonical_output) then
          close (ich2,iostat=canonical_status)
        else
          close (ich)
          close (ich2)
        end if
        call remove(ensemblefile)
        call remove(ensembleelog)
      end if
      call profiler%stop(1)
      call profiler%clear()
      deallocate (grads,calculations)
      if (allocated(mols)) deallocate (mols)
      if (allocated(molsnew)) deallocate (molsnew)
      if (allocated(canonical_xyz)) deallocate (canonical_xyz)
      if (allocated(canonical_eread)) deallocate (canonical_eread)
      if (allocated(canonical_comments)) deallocate (canonical_comments)
      if (allocated(canonical_frame_status)) deallocate (canonical_frame_status)
      return
    end if
    write (stdout,'(1x,a,i0,a)') 'Optimizer affinity completed: ',T, &
    & ' workers bound and restored.'
  end if

!>--- finalize progress printout
  call crest_oloop_pr_progress(env,nall,-1)

!>--- stop timer
  call profiler%stop(1)

!>--- prepare some summary printout
  percent = float(c)/float(nall)*100.0_wp
  write (atmp,'(f5.1,a)') percent,'% success)'
  write (stdout,'(">",1x,i0,a,i0,a,a)') c,' of ',nall,' structures successfully optimized (', &
  &     trim(adjustl(atmp))
  write (atmp,'(">",1x,a,i0,a)') 'Total runtime for ',nall,' optimizations:'
  call profiler%write_timing(stdout,1,trim(atmp),.true.)
  runtime = profiler%get(1)
  write (atmp,'(f16.3,a)') runtime/real(nall,wp),' sec'
  write (stdout,'(a,a,a)') '> Corresponding to approximately ',trim(adjustl(atmp)), &
  &                       ' per processed structure'

  if (canonical_output.and.ndump > 0) then
    canonical_message = ''
    canonical_artifact_bytes = 0_int64
    canonical_write_seconds = 0.0_wp
    canonical_workers = max(1,min(16,ndump,omp_get_max_threads()))
    canonical_team = 1
    canonical_fallbacks = 0_int64
    poststage_clock = omp_get_wtime()
!$omp parallel default(none) num_threads(canonical_workers) &
!$omp shared(ndump,nat,at,canonical_xyz,canonical_eread,canonical_comments) &
!$omp shared(canonical_frame_status,canonical_team) private(i,canonical_frame_fallbacks) &
!$omp reduction(+:canonical_fallbacks)
!$omp single
    canonical_team = omp_get_num_threads()
!$omp end single
!$omp do schedule(static)
    do i = 1,ndump
      call canonicalize_empty_refine(nat,canonical_xyz(:,:,i), &
      & canonical_comments(i),canonical_eread(i),canonical_frame_status(i), &
      & canonical_frame_fallbacks)
      canonical_fallbacks = canonical_fallbacks+canonical_frame_fallbacks
    end do
!$omp end do
!$omp end parallel
    canonical_iostat = 0
    do i = 1,ndump
      if (canonical_frame_status(i) /= 0) then
        canonical_iostat = canonical_frame_status(i)
        exit
      end if
    end do
    write (stdout,'(1x,a,f12.3,a)') &
    & 'Poststage empty-refinement parallel canonicalization wall time: ', &
    & omp_get_wtime()-poststage_clock,' sec'
    write (stdout,'(1x,a,i0,a,i0)') &
    & 'Poststage empty-refinement canonicalization workers requested/actual: ', &
    & canonical_workers,' / ',canonical_team
    write (stdout,'(1x,a,i0)') &
    & 'Poststage fixed-decimal near-boundary fallbacks: ',canonical_fallbacks

    if (canonical_iostat == 0) then
      poststage_clock = omp_get_wtime()
      call write_canonical_ensemble_fast(ensemblefile,at,canonical_xyz(:,:,1:ndump), &
      & canonical_comments(1:ndump),canonical_artifact_bytes,canonical_write_seconds, &
      & canonical_status,canonical_message)
      canonical_iostat = canonical_status
      write (stdout,'(1x,a,f12.3,a)') &
      & 'Poststage buffered C ordered artifact write wall time: ', &
      & omp_get_wtime()-poststage_clock,' sec'
      write (stdout,'(1x,a,i0,a,f12.3,a)') &
      & 'Poststage buffered C ordered artifact bytes: ',canonical_artifact_bytes, &
      & '; writer-internal time: ',canonical_write_seconds,' sec'
    end if

    if (canonical_iostat /= 0) then
      write (stdout,'(1x,a,i0)') &
      & 'Canonical empty-refinement serialization failed; iostat=',canonical_iostat
      if (len_trim(canonical_message) > 0) write(stdout,'(1x,a)') trim(canonical_message)
      env%iostatus_meta = status_failed
      close (ich2,iostat=canonical_status)
      call remove(ensemblefile)
      call remove(ensembleelog)
      call profiler%clear()
      deallocate (grads,calculations)
      if (allocated(mols)) deallocate (mols)
      if (allocated(molsnew)) deallocate (molsnew)
      if (allocated(canonical_xyz)) deallocate (canonical_xyz)
      if (allocated(canonical_eread)) deallocate (canonical_eread)
      if (allocated(canonical_comments)) deallocate (canonical_comments)
      if (allocated(canonical_frame_status)) deallocate (canonical_frame_status)
      return
    end if
  end if

!>--- close files (if they are open)
  canonical_status = 0
  if (dump) then
    if (canonical_output) then
      close (ich2,iostat=io)
      if (io /= 0 .and. canonical_status == 0) canonical_status = io
    else
      close (ich)
      close (ich2)
    end if
  end if

  if (canonical_output.and.canonical_status /= 0) then
    write (stdout,'(1x,a,i0)') &
    & 'Canonical empty-refinement artifact close failed; iostat=',canonical_status
    env%iostatus_meta = status_failed
    call poststage_buffer%clear()
    call remove(ensemblefile)
    call remove(ensembleelog)
    if (allocated(canonical_xyz)) deallocate (canonical_xyz)
    if (allocated(canonical_eread)) deallocate (canonical_eread)
    if (allocated(canonical_comments)) deallocate (canonical_comments)
    if (allocated(canonical_frame_status)) deallocate (canonical_frame_status)
    deallocate (grads)
    call profiler%clear()
    deallocate (calculations)
    if (allocated(mols)) deallocate (mols)
    if (allocated(molsnew)) deallocate (molsnew)
    return
  end if

  if (canonical_output.and.ndump > 0) then
    poststage_buffer%nat = nat
    poststage_buffer%nall = ndump
    allocate (poststage_buffer%at(nat),stat=io)
    if (io /= 0) then
      call poststage_buffer%clear()
      env%iostatus_meta = status_failed
      call remove(ensemblefile)
      call remove(ensembleelog)
      if (allocated(canonical_xyz)) deallocate (canonical_xyz)
      if (allocated(canonical_eread)) deallocate (canonical_eread)
      if (allocated(canonical_comments)) deallocate (canonical_comments)
      if (allocated(canonical_frame_status)) deallocate (canonical_frame_status)
      deallocate (grads)
      call profiler%clear()
      deallocate (calculations)
      if (allocated(mols)) deallocate (mols)
      if (allocated(molsnew)) deallocate (molsnew)
      return
    end if
    poststage_buffer%at = at
    if (ndump == nall) then
      call move_alloc(canonical_xyz,poststage_buffer%xyz)
      call move_alloc(canonical_eread,poststage_buffer%eread)
      call move_alloc(canonical_comments,poststage_buffer%comments)
    else
      allocate (poststage_buffer%xyz(3,nat,ndump), &
      & poststage_buffer%eread(ndump),poststage_buffer%comments(ndump),stat=io)
      if (io /= 0) then
        call poststage_buffer%clear()
        env%iostatus_meta = status_failed
        call remove(ensemblefile)
        call remove(ensembleelog)
        if (allocated(canonical_xyz)) deallocate (canonical_xyz)
        if (allocated(canonical_eread)) deallocate (canonical_eread)
        if (allocated(canonical_comments)) deallocate (canonical_comments)
        if (allocated(canonical_frame_status)) deallocate (canonical_frame_status)
        deallocate (grads)
        call profiler%clear()
        deallocate (calculations)
        if (allocated(mols)) deallocate (mols)
        if (allocated(molsnew)) deallocate (molsnew)
        return
      end if
      poststage_buffer%xyz = canonical_xyz(:,:,1:ndump)
      poststage_buffer%eread = canonical_eread(1:ndump)
      poststage_buffer%comments = canonical_comments(1:ndump)
    end if
    if (present(poststage_written)) poststage_written = .true.
  end if

  if (allocated(canonical_xyz)) deallocate (canonical_xyz)
  if (allocated(canonical_eread)) deallocate (canonical_eread)
  if (allocated(canonical_comments)) deallocate (canonical_comments)
  if (allocated(canonical_frame_status)) deallocate (canonical_frame_status)

  deallocate (grads)
  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  if (allocated(molsnew)) deallocate (molsnew)
  return
contains

  subroutine canonicalize_empty_refine(nat,xyz_canonical, &
  & comment_canonical,energy_canonical,status,fallbacks)
!********************************************************************************
!* Reproduce the fixed-decimal values produced by coord%append followed by an
!* empty crest_refine pass. Ordinary values use a numerical quantization kernel;
!* only near-halfway values enter the formatted-I/O compatibility fallback.
!********************************************************************************
    integer,intent(in) :: nat
    real(wp),intent(inout) :: xyz_canonical(3,nat)
    real(wp),intent(out) :: energy_canonical
    character(len=poststage_comment_length),intent(inout) :: comment_canonical
    integer,intent(out) :: status
    integer(int64),intent(out) :: fallbacks
    integer :: j,k,io
    real(wp) :: first_rounded,second_input,rounded
    real(wp) :: energy_first
    character(len=poststage_comment_length) :: optimizer_comment
    logical :: used_fallback,first_fallback

    status = 0
    fallbacks = 0_int64
    energy_canonical = 0.0_wp

    optimizer_comment = comment_canonical
    energy_first = grepenergy(optimizer_comment)
    write (comment_canonical,'(2x,f18.8)',iostat=io) energy_first
    if (io /= 0) then
      status = io
      return
    end if
    energy_canonical = grepenergy(comment_canonical)

    do j = 1,nat
      do k = 1,3
        call round_fixed_decimal_10(xyz_canonical(k,j)*bohr,first_rounded,io, &
        & first_fallback)
        if (io /= 0) then
          status = io
          return
        end if
        ! Preserve the two distinct legacy assignment/unit-conversion stages.
        ! Algebraically this is close to first_rounded, but the separate
        ! operations can matter at an F20.10 half-way boundary.
        second_input = (first_rounded/bohr)*bohr
        call round_fixed_decimal_10(second_input,rounded,io,used_fallback)
        if (io /= 0) then
          status = io
          return
        end if
        xyz_canonical(k,j) = rounded
        if (first_fallback) fallbacks = fallbacks+1_int64
        if (used_fallback) fallbacks = fallbacks+1_int64
      end do
    end do
  end subroutine canonicalize_empty_refine
end subroutine crest_oloop

!========================================================================================!
subroutine crest_oloop_pr_progress(env,total,current)
!*********************************************
!* subroutine crest_oloop_pr_progress
!* A subroutine to print and track progress of
!* concurrent geometry optimizations
!*********************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use iomod,only:to_str
  implicit none
  type(systemdata),intent(inout) :: env
  integer,intent(in) :: total,current
  real(wp) :: percent
  character(len=5) :: atmp
  real(wp),save :: increment
  real(wp),save :: progressbarrier

  percent = float(current)/float(total)*100.0_wp
  if (current == 0) then !> as a wrapper to start the printout
    progressbarrier = 0.0_wp
    if (env%niceprint) then
      percent = 0.0_wp
      call printprogbar(percent)
    end if
    increment = 10.0_wp
    if (total > 1000) increment = 7.5_wp
    if (total > 5000) increment = 5.0_wp
    if (total > 10000) increment = 2.5_wp
    if (total > 20000) increment = 1.0_wp

  else if (current <= total.and.current > 0) then !> the regular printout case
    if (env%niceprint) then
      call printprogbar(percent)

    else if (.not.env%legacy) then
      if (percent >= progressbarrier) then
        write (atmp,'(f5.1)') percent
        write (stdout,'(1x,a)',advance='no') '|>'//trim(adjustl(atmp))//'%'
        progressbarrier = progressbarrier+increment
        progressbarrier = min(progressbarrier,100.0_wp)
        flush (stdout)
      end if
    else
      write (stdout,'(1x,i0)',advance='no') current
      flush (stdout)
    end if

  else !> as a wrapper to finalize the printout
    if (.not.env%niceprint) then
      write (stdout,'(/,1x,a)') 'done.'
    else
      write (stdout,*)
    end if
  end if

end subroutine crest_oloop_pr_progress

!========================================================================================!
!========================================================================================!
!> Routines for parallel MDs
!========================================================================================!
!========================================================================================!

subroutine crest_search_multimd(env,mol,mddats,nsim)
!*****************************************************
!* subroutine crest_search_multimd
!* this runs #nsim MDs on the same structure (mol)
!*****************************************************
  use crest_parameters,only:wp,stdout,sep
  use iso_fortran_env,only:int64
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove,collect_stream_files_exact
  use omp_lib
  use crest_restartlog,only:trackrestart,restart_write_dummy
  use parallel_interface,only:crest_search_multimd2
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mol
  type(coord),allocatable :: moltmps(:),process_mols(:)
  integer :: i,j,io,ich,collector_status,path_length
  logical :: pr,ex,nested,use_tmd_threads
  logical :: all_gfnff,has_active_calc
  integer :: T,Tn,Trestore,Tnrestore,thread_save
  real(wp) :: percent
  real(wp) :: collector_seconds
  integer(int64) :: collector_bytes
  character(len=80) :: atmp
  character(len=1024) :: collector_message
  character(len=:),allocatable :: trajectory_paths(:)
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  real(wp) :: etmp
  real(wp),allocatable :: grdtmp(:,:,:)
  type(timer) :: profiler
!===========================================================!
!>--- decide wether to skip this call
  if (trackrestart(env)) then
    call restart_write_dummy('crest_dynamics.trj')
    return
  end if

!>--- check if we have any MD & calculation settings allocated
  if (.not.env%mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  ! All multi-trajectory same-input MTD work is delegated to the modern
  ! process-isolated scheduler used by crest_search_multimd2.  Keep only the
  ! single-trajectory direct path below.
  if (nsim > 1) then
    allocate(process_mols(nsim),source=mol)
    call crest_search_multimd2(env,process_mols,mddats,nsim)
    return
  end if

!>--- prepare calculation containers for parallelization (one per thread)
  use_tmd_threads = env%threadsmdsetmanual.and.env%ThreadsMD > 0
  thread_save = env%Threads
  if (use_tmd_threads) env%Threads = env%ThreadsMD

  all_gfnff = .true.
  has_active_calc = .false.
  do j = 1,env%calc%ncalculations
    if (env%calc%calcs(j)%active) then
      has_active_calc = .true.
      if (env%calc%calcs(j)%id /= jobtype%gfnff) all_gfnff = .false.
    end if
  end do
  all_gfnff = all_gfnff.and.has_active_calc
  call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  nested = env%omp_allow_nested
  if (all_gfnff) then
    if (.not.nested) Tn = 1
    write (stdout,'(1x,a,i0,a,i0)') &
    & 'GFN-FF MTD scheduling: parallel trajectories=',T,', cores per trajectory=',Tn
  end if

  allocate (grdtmp(3,mol%nat,T),source=0.0_wp)
  if (all_gfnff) then
    !> Create only the outer containers here.  Each bound outer worker performs
    !> the deep copies below, so its GFN-FF topology and workspace pages are
    !> first-touched in that worker's NUMA locality.
    allocate (calculations(T))
    allocate (moltmps(T))
    do j = 1,env%calc%ncalculations
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) io = makedir(trim(env%calc%calcs(j)%calcspace))
    end do
  else
    !> Preserve the existing setup path for mixed/external calculators.
    allocate (calculations(T),source=env%calc)
    allocate (moltmps(T),source=mol)
    do i = 1,T
      moltmps(i)%nat = mol%nat
      moltmps(i)%at = mol%at
      moltmps(i)%xyz = mol%xyz
      do j = 1,env%calc%ncalculations
        calculations(i)%calcs(j) = env%calc%calcs(j)
        ex = directory_exist(env%calc%calcs(j)%calcspace)
        if (.not.ex) io = makedir(trim(env%calc%calcs(j)%calcspace))
        write (atmp,'(a,"_",i0)') sep,i
        calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
        if (allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
        if (allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
        call calculations(i)%calcs(j)%printid(i,j)
      end do
      calculations(i)%pr_energies = .false.
      call engrad(moltmps(i),calculations(i),etmp,grdtmp(:,:,i),io)
    end do
  end if

  !>--- other settings
  pr = .false.
  call profiler%init(nsim)

  !>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mol,pr,percent,ich,nsim,moltmps,nested,Tn,grdtmp,all_gfnff) &
  !!$omp single
  !$omp private(vz,i,j,job,thread_id,io,ex,atmp,etmp)
  thread_id = omp_get_thread_num()
  job = thread_id+1
  if (all_gfnff) then
    calculations(job) = env%calc
    moltmps(job) = mol
    do j = 1,env%calc%ncalculations
      write (atmp,'(a,"_",i0)') sep,job
      calculations(job)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if (allocated(calculations(job)%calcs(j)%calcfile)) deallocate(calculations(job)%calcs(j)%calcfile)
      if (allocated(calculations(job)%calcs(j)%systemcall)) deallocate(calculations(job)%calcs(j)%systemcall)
      call calculations(job)%calcs(j)%printid(job,j)
    end do
    calculations(job)%pr_energies = .false.

    !> Preserve the existing one-initialization-evaluation-per-worker behavior.
    !> GFN-FF setup remains single-threaded inside each independent outer worker.
    call ompmklset(1)
    call engrad(moltmps(job),calculations(job),etmp,grdtmp(:,:,job),io)
  end if
  if (nested) call ompmklset(Tn)

  !$omp do
  do i = 1,nsim

    call initsignal()
    vz = i

    !!$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !$omp critical
    moltmps(job)%nat = mol%nat
    moltmps(job)%at = mol%at
    moltmps(job)%xyz = mol%xyz
    !$omp end critical
    !>--- startup printout (thread safe)
    call parallel_md_block_printout(mddats(vz),vz)

    !>--- the acutal MD call with timing
    call profiler%start(vz)
    !> Each worker calculator persists across independent MTD trajectories.
    !> Rebuild only its geometry-dependent GFN-FF HB/XB list at this boundary;
    !> retain topology, exact EEQ/frozen-host caches, D3, and all workspaces.
    call request_gfnff_hbond_update(calculations(job))
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)
    mddats(vz)%termination_status = io
    call profiler%stop(vz)

    !>--- finish printout (thread safe)
    call parallel_md_finish_printout(mddats(vz),vz,io,profiler)
    !!$omp end task
  end do
  !!$omp taskwait
  !$omp end parallel

  !>--- collect trajectories into one exact byte stream.  Construct the source
  !> list explicitly in numeric mddats order; the shared collector uses one
  !> writer and does not parallelize or reformat scientific trajectory data.
  collector_status = 0
  collector_message = ''
  collector_seconds = 0.0_wp
  collector_bytes = 0_int64
  path_length = 1
  do i = 1,nsim
    if (.not.allocated(mddats(i)%trajectoryfile)) then
      collector_status = 1
      write (collector_message,'(a,i0)') &
      & 'trajectory path is not allocated at numeric index ',i
      exit
    end if
    path_length = max(path_length,len(mddats(i)%trajectoryfile))
  end do
  if (collector_status == 0) then
    allocate (character(len=path_length) :: trajectory_paths(nsim), &
    &         stat=collector_status,errmsg=collector_message)
  end if
  if (collector_status == 0) then
    do i = 1,nsim
      trajectory_paths(i) = mddats(i)%trajectoryfile
    end do
    call collect_stream_files_exact(trajectory_paths,'crest_dynamics.trj', &
    & collector_bytes,collector_seconds,collector_status,collector_message)
  end if
  write (stdout,'(1x,a,f12.3,a)') &
  & 'Trajectory collector wall time: ',collector_seconds,' sec'
  write (stdout,'(1x,a,i0)') 'Trajectory collector bytes: ',collector_bytes
  if (collector_status /= 0) then
    write (stdout,'(1x,a,i0,2a)') 'Trajectory collector failed; iostat=', &
    & collector_status,': ',trim(collector_message)
    env%iostatus_meta = status_failed
  end if
  if (allocated(trajectory_paths)) deallocate (trajectory_paths)

  call profiler%clear()
  deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  if (use_tmd_threads.or.all_gfnff) then
    env%Threads = thread_save
    call new_ompautoset(env,'max',0,Trestore,Tnrestore)
  end if
  return
end subroutine crest_search_multimd

!========================================================================================!
subroutine crest_search_multimd_init(env,mol,mddat,nsim)
!*******************************************************
!* subroutine crest_search_multimd_init
!* This routine will initialize a copy of env%mddat
!* and save it to the local mddat. If we are about to
!* run RMSD metadynamics, the required number of
!* simulations (#nsim) is returned
!*******************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddat
  type(coord) :: mol
  integer,intent(inout) :: nsim
  integer :: i,io
  logical :: pr
!=======================================================!
  type(calcdata),target :: calc
  type(shakedata) :: shk

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)
  character(len=*),parameter :: mdir = 'MDFILES'
!======================================================!

  !>--- check if we have any MD & calculation settings allocated
  mddat = env%mddat
  if (.not.mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  !>--- init SHAKE?
  if (mddat%shake) then
    if (allocated(env%ref%wbo)) then
      shk%wbo = env%ref%wbo
    else
      calc = env%calc
      calc%calcs(1)%rdwbo = .true.
      allocate (grad(3,mol%nat),source=0.0_wp)
      call engrad(mol,calc,energy,grad,io)
      deallocate (grad)
      calc%calcs(1)%rdwbo = .false.

      shk%shake_mode = env%mddat%shk%shake_mode
      call move_alloc(calc%calcs(1)%wbo,shk%wbo)
    end if

    if (calc%nfreeze > 0) then
      shk%freezeptr => calc%freezelist
    else
      nullify (shk%freezeptr)
    end if

    shk%shake_mode = env%shake
    mddat%shk = shk
    call init_shake(mol%nat,mol%at,mol%xyz,mddat%shk,pr)
    mddat%nshake = mddat%shk%ncons
  end if
  !>--- complete real-time settings to steps
  call mdautoset(mddat,io)

  !>--- (optional)  MTD initialization
  if (nsim < 0) then
    mddat%simtype = type_mtd  !>-- set runtype to MTD

    call defaultGF(env)
    write (stdout,*) 'list of applied metadynamics Vbias parameters:'
    do i = 1,env%nmetadyn
      write (stdout,'(''$metadyn '',f10.5,f8.3,i5)') env%metadfac(i),env%metadexp(i)
    end do
    write (stdout,*)

    !>--- how many simulations
    nsim = env%nmetadyn
  end if

  return
end subroutine crest_search_multimd_init

!========================================================================================!
subroutine crest_search_multimd_init2(env,mddats,nsim)
  use crest_data
  use dynamics_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  integer,allocatable :: bias_indices(:),input_indices(:)
  integer :: i
  allocate(bias_indices(nsim),input_indices(nsim))
  do i=1,nsim
    bias_indices(i)=i
  end do
  input_indices=1
  call crest_search_multimd_init2_mapped(env,mddats,nsim,bias_indices,input_indices)
  deallocate(bias_indices,input_indices)
end subroutine crest_search_multimd_init2

!========================================================================================!
subroutine crest_search_multimd_init2_mapped(env,mddats,nsim,bias_indices,input_indices)
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  integer,intent(in) :: bias_indices(nsim),input_indices(nsim)
  integer :: i,io,j,bias_id
  logical :: ex
  type(mtdpot),allocatable :: mtds(:)
  character(len=80) :: atmp
  character(len=*),parameter :: mdir='MDFILES'

  ex=directory_exist(mdir)
  if (ex) call rmrf(mdir)
  io=makedir(mdir)
  do i=1,nsim
    if (input_indices(i)<1) error stop '**ERROR** invalid input structure index for MD job'
    mddats(i)%md_index=i
    mddats(i)%input_structure_id=input_indices(i)
    mddats(i)%bias_configuration_id=bias_indices(i)
    mddats(i)%termination_status=-1
    write(atmp,'(a,i0,a)') 'crest_',i,'.trj'
    mddats(i)%trajectoryfile=mdir//sep//trim(atmp)
    write(atmp,'(a,i0,a)') 'crest_',i,'.mdrestart'
    mddats(i)%restartfile=mdir//sep//trim(atmp)
  end do

  allocate(mtds(nsim))
  do i=1,nsim
    if (mddats(i)%simtype == type_mtd) then
      bias_id=mddats(i)%bias_configuration_id
      if (bias_id<1 .or. bias_id>env%nmetadyn) then
        error stop '**ERROR** invalid metadynamics bias index for MD job'
      end if
      mtds(i)%kpush=env%metadfac(bias_id)
      mtds(i)%alpha=env%metadexp(bias_id)
      mtds(i)%cvdump_fs=float(env%mddump)
      mtds(i)%mtdtype=cv_rmsd
      mtds(i)%com_bias=env%mtd_com_bias
      mtds(i)%com_factor=env%mtd_com_factor
      mtds(i)%com_width=env%mtd_com_width
      mtds(i)%com_mass_weighted=env%mtd_com_mass_weighted
      mddats(i)%npot=1
      allocate(mddats(i)%mtd(1),source=mtds(i))
      allocate(mddats(i)%cvtype(1),source=cv_rmsd)
      if (sum(env%includeRMSD) /= env%ref%nat) then
        if (.not.allocated(mddats(i)%mtd(1)%atinclude)) &
          allocate(mddats(i)%mtd(1)%atinclude(env%ref%nat),source=.true.)
        do j=1,env%ref%nat
          if (env%includeRMSD(j)/=1) mddats(i)%mtd(1)%atinclude(j)=.false.
        end do
      end if
    end if
  end do
  deallocate(mtds)
end subroutine crest_search_multimd_init2_mapped

!========================================================================================!
subroutine crest_search_multimd2(env,mols,mddats,nsim,poststage_out)
!*******************************************************************
!* subroutine crest_search_multimd2
!* this runs #nsim MDs on #nsim selected different structures (mols)
!*******************************************************************
  use crest_parameters,only:wp,stdout,sep
  use iso_fortran_env,only:int64
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove,collect_stream_files_exact
  use omp_lib
  use crest_restartlog,only:trackrestart,restart_write_dummy
  use mtd_process_scheduler,only:resolve_mtd_process_isolation,run_mtd_process_batch
  use crest_poststage_ensemble,only:poststage_ensemble
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(mddata),intent(inout) :: mddats(nsim)
  integer,intent(in) :: nsim
  type(coord),intent(in) :: mols(nsim)
  type(poststage_ensemble),intent(inout),optional :: poststage_out
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich,collector_status,path_length
  logical :: pr,ex,nested,use_tmd_threads,process_isolated
  logical :: all_gfnff,has_active_calc
  integer :: T,Tn,Trestore,Tnrestore,thread_save
  real(wp) :: percent
  real(wp) :: collector_seconds
  integer(int64) :: collector_bytes
  character(len=80) :: atmp
  character(len=1024) :: collector_message
  character(len=1024) :: process_message
  character(len=:),allocatable :: trajectory_paths(:)
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  integer :: process_status
  type(timer) :: profiler
!===========================================================!
  if (present(poststage_out)) call poststage_out%clear()
!>--- decide wether to skip this call
  if (trackrestart(env)) then
    call restart_write_dummy('crest_dynamics.trj')
    return
  end if

!>--- check if we have any MD & calculation settings allocated
  if (.not.env%mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
  use_tmd_threads = env%threadsmdsetmanual.and.env%ThreadsMD > 0
  thread_save = env%Threads
  if (use_tmd_threads) env%Threads = env%ThreadsMD

  all_gfnff = .true.
  has_active_calc = .false.
  do j = 1,env%calc%ncalculations
    if (env%calc%calcs(j)%active) then
      has_active_calc = .true.
      if (env%calc%calcs(j)%id /= jobtype%gfnff) all_gfnff = .false.
    end if
  end do
  all_gfnff = all_gfnff.and.has_active_calc
  call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  nested = env%omp_allow_nested
  if (all_gfnff) then
    if (.not.nested) Tn = 1
    write (stdout,'(1x,a,i0,a,i0)') &
    & 'GFN-FF MTD scheduling: parallel trajectories=',T,', cores per trajectory=',Tn
  end if

  !> The process scheduler is the fail-closed default for every multi-MTD
  !> batch in this routine.  All molecule/MTD preparation
  !> has already happened in the parent, and collection below remains unchanged.
  call resolve_mtd_process_isolation(process_isolated,process_status,process_message)
  if (process_status /= status_normal) then
    write(stdout,'(1x,a)') 'Process-isolated MTD option failed closed: '//trim(process_message)
    env%iostatus_meta = process_status
    if (use_tmd_threads.or.all_gfnff) then
      env%Threads = thread_save
      call new_ompautoset(env,'max',0,Trestore,Tnrestore)
    end if
    return
  end if
  if (process_isolated .and. nsim > 1) then
    if (present(poststage_out)) then
      call run_mtd_process_batch(env,mols,mddats,nsim,T,Tn,process_status, &
      & process_message,poststage_out)
    else
      call run_mtd_process_batch(env,mols,mddats,nsim,T,Tn,process_status,process_message)
    end if
    if (process_status /= status_normal) then
      write(stdout,'(1x,a)') 'Process-isolated MTD failed closed: '//trim(process_message)
      env%iostatus_meta = process_status
      if (use_tmd_threads.or.all_gfnff) then
        env%Threads = thread_save
        call new_ompautoset(env,'max',0,Trestore,Tnrestore)
      end if
      return
    end if
    goto 800
  end if
  !> Native/threaded MTD retains the historical file handoff.  The optional
  !> carrier was cleared on entry, so its caller will deterministically fall
  !> back to crest_dynamics.trj after this branch completes.

  if (all_gfnff) then
    !> Leave each derived-type element unallocated until its bound outer worker
    !> deep-copies it below.  This first-touches private topology/workspace pages
    !> in the same NUMA locality that will run the trajectory.
    allocate (calculations(T))
    allocate (moltmps(T))
    do j = 1,env%calc%ncalculations
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) io = makedir(trim(env%calc%calcs(j)%calcspace))
    end do
  else
    !> Preserve the existing setup path for mixed/external calculators.
    allocate (calculations(T),source=env%calc)
    allocate (moltmps(T),source=mols(1))
    do i = 1,T
      do j = 1,env%calc%ncalculations
        calculations(i)%calcs(j) = env%calc%calcs(j)
        ex = directory_exist(env%calc%calcs(j)%calcspace)
        if (.not.ex) io = makedir(trim(env%calc%calcs(j)%calcspace))
        write (atmp,'(a,"_",i0)') sep,i
        calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
        if (allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
        if (allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
        call calculations(i)%calcs(j)%printid(i,j)
      end do
      calculations(i)%pr_energies = .false.
    end do
  end if

!>--- other settings
  pr = .false.
  call profiler%init(nsim)

!>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mols,pr,percent,ich,moltmps,profiler,nested,Tn,all_gfnff) &
  !$omp private(i,j,vz,job,thread_id,io,ex,atmp)
  thread_id = omp_get_thread_num()
  job = thread_id+1
  if (all_gfnff) then
    calculations(job) = env%calc
    moltmps(job) = mols(1)
    do j = 1,env%calc%ncalculations
      write (atmp,'(a,"_",i0)') sep,job
      calculations(job)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if (allocated(calculations(job)%calcs(j)%calcfile)) deallocate(calculations(job)%calcs(j)%calcfile)
      if (allocated(calculations(job)%calcs(j)%systemcall)) deallocate(calculations(job)%calcs(j)%systemcall)
      call calculations(job)%calcs(j)%printid(job,j)
    end do
    calculations(job)%pr_energies = .false.
  end if
  if (nested) call ompmklset(Tn)
  !$omp barrier
  !$omp single
  do i = 1,nsim

    call initsignal()
    vz = i

    !$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !$omp critical
    moltmps(job)%nat = mols(vz)%nat
    moltmps(job)%at = mols(vz)%at
    moltmps(job)%xyz = mols(vz)%xyz
    !$omp end critical
    !>--- startup printout (thread safe)
    call parallel_md_block_printout(mddats(vz),vz)

    !>--- the acutal MD call with timing
    call profiler%start(vz)
    !> Each worker calculator persists across independent multi-input MTD
    !> trajectories. Rebuild only its geometry-dependent GFN-FF HB/XB list at
    !> this boundary; retain topology, exact EEQ/frozen-host caches, D3, and all
    !> workspaces.
    call request_gfnff_hbond_update(calculations(job))
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)
    mddats(vz)%termination_status = io
    call profiler%stop(vz)

    !>--- finish printout (thread safe)
    call parallel_md_finish_printout(mddats(vz),vz,io,profiler)
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel

!>--- collect trajectories into one exact byte stream.  Construct the source
!> list explicitly in numeric mddats order; the shared collector uses one
!> writer and does not parallelize or reformat scientific trajectory data.
800 continue
  collector_status = 0
  collector_message = ''
  collector_seconds = 0.0_wp
  collector_bytes = 0_int64
  path_length = 1
  do i = 1,nsim
    if (.not.allocated(mddats(i)%trajectoryfile)) then
      collector_status = 1
      write (collector_message,'(a,i0)') &
      & 'trajectory path is not allocated at numeric index ',i
      exit
    end if
    path_length = max(path_length,len(mddats(i)%trajectoryfile))
  end do
  if (collector_status == 0) then
    allocate (character(len=path_length) :: trajectory_paths(nsim), &
    &         stat=collector_status,errmsg=collector_message)
  end if
  if (collector_status == 0) then
    do i = 1,nsim
      trajectory_paths(i) = mddats(i)%trajectoryfile
    end do
    call collect_stream_files_exact(trajectory_paths,'crest_dynamics.trj', &
    & collector_bytes,collector_seconds,collector_status,collector_message)
  end if
  write (stdout,'(1x,a,f12.3,a)') &
  & 'Trajectory collector wall time: ',collector_seconds,' sec'
  write (stdout,'(1x,a,i0)') 'Trajectory collector bytes: ',collector_bytes
  if (collector_status /= 0) then
    write (stdout,'(1x,a,i0,2a)') 'Trajectory collector failed; iostat=', &
    & collector_status,': ',trim(collector_message)
    env%iostatus_meta = status_failed
    if (present(poststage_out)) call poststage_out%clear()
  end if
  if (allocated(trajectory_paths)) deallocate (trajectory_paths)

  call profiler%clear()
  if (allocated(calculations)) deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  if (use_tmd_threads.or.all_gfnff) then
    env%Threads = thread_save
    call new_ompautoset(env,'max',0,Trestore,Tnrestore)
  end if
  return
end subroutine crest_search_multimd2

!========================================================================================!
subroutine parallel_md_block_printout(MD,vz)
!***********************************************
!* subroutine parallel_md_block_printout
!* This will print information about the MD/MTD
!* simulation. The execution is omp threadsave
!***********************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod,only:to_str
  implicit none
  type(mddata),intent(in) :: MD
  integer,intent(in) :: vz
  character(len=40) :: atmp
  integer :: il
  !$omp critical

  if (MD%simtype == type_md) then
    write (atmp,'(a,1x,i3)') 'starting MD',vz
  else if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd_static) then
      write (atmp,'(a,1x,i3)') 'starting static MTD',vz
    else
      write (atmp,'(a,1x,i4)') 'starting MTD',vz
    end if
  end if
  il = (44-len_trim(atmp))/2
  write (stdout,'(2x,a,1x,a,1x,a)') repeat(':',il),trim(atmp),repeat(':',il)

  write (stdout,'(2x,"|   MD simulation time   :",f8.1," ps       |")') MD%length_ps
  write (stdout,'(2x,"|   target T             :",f8.1," K        |")') MD%tsoll
  write (stdout,'(2x,"|   timestep dt          :",f8.1," fs       |")') MD%tstep
  write (stdout,'(2x,"|   dump interval(trj)   :",f8.1," fs       |")') MD%dumpstep
  if (MD%shake.and.MD%shk%shake_mode > 0) then
    if (MD%shk%shake_mode == 2) then
      write (stdout,'(2x,"|   SHAKE algorithm      :",a5," (all bonds) |")') to_str(MD%shake)
    else
      write (stdout,'(2x,"|   SHAKE algorithm      :",a5," (H only) |")') to_str(MD%shake)
    end if
  end if
  if (allocated(MD%active_potentials)) then
    write (stdout,'(2x,"|   active potentials    :",i4," potential    |")') size(MD%active_potentials,1)
  end if
  if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd) then
      write (stdout,'(2x,"|   dump interval(Vbias) :",f8.2," ps       |")') &
          & MD%mtd(1)%cvdump_fs/1000.0_wp
    end if
    write (stdout,'(2x,"|   Vbias prefactor (k)  :",f8.4," Eh       |")') &
      &  MD%mtd(1)%kpush
    if (MD%cvtype(1) == cv_rmsd.or.MD%cvtype(1) == cv_rmsd_static) then
      write (stdout,'(2x,"|   Vbias exponent (α)   :",f8.4," bohr⁻²   |")') MD%mtd(1)%alpha
    else
      write (stdout,'(2x,"|   Vbias exponent (α)   :",f8.4,"          |")') MD%mtd(1)%alpha
    end if
  end if

  !$omp end critical

end subroutine parallel_md_block_printout

subroutine parallel_md_finish_printout(MD,vz,io,profiler)
!*******************************************
!* subroutine parallel_md_finish_printout
!* This will print information termination
!* info about the MD/MTD simulation
!*******************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  implicit none
  type(mddata),intent(in) :: MD
  integer,intent(in) :: vz,io
  type(timer),intent(inout) :: profiler
  character(len=40) :: atmp
  character(len=80) :: btmp

  !$omp critical

  if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd_static) then
      write (atmp,'(a)') '*sMTD'
    else
      write (atmp,'(a)') '*MTD'
    end if
  else
    write (atmp,'(a)') '*MD'
  end if
  if (io == 0) then
    write (btmp,'(a,1x,i3,a)') trim(atmp),vz,' completed successfully'
  else
    write (btmp,'(a,1x,i3,a)') trim(atmp),vz,' terminated EARLY'
  end if
  call profiler%write_timing(stdout,vz,trim(btmp))

  !$omp end critical

end subroutine parallel_md_finish_printout
!========================================================================================!
!========================================================================================!
