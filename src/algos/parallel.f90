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
    subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump,customcalc)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
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
    end subroutine crest_oloop
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
subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump,customcalc)
!*******************************************************************************
!* subroutine crest_oloop
!* This subroutine performs concurrent geometry optimizations
!* for the given ensemble. Inputs xyz and eread are overwritten
!* env        - contains parallelization and other program settings
!* dump       - decides on whether to dump an ensemble file
!*              WARNING: the ensemble file will NOT be in the same order
!*              as the input xyz array. However, the overwritten xyz will be! 
!* customcalc - customized (optional) calculation level data
!*
!* IMPORTANT: xyz should be in Bohr(!) for this routine
!******************************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
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

  type(coord),allocatable :: mols(:)
  type(coord),allocatable :: molsnew(:)
  type(coord) :: canonical_mol
  integer :: i,j,k,l,io,ich,ich2,c,z,job_id,zcopy,ninitialized,ngfnff
  logical :: pr,wr,ex
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grads(:,:,:)
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent,runtime
  type(calcdata),pointer :: mycalc
  type(calcdata),target :: local_template
  type(timer) :: profiler
  integer :: T,Tn  !> threads and threads per core
  logical :: nested

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
    open (newunit=ich,file=ensemblefile)
    open (newunit=ich2,file=ensembleelog)
  end if
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
  !$omp shared(env,calculations,nat,nall,at,xyz,eread,grads,c,k,z,pr,wr,dump) &
  !$omp shared(ich,ich2,mols,molsnew, nested,Tn)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(j,job,energy,io,atmp,gnorm,thread_id,zcopy)
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
        call molsnew(job)%append(ich)
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
  !$omp end single
  !$omp end parallel

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

!>--- close files (if they are open)
  if (dump) then
    close (ich)
    close (ich2)
  end if

  deallocate (grads)
  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  if (allocated(molsnew)) deallocate (molsnew)
  return
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
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  use crest_restartlog,only:trackrestart,restart_write_dummy
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mol
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex,nested,use_tmd_threads
  logical :: all_gfnff,has_active_calc,gfnff_thread_cap
  integer :: T,Tn,Trestore,Tnrestore,thread_save
  real(wp) :: percent
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  real(wp) :: etmp
  real(wp),allocatable :: grdtmp(:,:)
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
  gfnff_thread_cap = all_gfnff.and.nsim > 0.and.env%Threads > nsim
  if (gfnff_thread_cap) env%Threads = nsim
  call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  nested = env%omp_allow_nested
  if (all_gfnff) then
    Tn = 1
    write (stdout,'(1x,a,i0,a,i0)') &
    & 'GFN-FF MTD scheduling: parallel trajectories=',T,', cores per trajectory=',1
  end if

  allocate (calculations(T),source=env%calc)
  allocate (moltmps(T),source=mol)
  allocate (grdtmp(3,mol%nat),source=0.0_wp)
  if (all_gfnff) call ompmklset(1)
  do i = 1,T
    moltmps(i)%nat = mol%nat
    moltmps(i)%at = mol%at
    moltmps(i)%xyz = mol%xyz
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
    !>--- initialize the calculations
    call engrad(moltmps(i),calculations(i),etmp,grdtmp,io)
  end do
  if (all_gfnff) call ompmklset(T)

  !>--- other settings
  pr = .false.
  call profiler%init(nsim)

  !>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mol,pr,percent,ich, nsim, moltmps, nested,Tn) &
  !!$omp single
  !$omp private(vz,i,job,thread_id,io,ex)
  !$omp do
  do i = 1,nsim

    call initsignal()
    vz = i

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

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
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)
    mddats(vz)%termination_status = io
    call profiler%stop(vz)

    !>--- finish printout (thread safe)
    call parallel_md_finish_printout(mddats(vz),vz,io,profiler)
    !!$omp end task
  end do
  !!$omp taskwait
  !$omp end parallel

  !>--- collect trajectories into one
  call collect(nsim,mddats)

  call profiler%clear()
  deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  if (use_tmd_threads.or.gfnff_thread_cap) then
    env%Threads = thread_save
    call new_ompautoset(env,'max',0,Trestore,Tnrestore)
  end if
  return
contains
  subroutine collect(n,mddats)
    implicit none
    integer :: n
    type(mddata) :: mddats(n)
    logical :: ex
    integer :: i,io,ich,ich2
    character(len=:),allocatable :: atmp
    character(len=256) :: btmp
    open (newunit=ich,file='crest_dynamics.trj')
    do i = 1,n
      atmp = mddats(i)%trajectoryfile
      inquire (file=atmp,exist=ex)
      if (ex) then
        open (newunit=ich2,file=atmp)
        io = 0
        do while (io == 0)
          read (ich2,'(a)',iostat=io) btmp
          if (io == 0) then
            write (ich,'(a)') trim(btmp)
          end if
        end do
        close (ich2)
      end if
    end do
    close (ich)
    return
  end subroutine collect
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
subroutine crest_search_multimd2(env,mols,mddats,nsim)
!*******************************************************************
!* subroutine crest_search_multimd2
!* this runs #nsim MDs on #nsim selected different structures (mols)
!*******************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  use crest_restartlog,only:trackrestart,restart_write_dummy
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mols(nsim)
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex,nested,use_tmd_threads
  logical :: all_gfnff,has_active_calc,gfnff_thread_cap
  integer :: T,Tn,Trestore,Tnrestore,thread_save
  real(wp) :: percent
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
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
  gfnff_thread_cap = all_gfnff.and.nsim > 0.and.env%Threads > nsim
  if (gfnff_thread_cap) env%Threads = nsim
  call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  nested = env%omp_allow_nested
  if (all_gfnff) then
    Tn = 1
    write (stdout,'(1x,a,i0,a,i0)') &
    & 'GFN-FF MTD scheduling: parallel trajectories=',T,', cores per trajectory=',1
  end if

  allocate (calculations(T),source=env%calc)
  allocate (moltmps(T),source=mols(1))
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
  end do

!>--- other settings
  pr = .false.
  call profiler%init(nsim)

!>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mols,pr,percent,ich, moltmps,profiler, nested,Tn)
  !$omp single
  do i = 1,nsim

    call initsignal()
    vz = i

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

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

!>--- collect trajectories into one
  call collect(nsim,mddats)

  call profiler%clear()
  deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  if (use_tmd_threads.or.gfnff_thread_cap) then
    env%Threads = thread_save
    call new_ompautoset(env,'max',0,Trestore,Tnrestore)
  end if
  return
contains
  subroutine collect(n,mddats)
    implicit none
    integer :: n
    type(mddata) :: mddats(n)
    logical :: ex
    integer :: i,io,ich,ich2
    character(len=:),allocatable :: atmp
    character(len=256) :: btmp
    open (newunit=ich,file='crest_dynamics.trj')
    do i = 1,n
      atmp = mddats(i)%trajectoryfile
      inquire (file=atmp,exist=ex)
      if (ex) then
        open (newunit=ich2,file=atmp)
        io = 0
        do while (io == 0)
          read (ich2,'(a)',iostat=io) btmp
          if (io == 0) then
            write (ich,'(a)') trim(btmp)
          end if
        end do
        close (ich2)
      end if
    end do
    close (ich)
    return
  end subroutine collect
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
