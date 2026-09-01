!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

subroutine crest_search_imtdgc(env,tim)
!*******************************************************************
!* This is the re-implementation of CREST's iMTD-GC default workflow
!* 
!* Compared to the legacy implementation, this version
!* is separated from the entropy algo to keep things clean
!* The entropy algo (sMTD-iMTD) can be found in search_entropy.f90
!*******************************************************************
  use crest_parameters, only: wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod
  use utilities
  use cregen_interface
  use crest_multilevel_interface,only:crest_multilevel_oloop
  use parallel_interface,only:crest_search_multimd2
  use crest_poststage_ensemble,only:poststage_ensemble
  use nci_input_ensemble,only:load_nci_input_ensemble,replace_and_validate_nci_wall, &
  & write_nci_input_starts,write_nci_mtd_jobs
  use external_rerank_restart,only:external_rerank_state,external_rerank_checkpoint_exists, &
  & read_external_rerank_checkpoint,restore_external_rerank_state, &
  & write_external_rerank_checkpoint,validate_external_rerank_seed, &
  & remove_external_rerank_checkpoint
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  type(coord),allocatable :: nci_input_mols(:),mtd_mols(:)
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer,allocatable :: bias_indices(:),input_indices(:)
  integer :: nsim,nallout,ninputs,nbias,nbias_check
  integer :: iinput,ibias,ijob

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,ex
  character(len=80) :: atmp,btmp,str
  logical :: multilevel(6)
  logical :: start,lower
  type(external_rerank_state) :: erstate
  type(poststage_ensemble) :: mtd_poststage
  logical :: special_external_restart
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,13x,"│")') "CREST iMTD-GC SAMPLING"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!>--- saftey termination
  if(mol%nat .le. 2)then
     call catchdiatomic(env)
    return
  endif

!>--- Interpret a positional multi-frame XYZ only for the first NCI MTD stage.
  ninputs=1
  if (env%NCI) then
    call load_nci_input_ensemble(env,mol,nci_input_mols,ninputs)
    if (.not.allocated(nci_input_mols)) allocate(nci_input_mols(1),source=mol)
    call replace_and_validate_nci_wall(env,nci_input_mols)
    if (ninputs>1) call write_nci_input_starts(nci_input_mols)
  end if

!>--- sets the MD length according to a flexibility measure
  call md_length_setup(env) 
!>--- create the MD calculator saved to env
  call env_to_mddat(env)

!>--- staged external-rerank continuation or initial checkpoint validation
  special_external_restart=env%external_rerank_restart
  if(special_external_restart) then
    call read_external_rerank_checkpoint(erstate)
    call restore_external_rerank_state(env,erstate)
    call validate_external_rerank_seed(env,erstate)
  else if(env%external_rerank) then
    if(.not.allocated(env%input_settings_file)) then
      error stop '**ERROR** external reranking requires a TOML input file'
    end if
    if(env%Maxrestart<2) error stop '**ERROR** external reranking requires at least two MTD iterations'
    if(external_rerank_checkpoint_exists()) then
      error stop '**ERROR** crest.external_rerank already exists; use --restart or remove it'
    end if
  end if

  if (env%performMTD.and..not.special_external_restart) then
!>--- (optional) calculate a short 1ps test MTD to check settings
   call tim%start(1,'Trial metadynamics (MTD)')
   call trialmd(env)  
   call tim%stop(1)
   if(env%iostatus_meta .ne. 0 ) return
  end if

!===========================================================!
!>--- Start mainloop 
  env%nreset = 0
  start = .true.
  if(special_external_restart) then
    start=.false.
    env%nreset=0
  end if
  MAINLOOP : do
    call printiter
    if (.not. start) then
!>--- clean working files but preserve .cre_*.xyz iteration archives
      call clean_V2i
      if(.not.special_external_restart) env%nreset=env%nreset+1
    else 
!>--- at the beginning, wipe directory clean
      call V2cleanup(.false.)
    end if
!===========================================================!
!>--- Meta-dynamics loop
  mtdloop: do i = 1,env%Maxrestart

    call mtd_poststage%clear()

    if(special_external_restart.and.i<=erstate%completed_mtd) cycle mtdloop

    write(stdout,*)
    write(stdout,'(1x,a)') '------------------------------'
    write(stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
    write(stdout,'(1x,a)') '------------------------------'

    !> Always seed this MTD iteration from the currently selected reference.
    call env%ref%to(mol)
    if (start .and. i==1 .and. ninputs>1) then
      nbias=-1
      call crest_search_multimd_init(env,nci_input_mols(1),mddat,nbias)
      nsim=ninputs*nbias
      allocate(mddats(nsim),mtd_mols(nsim))
      allocate(bias_indices(nsim),input_indices(nsim))
      ijob=0
      do iinput=1,ninputs
        if (iinput>1) then
          nbias_check=nbias
          call crest_search_multimd_init(env,nci_input_mols(iinput),mddat,nbias_check)
          if (nbias_check/=nbias) error stop 'NCI multi-input bias count changed.'
          mddat%simtype=type_mtd
        end if
        do ibias=1,nbias
          ijob=ijob+1
          mtd_mols(ijob)=nci_input_mols(iinput)
          mddats(ijob)=mddat
          input_indices(ijob)=iinput
          bias_indices(ijob)=ibias
        end do
      end do
      if (ijob/=nsim) error stop 'NCI multi-input job construction failed.'
      call crest_search_multimd_init2_mapped(env,mddats,nsim,bias_indices,input_indices)
      call write_nci_mtd_jobs(mddats)
      write(stdout,'(1x,a,i0,a,i0,a,i0,a)') 'NCI first MTD batch: ',ninputs, &
        ' inputs x ',nbias,' biases = ',nsim,' trajectories'
      call tim%start(2,'Metadynamics (MTD)')
      call crest_search_multimd2(env,mtd_mols,mddats,nsim, &
      & poststage_out=mtd_poststage)
      call tim%stop(2)
      call write_nci_mtd_jobs(mddats)
      deallocate(mtd_mols,bias_indices,input_indices)
    else
      nsim=-1
      call crest_search_multimd_init(env,mol,mddat,nsim)
      allocate(mddats(nsim),source=mddat)
      call crest_search_multimd_init2(env,mddats,nsim)
      call tim%start(2,'Metadynamics (MTD)')
      call crest_search_multimd(env,mol,mddats,nsim)
      call tim%stop(2)
    end if
    if (env%iostatus_meta /= 0) then
      if (allocated(mddats)) deallocate(mddats)
      return
    end if
!>--- a file called crest_dynamics.trj should have been written
    ensnam = 'crest_dynamics.trj'
!>--- deallocate for next iteration
    if(allocated(mddats))deallocate(mddats)

!==========================================================!
!>--- Reoptimization of trajectories
    call tim%start(3,'Geometry optimization')
    call optlev_to_multilev(env%optlev,multilevel)
    if (mtd_poststage%valid()) then
      call crest_multilevel_oloop(env,ensnam,multilevel,input_buffer=mtd_poststage)
    else
      call crest_multilevel_oloop(env,ensnam,multilevel)
    end if
    call tim%stop(3)
    if(env%iostatus_meta .ne. 0 ) return

!>--- save the CRE under a backup name
    call checkname_xyz(crefile,atmp,str)
    call checkname_xyz('.cre',str,btmp)
    call rename(atmp,btmp)
!>--- save cregen output
    call checkname_tmp('cregen',atmp,btmp)
    call rename('cregen.out.tmp',btmp)

!=========================================================!
!>--- cleanup after first iteration and prepare next
    if (i .eq. 1 .and. start) then
      start = .false.
!>-- obtain a first lowest energy as reference
      env%eprivious = env%elowest
!>-- select the MTD bias count for the next iteration
      if (.not. env%readbias .and.  env%runver .ne. 33 .and. &
      &   env%runver .ne. 787878 ) then
        if (env%NCI) then
          if (env%nmetadyn /= 6) then
            error stop '**ERROR** standard NCI first iteration did not contain 6 MTD biases'
          end if
          select case (env%nci_next_iteration_trajectories)
          case (4)
            env%nmetadyn = 4
            write(stdout,'(1x,a)') 'NCI next-iteration trajectories: 4 (reduced from 6).'
          case (6)
            write(stdout,'(1x,a)') 'NCI next-iteration trajectories: 6 (full standard bias set).'
          case default
            error stop '**ERROR** nci_next_iteration_trajectories must be 4 or 6'
          end select
        else
          !> Legacy behavior for non-NCI conformer searches.
          env%nmetadyn = env%nmetadyn - 2
        end if
      end if
!>-- the cleanup 
      call clean_V2i
      if(env%external_rerank) then
        call write_external_rerank_checkpoint(env,trim(str),i)
        env%keepModef = .true. ! preserve iteration archives until explicit restart
        write(stdout,'(/,1x,a)') 'External-rerank checkpoint reached after iteration 1.'
        write(stdout,'(1x,a)') 'Rank the source conformers externally, copy the selected original finite geometry to'
        write(stdout,'(3x,a)') 'crest-best-external.xyz'
        write(stdout,'(1x,a)') 'and continue with:'
        write(stdout,'(3x,a)') 'crest crest-best-external.xyz --restart'
        return
      end if
!>-- and always do two cycles of MTDs
      cycle mtdloop 
    endif
!=========================================================!
!>--- Check for lowest energy
    call elowcheck(lower,env)
    if (.not. lower) then
      exit mtdloop
    end if
  enddo mtdloop
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
  write(stdout,*)
  write (stdout,'(''========================================'')')
  write (stdout,'(''           MTD Simulations done         '')')
  write (stdout,'(''========================================'')')
  write (stdout,'(1x,''Collecting ensmbles.'')')
!>-- collecting all ensembles saved as ".cre_*.xyz"
  call collectcre(env)                      
  call newcregen(env,0)
  call checkname_xyz(crefile,atmp,btmp)
!>--- remaining number of structures
  call remaining_in(atmp,env%ewin,nallout) 

!=========================================================!
!>--- (optional) Perform additional MDs on the lowest conformers
  if (env%rotamermds) then
    call tim%start(4,'Molecular dynamics (MD)')
    call crest_rotamermds(env,conformerfile)
    call tim%stop(4)
    if(env%iostatus_meta .ne. 0 ) return

!>--- Reoptimization of trajectories
    call checkname_xyz(crefile,atmp,btmp)
    write(stdout,'('' Appending file '',a,'' with new structures'')')trim(atmp)
    ensnam = 'crest_dynamics.trj'
    call appendto(ensnam,trim(atmp))
    call tim%start(3,'Geometry optimization')
    call crest_multilevel_wrap(env,trim(atmp),-1)
    call tim%stop(3)
    if(env%iostatus_meta .ne. 0 ) return

    call elowcheck(lower,env)
    if (lower) then
      call checkname_xyz(crefile,atmp,str)
      call checkname_xyz('.cre',str,btmp)
      call rename(atmp,btmp)
      cycle MAINLOOP
    end if
  end if

!=========================================================!
!>--- (optional) Perform GC step
    if (env%performCross) then
      call tim%start(5,'Genetic crossing (GC)')
      call crest_newcross3(env)
      call tim%stop(5)
      if(env%iostatus_meta .ne. 0 ) return

      call confg_chk3(env)
      call elowcheck(lower,env)
      if (lower) then
        call checkname_xyz(crefile,atmp,str)
        call checkname_xyz('.cre',str,btmp)
        call rename(atmp,btmp)
        if (env%iterativeV2) cycle MAINLOOP
      end if
    end if

!==========================================================!
!>--- exit mainloop
   exit MAINLOOP
  enddo MAINLOOP

!==========================================================!
!>--- final ensemble optimization
    write (stdout,'(/)')
    write (stdout,'(3x,''================================================'')')
    write (stdout,'(3x,''|           Final Geometry Optimization        |'')')
    write (stdout,'(3x,''================================================'')')
    call tim%start(3,'Geometry optimization')
    call checkname_xyz(crefile,atmp,str)
    call crest_multilevel_wrap(env,trim(atmp),0) 
    call tim%stop(3)                 
    if(env%iostatus_meta .ne. 0 ) return

!==========================================================!
!>--- final ensemble sorting
!  call newcregen(env,0) 
!> this is actually done within the last crest_multilevel_
!> call, so I comment it out here

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
    write (stdout,'(/)')
    call smallhead('Final Ensemble Information')
    call V2terminating()
    if(special_external_restart) call remove_external_rerank_checkpoint()

!==========================================================!
  return
end subroutine crest_search_imtdgc

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine crest_multilevel_wrap(env,ensnam,level)
!*************************************************
!* wrapper for the multilevel_oloop to select
!* only a single optimization level
!*************************************************
  use crest_parameters, only: wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use crest_multilevel_interface,only:crest_multilevel_oloop
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: ensnam
  integer,intent(in) :: level
  logical :: multilevel(6)
  integer :: k
  multilevel = .false.
  select case(level)
  case( 1: ) !> explicit selection (level is a positie integer)
    k = min(level,6) 
    multilevel(k) =.true.
  case default 
  !>-- map global variable to multilevel selection (level is 0 or negative)
    k = optlevmap_alt(env%optlev) + level
    k = max(1,k)
    multilevel(k) =.true.
  end select
  call crest_multilevel_oloop(env,ensnam,multilevel)
end subroutine crest_multilevel_wrap

!========================================================================================!
subroutine crest_multilevel_oloop(env,ensnam,multilevel_in,input_buffer)
!*******************************************************
!* multilevel optimization loop.
!* construct consecutive optimizations starting with
!* crude thresholds to very tight ones
!*******************************************************
  use crest_parameters, only: wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use utilities
  use crest_restartlog
  use parallel_interface
  use crest_poststage_ensemble,only:poststage_ensemble
  use crest_poststage_sort,only:sort_and_check_poststage
  use omp_lib,only:omp_get_wtime
  implicit none
  type(systemdata),intent(inout) :: env 
  character(len=*),intent(in) :: ensnam
  logical,intent(in) :: multilevel_in(6)
  type(poststage_ensemble),intent(inout),optional :: input_buffer
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,pr
  character(len=128) :: inpnam,outnam
  integer :: i,l,k,T,Tn
  real(wp) :: ewinbackup,rthrbackup
  real(wp) :: hlowbackup
  integer :: microbackup
  integer :: optlevelbackup
  logical :: multilevel(6)
  type(constraint),allocatable :: auto_nci_walls(:)
  type(poststage_ensemble) :: optimized_poststage,sorted_poststage
  real(wp) :: poststage_clock
  logical :: poststage_written,empty_refine_eligible

  interface 
    subroutine crest_refine(env,input,output)
      use crest_data
      implicit none
      type(systemdata),intent(inout) :: env
      character(len=*),intent(in) :: input
      character(len=*),intent(in),optional :: output
    end subroutine crest_refine
  end interface

!>--- save backup thresholds before any path can branch to the common cleanup
  ewinbackup     = env%ewin
  rthrbackup     = env%rthr
  optlevelbackup = env%calc%optlev
  hlowbackup     = env%calc%hlow_opt
  microbackup    = env%calc%micro_opt

!>--- NCI container walls are an MTD-only potential.  Remove only walls
!>    generated by the NCI setup and restore them at the single exit below.
  call env%calc%detach_auto_nci_walls(auto_nci_walls)
  if (env%NCI.and..not.env%legacy) then
    if (size(auto_nci_walls) /= 1) then
      write(stdout,'(1x,a,i0)') &
      & '**ERROR** expected one automatic NCI wall before optimization; found ', &
      & size(auto_nci_walls)
      env%iostatus_meta = status_failed
      goto 900
    end if
  end if
  if (size(auto_nci_walls) > 0) then
    write(stdout,'(1x,a,i0,a,i0)') 'Optimization NCI walls removed: ', &
    & size(auto_nci_walls),'; remaining calculator constraints: ', &
    & env%calc%nconstraints
  end if

!>--- set multilevels, or enforce just one
  multilevel(:) = .false.
  if(env%multilevelopt)then
     multilevel(:) = multilevel_in(:)    
  else
    k = optlevmap_alt(env%optlev)
    multilevel(k) = .true.
  endif

  pr = .false.
  l = count(multilevel)
  if( l > 1 )then
  pr = .true.
  write(stdout,*)
  write(stdout,'(1x,a)') '======================================'
  write(stdout,'(1x,a)') '|  Multilevel Ensemble Optimization  |'
  write(stdout,'(1x,a)') '======================================'
  endif
  
!>--- read ensemble or consume the exact concurrently parsed worker buffer
  poststage_clock = omp_get_wtime()
  if (present(input_buffer)) then
    if (.not.input_buffer%valid()) then
      write(stdout,*) '**ERROR** invalid in-memory poststage ensemble for ',trim(ensnam)
      env%iostatus_meta = status_failed
      goto 900
    end if
    nat = input_buffer%nat
    nall = input_buffer%nall
    call move_alloc(input_buffer%at,at)
    call move_alloc(input_buffer%xyz,xyz)
    call move_alloc(input_buffer%eread,eread)
    call input_buffer%clear()
    write(stdout,'(1x,a,i0,a)') 'Poststage in-memory input accepted: ',nall,' structures'
  else
    call rdensembleparam(ensnam,nat,nall)
    if (nall .lt. 1) then
      write(stdout,*) '**ERROR** empty ensemble file ',trim(ensnam)
      env%iostatus_meta = status_failed
      goto 900
    endif
    allocate (xyz(3,nat,nall),at(nat),eread(nall))
    call rdensemble(ensnam,nat,nall,at,xyz,eread)
  end if
  write(stdout,'(1x,a,f12.3,a)') 'Poststage initial ensemble ingest wall time: ', &
  & omp_get_wtime()-poststage_clock,' sec'
!>--- track ensemble for restart
  poststage_clock = omp_get_wtime()
  call trackensemble(ensnam,nat,nall,at,xyz,eread)
  write(stdout,'(1x,a,f12.3,a)') 'Poststage initial restart snapshot wall time: ', &
  & omp_get_wtime()-poststage_clock,' sec'
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

  write(stdout,'(1x,a,i0,a,a,a)')'Optimizing all ',nall, &
  & ' structures from file "',trim(ensnam),'" ...'

!>--- sequential optimizations of ensembles
  dump = .true. !> optimized structures will be written to crest_ensemble.xyz
  do i=1,6
    if(multilevel(i))then
     !>--- set threads
       call new_ompautoset(env,'auto',nall,T,Tn)
     !>--- set optimization parameters
       call set_multilevel_options(env,i,.true.)
     !>--- run parallel optimizations
       call optimized_poststage%clear()
       poststage_written = .false.
       empty_refine_eligible = .not.env%refine_presort.and. &
       & .not.allocated(env%refine_queue)
       if (empty_refine_eligible) then
         call crest_oloop(env,nat,nall,at,xyz,eread,dump, &
         & poststage_buffer=optimized_poststage, &
         & poststage_written=poststage_written)
       else
         call crest_oloop(env,nat,nall,at,xyz,eread,dump)
       end if
       if (env%iostatus_meta /= status_normal) then
         deallocate(eread,at,xyz)
         env%ewin           = ewinbackup
         env%rthr           = rthrbackup
         env%calc%optlev    = optlevelbackup
         env%calc%hlow_opt  = hlowbackup
         env%calc%micro_opt = microbackup
         goto 900
       end if
       deallocate(eread,at,xyz)
     !>--- rename ensemble and sort
       call checkname_xyz(crefile,inpnam,outnam)
       call rename(ensemblefile,trim(inpnam))
     !>--- check for empty ensemble content
       if (poststage_written) then
         if (.not.optimized_poststage%valid()) then
           write(stdout,*) '**ERROR** invalid canonical optimizer poststage buffer'
           env%iostatus_meta = status_failed
           goto 900
         end if
         nat = optimized_poststage%nat
         nall = optimized_poststage%nall
       else
         poststage_written = .false.
         call optimized_poststage%clear()
         call rdensembleparam(trim(inpnam),nat,nall)
       end if
       if (nall .lt. 1) then
         write(stdout,*) '**ERROR** empty ensemble file',trim(inpnam)
         env%iostatus_meta = status_failed
         goto 900
       endif

       write(stdout,*)
     !==========================================================!
     !>-- dedicated ensemble refinement step (overwrites inpnam)
      poststage_clock = omp_get_wtime()
      if (.not.poststage_written) call crest_refine(env,trim(inpnam))
      write(stdout,'(1x,a,f12.3,a)') 'Poststage refinement wall time: ', &
      & omp_get_wtime()-poststage_clock,' sec'
      if (env%iostatus_meta /= status_normal) then
        env%ewin           = ewinbackup
        env%rthr           = rthrbackup
        env%calc%optlev    = optlevelbackup
        env%calc%hlow_opt  = hlowbackup
        env%calc%micro_opt = microbackup
        goto 900
      end if
     !==========================================================!

     !>--- CREGEN sorting.  The empty-refinement fast path consumes the exact
     !>    canonical optimizer representation directly and returns the exact
     !>    canonical CREGEN writer representation.  All other modes retain the
     !>    legacy file-based path.
       if (poststage_written) then
         call sorted_poststage%clear()
         call sort_and_check_poststage(env,trim(inpnam),optimized_poststage, &
         & sorted_poststage)
       else
         call sort_and_check(env,trim(inpnam))
       end if
       call optimized_poststage%clear()
       if (env%iostatus_meta /= status_normal) goto 900
       call checkname_xyz(crefile,inpnam,outnam)
       if (poststage_written) then
         if (.not.sorted_poststage%valid()) then
           write(stdout,*) '**ERROR** invalid canonical CREGEN poststage buffer'
           env%iostatus_meta = status_failed
           goto 900
         end if
         nat = sorted_poststage%nat
         nall = sorted_poststage%nall
         call move_alloc(sorted_poststage%at,at)
         call move_alloc(sorted_poststage%xyz,xyz)
         call move_alloc(sorted_poststage%eread,eread)
         call sorted_poststage%clear()
       else
       !>--- check for empty ensemble content (again)
         call rdensembleparam(trim(inpnam),nat,nall)
         if (nall .lt. 1) then
           write(stdout,*) '**ERROR** empty ensemble file',trim(inpnam)
           env%iostatus_meta = status_failed
           goto 900
         endif
       !>--- read new ensemble for next iteration
         allocate (xyz(3,nat,nall),at(nat),eread(nall))
         call rdensemble(trim(inpnam),nat,nall,at,xyz,eread)
       end if
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
     !>--- Important: crest_oloop requires coordinates in Bohrs
       xyz = xyz / bohr
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
     !>--- restore default sorting thresholds
       env%ewin        = ewinbackup
       env%rthr        = rthrbackup
       env%calc%optlev = optlevelbackup
       env%calc%hlow_opt  = hlowbackup
       env%calc%micro_opt = microbackup
    endif
  enddo

900 continue
  call optimized_poststage%clear()
  call sorted_poststage%clear()
  if(allocated(eread)) deallocate(eread)
  if(allocated(at))  deallocate(at)
  if(allocated(xyz)) deallocate(xyz)
  env%ewin           = ewinbackup
  env%rthr           = rthrbackup
  env%calc%optlev    = optlevelbackup
  env%calc%hlow_opt  = hlowbackup
  env%calc%micro_opt = microbackup
  if (allocated(auto_nci_walls)) then
    call env%calc%attach_auto_nci_walls(auto_nci_walls)
    if (size(auto_nci_walls) > 0) then
      write(stdout,'(1x,a,i0,a,i0)') 'Optimization NCI walls restored: ', &
      & size(auto_nci_walls),'; total calculator constraints: ', &
      & env%calc%nconstraints
    end if
    deallocate(auto_nci_walls)
  end if
  return
contains
  subroutine set_multilevel_options(env,i,pr)
    implicit none
    type(systemdata) :: env
    integer,intent(in) :: i
    logical,intent(in) :: pr 

    env%calc%hlow_opt  = env%hlowopt
    env%calc%micro_opt = nint(env%microopt)

    select case( i )  
    case( 1 )
     if(pr) call smallhead('crude pre-optimization')
     env%calc%optlev =  -3
     !> larger thresholds
     env%rthr = env%rthr * 2.0d0
     env%ewin = aint(env%ewin * 2.0d0)
    case( 2 )
     if(pr) call smallhead('optimization with very loose thresholds')
     env%calc%optlev =  -2
     env%rthr = env%rthr *1.5d0
     env%ewin = aint(env%ewin * 2.0d0)
    case( 3 )
     if(pr) call smallhead('optimization with loose thresholds')
     env%calc%optlev =  -1
      env%ewin = aint(env%ewin*(10.0d0/6.0d0))
    case( 4 )
     if(pr) call smallhead('optimization with regular thresholds')
     env%calc%optlev =  0
    case( 5 )
     if(pr) call smallhead('optimization with tight thresholds')
     env%calc%optlev =  1
    case( 6 )
     if(pr) call smallhead('optimization with very tight thresholds')
     env%calc%optlev =  2
    case default
     if(pr) call smallhead('optimization with default thresholds')
     env%ewin        = 6.0_wp
     env%rthr        = 0.125_wp
     env%calc%optlev = 0
    end select

    call print_opt_data(env%calc, stdout)

  end subroutine set_multilevel_options
end subroutine crest_multilevel_oloop

!========================================================================================!
subroutine crest_rotamermds(env,ensnam)
!***********************************************************
!* set up and perform several MDs at different temperatures
!* on the lowest few conformers
!***********************************************************
  use crest_parameters, only: wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use parallel_interface,only:crest_search_multimd2
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: ensnam

  integer :: nsim
  type(mddata) :: mddat
  type(mddata),allocatable :: mddats(:)
  type(coord) :: mol 
  type(coord),allocatable :: mols(:)
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  integer :: nstrucs,i,j,k,io
  real(wp) :: temp,newtemp
  character(len=80) :: atmp
  
!>--- coord setup
  call env%ref%to(mol)
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write(stdout,*) '**ERROR** empty ensemble file',trim(ensnam)
    env%iostatus_meta = status_failed
    return
  endif

!>--- determine how many MDs need to be run and setup
  call adjustnormmd(env)
  nstrucs = min(nall, env%nrotammds)
  nsim = nstrucs * env%temps 
  call crest_search_multimd_init(env,mol,mddat,nsim)
  allocate (mddats(nsim), source=mddat)
  call crest_search_multimd_init2(env,mddats,nsim)
!>--- adjust T's and runtimes
  k = 0
  do i=1,env%temps
    !> each T block 100K higher
    temp = env%nmdtemp + (i-1)*100.0_wp
    do j=1,nstrucs
      k= k + 1
      mddats(k)%tsoll = temp
      !> reduce runtime by 50% compared to MTDs
      mddats(k)%length_ps = mddats(k)%length_ps * 0.5_wp
      call mdautoset(mddats(k),io)
    enddo
  enddo 

!>--- read ensemble and prepare mols
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: mols must be in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  allocate(mols(nsim), source=mol) 
  k = 0
  do i=1,env%temps
    do j=1,nstrucs
      k = k + 1
      mols(k)%at = at
      mols(k)%xyz(:,:) = xyz(:,:,j)
    enddo
  enddo
  deallocate(eread,at,xyz)
  
!>--- print what we are doing
  write(stdout,*)
  write(atmp,'(''Additional regular MDs on lowest '',i0,'' conformer(s)'')')nstrucs
  call smallheadline(trim(atmp)) 

!>--- and finally, run the MDs
  call crest_search_multimd2(env,mols,mddats,nsim)

  if(allocated(mols))deallocate(mols)
  if(allocated(mddats))deallocate(mddats) 
  return
end subroutine crest_rotamermds

!========================================================================================!
subroutine crest_newcross3(env)
!*********************************************************
!* wrapper for the conformational crossing
!* takes the latest crest_rotamers_*, crosses structures
!* and writes the optimized ones back to the file
!*********************************************************
  use crest_parameters
  use crest_data
  use iomod
  use utilities
  use crest_multilevel_interface,only:crest_multilevel_oloop
  implicit none
  type(systemdata) :: env  
  real(wp) :: ewinbackup
  integer  :: i,imax,tmpconf,nremain
  character(len=128) :: inpnam,outnam,refnam
  character(len=512) :: thispath,tmppath
  logical :: multilevel(6)
  real(wp),allocatable :: backupthr(:)

  multilevel = .false.
  call getcwd(thispath)

  do i = 1,1  !>-- technically it would be possible to repeat the crossing
!>-- determine max number of new structures
    imax = min(nint(env%mdtime*50.0d0),5000)
    if (env%setgcmax) then
      imax = nint(env%gcmax)
    else if(imax<0)then
      imax=5000
    end if
    if (env%quick) then
      imax = nint(float(imax)*0.5d0)
    end if

!>-- call the crossing routine
    call checkname_xyz(crefile,refnam,tmppath)
    call touch(trim(tmppath)) 
    call crest_crossing(env,imax,trim(refnam),env%gcmaxparent)
    if (imax .lt. 1) then
      call remove(trim(tmppath))
      return
      exit
    end if

!>-- optimize ensemble
    if (env%gcmultiopt) then !>-- optionally split into two steps
      call optlev_to_multilev(env%optlev,multilevel)
    else
      multilevel(4) = .true.
    end if
    call crest_multilevel_oloop(env,'confcross.xyz',multilevel)
    if(env%iostatus_meta .ne. 0 ) return

!>-- append optimized crossed structures and original to a single file
    call checkname_xyz(crefile,inpnam,outnam)
    write(stdout,'(a,a)')'appending new structures to ',trim(refnam)
    call appendto(trim(inpnam),trim(refnam))
    do while(trim(inpnam).ne.trim(refnam))
      call remove(trim(inpnam))
      call checkname_xyz(crefile,inpnam,outnam)
    enddo
  end do
end subroutine crest_newcross3
