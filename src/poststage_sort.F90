!================================================================================!
! Exact array handoff for the multilevel post-optimizer CREGEN path.
!================================================================================!
module crest_poststage_sort
  use crest_parameters,only:wp,stdout
  use crest_data,only:systemdata,status_normal,status_failed,crefile
  use crest_poststage_ensemble,only:poststage_ensemble
  use cregen_interface,only:newcregen
  use crest_restartlog,only:restart_write_dummy
  use strucrd,only:rdensembleparam,rdensemble,grepenergy
  use iomod,only:remove
  use utilities,only:checkname_xyz
  use omp_lib,only:omp_get_wtime
  implicit none
  private
  public :: sort_and_check_poststage

contains

  subroutine sort_and_check_poststage(env,filename,input_buffer,output_buffer)
    type(systemdata),intent(inout) :: env
    character(len=*),intent(in) :: filename
    type(poststage_ensemble),intent(inout) :: input_buffer
    type(poststage_ensemble),intent(inout) :: output_buffer

    real(wp),parameter :: increase = 1.5_wp
    real(wp),parameter :: nthr = 0.05_wp
    real(wp) :: ewin,clock_start
    integer :: nallin,nallout,nallthr,T,Tn
    character(len=80) :: inpnam,outnam
    logical :: produced,input_exists

    call output_buffer%clear()
    if (.not.input_buffer%valid()) then
      write(stdout,*) '**ERROR** invalid poststage buffer before CREGEN'
      env%iostatus_meta = status_failed
      return
    end if
    nallin = input_buffer%nall
    env%nat = input_buffer%nat
    ewin = env%ewin
    call checkname_xyz(crefile,inpnam,outnam)
    if (trim(inpnam) /= trim(filename)) then
      write(stdout,*) '**ERROR** poststage CREGEN input-name mismatch: ', &
      & trim(filename),' /= ',trim(inpnam)
      env%iostatus_meta = status_failed
      return
    end if
    inquire(file=trim(filename),exist=input_exists)
    if (.not.input_exists) then
      write(stdout,*) '**ERROR** missing canonical optimizer artifact ',trim(filename)
      env%iostatus_meta = status_failed
      return
    end if

    clock_start = omp_get_wtime()
    call new_ompautoset(env,'max',0,T,Tn)
    call newcregen(env,0,input_buffer=input_buffer,output_buffer=output_buffer, &
    & memory_produced=produced,preserve_input=.true.)
    if (env%iostatus_meta /= status_normal) return
    call restart_write_dummy(trim(outnam))
    if (.not.produced) then
      call read_poststage_file(trim(outnam),output_buffer,env)
      if (env%iostatus_meta /= status_normal) return
    end if
    nallout = output_buffer%nall
    write(stdout,'(1x,i0,'' structures remain within '',f8.2,'' kcal/mol window'')') &
    & nallout,ewin

    if (.not.env%entropic) then
      nallthr = nint(float(nallin)*nthr)
      if (nallout < nallthr) then
        write(stdout,'(1x,''This is less than '',i0,''% of the initial '',i0,'' structures.'')') &
        & nint(nthr*100.0_wp),nallin
        write(stdout,'(1x,''Increasing energy window to include more...'')')
        call remove(trim(outnam))
        call output_buffer%clear()
        env%ewin = ewin*increase
        call new_ompautoset(env,'max',0,T,Tn)
        ! Consume the preserved full optimizer ensemble directly.  This avoids
        ! rereading and reconstructing the multi-gigabyte ordered artifact.
        call newcregen(env,0,input_buffer=input_buffer,output_buffer=output_buffer, &
        & memory_produced=produced)
        if (env%iostatus_meta /= status_normal) return
        call restart_write_dummy(trim(outnam))
        if (.not.produced) then
          call read_poststage_file(trim(outnam),output_buffer,env)
          if (env%iostatus_meta /= status_normal) return
        end if
        nallout = output_buffer%nall
        write(stdout,'(1x,i0,'' structures remain within '',f8.2,'' kcal/mol window'')') &
        & nallout,ewin*increase
      end if
    end if
    call input_buffer%clear()
    write(stdout,'(1x,a,f12.3,a)') 'Poststage CREGEN total wall time: ', &
    & omp_get_wtime()-clock_start,' sec'
  end subroutine sort_and_check_poststage

  subroutine read_poststage_file(filename,buffer,env)
    character(len=*),intent(in) :: filename
    type(poststage_ensemble),intent(inout) :: buffer
    type(systemdata),intent(inout) :: env
    integer :: nat,nall,i,io,unit
    logical :: exists

    call buffer%clear()
    inquire(file=filename,exist=exists)
    if (.not.exists) then
      write(stdout,*) '**ERROR** missing CREGEN output ',trim(filename)
      env%iostatus_meta = status_failed
      return
    end if
    open(newunit=unit,file=filename,status='old',action='read',iostat=io)
    if (io /= 0) then
      env%iostatus_meta = status_failed
      return
    end if
    read(unit,*,iostat=io) nat
    close(unit)
    if (io /= 0) then
      write(stdout,*) '**ERROR** empty CREGEN output ',trim(filename)
      env%iostatus_meta = status_failed
      return
    end if
    call rdensembleparam(filename,nat,nall)
    if (nat < 1 .or. nall < 1) then
      env%iostatus_meta = status_failed
      return
    end if
    allocate(buffer%at(nat),buffer%xyz(3,nat,nall),buffer%comments(nall), &
    & buffer%eread(nall),stat=io)
    if (io /= 0) then
      call buffer%clear()
      env%iostatus_meta = status_failed
      return
    end if
    call rdensemble(filename,nat,nall,buffer%at,buffer%xyz,buffer%comments)
    do i = 1,nall
      buffer%eread(i) = grepenergy(buffer%comments(i))
    end do
    buffer%nat = nat
    buffer%nall = nall
  end subroutine read_poststage_file

end module crest_poststage_sort
