!================================================================================!
! Fail-closed process-isolated scheduler for prepared GFN-FF MTD state.
!
! Process scheduling is runtime N workers x K calculator threads per worker.
! Multi-MTD execution uses this scheduler by default; workers never parse the
! original CREST input.  The parent retains preparation, mapping, collection,
! optimization, and reranking.
!================================================================================!
module mtd_process_scheduler
  use iso_c_binding,only:c_char,c_int,c_null_char
  use iso_fortran_env,only:int8,int32,int64,iostat_end
  use ieee_arithmetic,only:ieee_is_finite
  use crest_parameters,only:wp,stdout
  use crest_data,only:systemdata,status_normal,status_failed,status_config
  use crest_calculator,only:calcdata,jobtype,engrad,engrad_total,prepare_gfnff_topology, &
  & request_gfnff_hbond_update
  use gfnff_api,only:gfnff_write_topology_restart
  use api_helpers,only:gfnff_fragment_ids
  use strucrd,only:coord,i2e
  use crest_poststage_ensemble,only:poststage_ensemble, &
  & read_canonical_trajectory_slice
  use dynamics_module,only:mddata,dynamics,type_mtd,cv_rmsd
  use iomod,only:makedir,directory_exist
  use omp_lib,only:omp_set_dynamic,omp_set_max_active_levels,omp_set_num_threads, &
  & omp_get_dynamic,omp_get_max_threads,omp_get_max_active_levels, &
  & omp_get_num_places,omp_get_num_threads,omp_get_thread_num,omp_get_place_num, &
  & omp_get_place_num_procs,omp_get_place_proc_ids,omp_get_wtime
  use mtd_process_capsule,only:mtd_capsule_meta,write_mtd_capsule,read_mtd_capsule, &
  & verify_mtd_capsule_roundtrip,compare_binary_files_exact,supported_process_nci_walls
  implicit none
  private

  ! Preserve the pre-existing optimizer-affinity arm only for its historical
  ! validated full-node process batch.  These are optimizer-reference constants,
  ! not limits on process-MTD worker count or threads per trajectory.
  integer,parameter :: optimizer_reference_workers = 48
  integer,parameter :: optimizer_reference_threads = 4
  integer,parameter :: optimizer_reference_cpus = 192
  integer,parameter :: affinity_report_limit = 16
  real(wp),parameter :: production_thread_oracle_atol = 1.0e-12_wp
  integer(int32),parameter :: result_version = 2_int32
  character(len=32),parameter :: result_magic = 'CREST_MTD_PROCESS_RESULT_V2     '
  character(len=*),parameter :: optin_name = 'CREST_EXPERIMENTAL_MTD_PROCESS_ISOLATION'
  character(len=*),parameter :: thread_oracle_optin_name = &
  & 'CREST_VALIDATE_PROCESS_MTD_THREADS'
  character(len=*),parameter :: worker_flag = '--crest-internal-mtd-worker-v1'
  integer,save :: process_batch_counter = 0
  integer,save :: process_session_parent_pid = 0
  character(len=:),allocatable,save :: process_session_directory
  character(len=:),allocatable,save :: process_session_topology

  public :: resolve_mtd_process_isolation
  public :: run_mtd_process_batch
  public :: mtd_process_worker_dispatch

  interface
    subroutine openblasset(threads)
      integer,intent(in) :: threads
    end subroutine openblasset
#ifdef WITH_OPENBLAS
    integer function openblas_get_num_threads()
    end function openblas_get_num_threads
#endif

    integer(c_int) function c_parent_pid() bind(C,name='crest_mtd_parent_pid')
      import :: c_int
    end function c_parent_pid

    integer(c_int) function c_current_cpu() bind(C,name='crest_mtd_current_cpu')
      import :: c_int
    end function c_current_cpu

    integer(c_int) function c_get_thread_affinity(cpu_capacity,cpu_count,cpus) &
    & bind(C,name='crest_mtd_get_thread_affinity')
      import :: c_int
      integer(c_int),value :: cpu_capacity
      integer(c_int),intent(out) :: cpu_count,cpus(*)
    end function c_get_thread_affinity

    integer(c_int) function c_set_optimizer_affinity_process_batch_ready(ready) &
    & bind(C,name='crest_optimizer_affinity_set_process_batch_ready')
      import :: c_int
      integer(c_int),value :: ready
    end function c_set_optimizer_affinity_process_batch_ready

    integer(c_int) function c_secure_directory(path) bind(C,name='crest_mtd_secure_directory')
      import :: c_char,c_int
      character(kind=c_char),intent(in) :: path(*)
    end function c_secure_directory

    integer(c_int) function c_validate_layout(worker_count,thread_count) &
    & bind(C,name='crest_mtd_validate_layout')
      import :: c_int
      integer(c_int),value :: worker_count,thread_count
    end function c_validate_layout

    integer(c_int) function c_spawn_worker(capsule,workdir,worker_index,worker_count, &
    & thread_count,parent_pid,pid_out) bind(C,name='crest_mtd_spawn_worker')
      import :: c_char,c_int
      character(kind=c_char),intent(in) :: capsule(*)
      character(kind=c_char),intent(in) :: workdir(*)
      integer(c_int),value :: worker_index,worker_count,thread_count,parent_pid
      integer(c_int),intent(out) :: pid_out
    end function c_spawn_worker

    integer(c_int) function c_prepare_worker(workdir,worker_index,worker_count, &
    & thread_count,parent_pid,cpu_list) bind(C,name='crest_mtd_prepare_worker')
      import :: c_char,c_int
      character(kind=c_char),intent(in) :: workdir(*)
      integer(c_int),value :: worker_index,worker_count,thread_count,parent_pid
      character(kind=c_char),intent(in) :: cpu_list(*)
    end function c_prepare_worker

    integer(c_int) function c_wait_worker(pid_value,exit_code,term_signal) &
    & bind(C,name='crest_mtd_wait_worker')
      import :: c_int
      integer(c_int),value :: pid_value
      integer(c_int),intent(out) :: exit_code,term_signal
    end function c_wait_worker

    integer(c_int) function c_kill_worker(pid_value) bind(C,name='crest_mtd_kill_worker')
      import :: c_int
      integer(c_int),value :: pid_value
    end function c_kill_worker

    integer(c_int) function c_reap_worker_bounded(pid_value,grace_ms,exit_code,term_signal) &
    & bind(C,name='crest_mtd_reap_worker_bounded')
      import :: c_int
      integer(c_int),value :: pid_value,grace_ms
      integer(c_int),intent(out) :: exit_code,term_signal
    end function c_reap_worker_bounded
  end interface

contains

  subroutine resolve_mtd_process_isolation(requested,status_out,message)
    logical,intent(out) :: requested
    integer,intent(out) :: status_out
    character(len=*),intent(out) :: message
    character(len=32) :: value
    integer :: length,status
    requested = .true.
    status_out = status_normal
    message = ''
    value = ''
    call get_environment_variable(optin_name,value,length=length,status=status,trim_name=.true.)
    if (status == 1) return
    if (status /= 0) then
      status_out = status_config
      message = optin_name//' is present but unreadable or longer than 32 bytes'
      return
    end if
    if (length /= 1 .or. value(1:1) /= '1') then
      status_out = status_config
      message = optin_name//' is deprecated; if present, the only accepted value is 1'
      return
    end if
  end subroutine resolve_mtd_process_isolation

  subroutine run_mtd_process_batch(env,mols,mddats,nsim,outer_threads,inner_threads, &
  & status,message,poststage_out)
    type(systemdata),intent(inout) :: env
    integer,intent(in) :: nsim,outer_threads,inner_threads
    type(coord),intent(in) :: mols(nsim)
    type(mddata),intent(inout) :: mddats(nsim)
    integer,intent(out) :: status
    character(len=*),intent(out) :: message
    type(poststage_ensemble),intent(inout),optional :: poststage_out
    type(mtd_capsule_meta) :: meta
    type(mddata) :: md_copy
    type(calcdata) :: capsule_calc
    integer(c_int),allocatable :: pids(:)
    integer(c_int) :: exit_code,term_signal,cstat
    logical,allocatable :: reaped(:)
    logical :: identical,topology_exists
    integer :: parent_pid,i,io,frames,expected_frames,parser_threads
    integer :: requested_parser_threads
    integer,allocatable :: expected_frame_counts(:)
    integer,allocatable :: slice_first(:),slice_last(:)
    integer,allocatable :: trajectory_status(:)
    integer(int64) :: trajectory_bytes,worker_engrad_calls,batch_engrad_calls
    integer(int64),allocatable :: trajectory_bytes_by_worker(:)
    integer(int64) :: total_frames64
    real(wp) :: worker_wall
    real(wp),allocatable :: worker_walls(:)
    real(wp) :: parser_start,parser_finish
    real(wp),allocatable :: trajectory_read_seconds(:)
    real(wp),allocatable :: trajectory_decode_seconds(:)
    character(len=4096) :: cwd_buffer
    character(len=64) :: parser_thread_setting
    character(len=1024) :: iomessage
    character(len=1024),allocatable :: trajectory_messages(:)
    character(len=:),allocatable :: cwd,auditroot,sessiondir,batchdir,workdir,capsule
    character(len=:),allocatable :: topology_master
    character(len=:),allocatable :: trajectory,restart_file,topology,result_file,calcdir
    character(kind=c_char),allocatable :: c_capsule(:),c_workdir(:)

    status = status_normal
    message = ''
    if (present(poststage_out)) call poststage_out%clear()
    allocate(pids(nsim),reaped(nsim),expected_frame_counts(nsim),slice_first(nsim), &
    & slice_last(nsim),trajectory_status(nsim),trajectory_bytes_by_worker(nsim), &
    & worker_walls(nsim),trajectory_read_seconds(nsim),trajectory_decode_seconds(nsim), &
    & trajectory_messages(nsim))
    pids = 0_c_int
    reaped = .false.
    cstat = c_set_optimizer_affinity_process_batch_ready(0_c_int)
    if (cstat /= 0_c_int) then
      status = status_config
      message = 'failed to clear optimizer-affinity process-batch state'
      return
    end if
    pids = 0_c_int
    reaped = .false.
    batch_engrad_calls = 0_int64
    ! shakedata%freezeptr has no default => null() initialization upstream.
    ! SHAKE is not serialized by this process capsule, so define its semantically-dead
    ! association status before any ASSOCIATED() check or capsule assignment.
    if (nsim > 0) then
      do i = 1,nsim
        nullify(mddats(i)%shk%freezeptr)
        ! Native dynamics() calls this before reading maxblock or any other
        ! fallback field.  Do the same before serialization; the child repeats
        ! this idempotently, so its resolved scientific settings are unchanged.
        call mddats(i)%defaults()
      end do
    end if
    call validate_process_configuration(env,mols,mddats,nsim,outer_threads, &
    & inner_threads,status,message)
    if (status /= status_normal) return

    cwd_buffer = ''
    call getcwd(cwd_buffer)
    cwd = trim(cwd_buffer)
    if (len_trim(cwd) == 0) then
      status = status_config
      message = 'cannot determine parent working directory for process MTD'
      return
    end if
    parent_pid = int(c_parent_pid())
    if (parent_pid < 2) then
      status = status_config
      message = 'invalid process ID for process-isolated MTD parent'
      return
    end if
    auditroot = trim(cwd)//'/crest_process_mtd_audit'
    if (.not.directory_exist(auditroot)) then
      io = makedir(auditroot)
      if (io /= 0) then
        status = status_config
        message = 'cannot create persistent process MTD audit root: '//trim(auditroot)
        return
      end if
    end if
    call secure_directory(auditroot,io,iomessage)
    if (io /= 0) then
      status = status_config
      message = 'process MTD audit root is not a private 0700 directory: '//trim(iomessage)
      return
    end if
    if (process_session_parent_pid == 0) then
      sessiondir = trim(auditroot)//'/process_session_'//integer_string(parent_pid)
      if (directory_exist(sessiondir)) then
        status = status_config
        message = 'refusing to reuse existing process MTD session directory: '//trim(sessiondir)
        return
      end if
      io = makedir(sessiondir)
      if (io /= 0) then
        status = status_config
        message = 'cannot create private process MTD session directory: '//trim(sessiondir)
        return
      end if
      call secure_directory(sessiondir,io,iomessage)
      if (io /= 0) then
        status = status_config
        message = 'process MTD session is not a private 0700 directory: '//trim(iomessage)
        return
      end if
      topology_master = trim(sessiondir)//'/canonical_gfnff_topology.input'
      call gfnff_write_topology_restart(mols(1)%nat,env%calc%calcs(1)%ff_dat, &
      & topology_master,io)
      if (io /= 0) then
        status = status_config
        message = 'cannot export initialized canonical GFN-FF topology; iostat='// &
        & integer_string(io)
        return
      end if
      inquire(file=topology_master,exist=topology_exists,iostat=io,iomsg=iomessage)
      if (io /= 0 .or. .not.topology_exists) then
        status = status_config
        message = 'exported canonical GFN-FF topology is absent: '//trim(iomessage)
        return
      end if
      process_session_parent_pid = parent_pid
      process_session_directory = sessiondir
      process_session_topology = topology_master
      process_batch_counter = 0
    else
      if (process_session_parent_pid /= parent_pid) then
        status = status_config
        message = 'process MTD parent process ID changed within one invocation'
        return
      end if
      if (.not.allocated(process_session_directory)) then
        status = status_config
        message = 'process MTD session directory state is unavailable'
        return
      end if
      if (.not.allocated(process_session_topology)) then
        status = status_config
        message = 'process MTD canonical topology state is unavailable'
        return
      end if
      sessiondir = process_session_directory
      topology_master = process_session_topology
      if (.not.directory_exist(sessiondir)) then
        status = status_config
        message = 'process MTD session directory disappeared: '//trim(sessiondir)
        return
      end if
      call secure_directory(sessiondir,io,iomessage)
      if (io /= 0) then
        status = status_config
        message = 'process MTD session is not a private 0700 directory: '//trim(iomessage)
        return
      end if
      inquire(file=topology_master,exist=topology_exists,iostat=io,iomsg=iomessage)
      if (io /= 0 .or. .not.topology_exists) then
        status = status_config
        message = 'canonical GFN-FF topology disappeared from process MTD session'
        return
      end if
    end if
    if (process_batch_counter == huge(process_batch_counter)) then
      status = status_config
      message = 'process MTD batch counter overflow'
      return
    end if
    process_batch_counter = process_batch_counter+1
    batchdir = trim(sessiondir)//'/process_batch_'//zero_pad4(process_batch_counter)
    if (directory_exist(batchdir)) then
      status = status_config
      message = 'refusing to reuse existing process MTD batch directory: '//trim(batchdir)
      return
    end if
    io = makedir(batchdir)
    if (io /= 0) then
      status = status_config
      message = 'cannot create private process MTD batch directory: '//trim(batchdir)
      return
    end if
    call secure_directory(batchdir,io,iomessage)
    if (io /= 0) then
      status = status_config
      message = 'process MTD batch is not a private 0700 directory: '//trim(iomessage)
      return
    end if
    ! Native threaded workers inherit the calculator state left by trial MTD,
    ! while an exec'd process cannot safely inherit allocatable Fortran/OpenMP
    ! workspaces.  Make that semantic boundary explicit with two distinct
    ! contracts.  First prove inherited-versus-fresh history semantics bitwise
    ! at deterministic T1.  Then prove a separately fresh T1 calculator versus
    ! the production TK path under the preregistered absolute 1e-12 envelope.
    ! Topology, fragment, freeze, frozen-gradient, and mutable HB/XB state remain
    ! exact in both contracts.  Only the proven-clean settings object is
    ! serialized; ff_dat is never serialized.
    call prepare_canonical_capsule_calculator(env%calc,mols,mddats,topology_master, &
    & inner_threads,capsule_calc,status,iomessage)
    if (status /= status_normal) then
      message = 'post-trial versus canonical GFN-FF oracle failed: '//trim(iomessage)
      return
    end if

    write(stdout,'(1x,a)') 'Process-isolated MTD scheduler enabled (capsule v5, runtime N x K).'
    write(stdout,'(1x,a,i0,a,i0)') 'Process MTD scheduling: workers=',nsim, &
    & ', GFN-FF threads per worker=',inner_threads
    write(stdout,'(1x,a)') 'Process MTD RNG policy: stochastic; samerand is forbidden.'

    ! Prepare every private input before starting any child.  A malformed or
    ! non-round-tripping capsule therefore fails before any dynamics executes.
    do i = 1,nsim
      workdir = trim(batchdir)//'/worker_'//zero_pad4(i)
      calcdir = trim(workdir)//'/CALC'
      io = makedir(workdir)
      if (io == 0) io = makedir(calcdir)
      if (io /= 0) then
        status = status_config
        message = 'cannot create private directory for MTD worker '//integer_string(i)
        return
      end if
      call secure_directory(workdir,io,iomessage)
      if (io == 0) call secure_directory(calcdir,io,iomessage)
      if (io /= 0) then
        status = status_config
        message = 'worker directory is not private 0700 for MTD worker '//integer_string(i)
        return
      end if
      capsule = trim(workdir)//'/input.capsule'
      topology = trim(workdir)//'/gfnff_topo'
      result_file = trim(workdir)//'/result.bin'
      trajectory = trim(workdir)//'/crest_'//integer_string(i)//'.trj'
      restart_file = trim(workdir)//'/crest_'//integer_string(i)//'.mdrestart'
      call copy_file_exact(topology_master,topology,io,iomessage)
      if (io /= 0) then
        status = status_config
        message = 'private topology copy failed for worker '//integer_string(i)//': '//trim(iomessage)
        return
      end if
      call compare_binary_files_exact(topology_master,topology,identical,io,iomessage)
      if (io /= 0 .or. .not.identical) then
        status = status_config
        message = 'private topology pre-spawn verification failed for worker '// &
        & integer_string(i)//': '//trim(iomessage)
        return
      end if

      md_copy = mddats(i)
      nullify(md_copy%shk%freezeptr)
      md_copy%trajectoryfile = trajectory
      md_copy%restartfile = restart_file
      md_copy%termination_status = -1
      meta%worker_index = i
      meta%parent_pid = parent_pid
      meta%inner_threads = inner_threads
      meta%expected_atoms = mols(i)%nat
      meta%workdir = workdir
      meta%topology_file = topology
      meta%result_file = result_file
      call write_mtd_capsule(capsule,meta,mols(i),md_copy,capsule_calc,io,iomessage)
      if (io /= 0) then
        status = status_config
        message = 'capsule write failed for worker '//integer_string(i)//': '//trim(iomessage)
        return
      end if
      call verify_mtd_capsule_roundtrip(capsule,meta,mols(i),md_copy,capsule_calc,io,iomessage)
      if (io /= 0) then
        status = status_config
        message = 'capsule roundtrip failed for worker '//integer_string(i)//': '//trim(iomessage)
        return
      end if
      mddats(i) = md_copy
      write(stdout,'(2x,a,i0,a,i0,a,i0)') 'prepared process MTD ',i, &
      & ' (input=',mddats(i)%input_structure_id,', bias=',mddats(i)%bias_configuration_id,')'
    end do

    do i = 1,nsim
      workdir = trim(batchdir)//'/worker_'//zero_pad4(i)
      capsule = trim(workdir)//'/input.capsule'
      call to_c_string(capsule,c_capsule)
      call to_c_string(workdir,c_workdir)
      cstat = c_spawn_worker(c_capsule,c_workdir,int(i,c_int),int(nsim,c_int), &
      & int(inner_threads,c_int),int(parent_pid,c_int),pids(i))
      if (cstat /= 0_c_int) then
        call terminate_and_reap_workers(pids,reaped,i-1)
        status = status_failed
        message = 'posix_spawn failed for MTD worker '//integer_string(i)// &
        & ' with errno '//integer_string(int(cstat))
        return
      end if
    end do

    do i = 1,nsim
      exit_code = -1_c_int
      term_signal = -1_c_int
      cstat = c_wait_worker(pids(i),exit_code,term_signal)
      reaped(i) = cstat == 0_c_int
      if (cstat /= 0_c_int .or. exit_code /= 0_c_int .or. term_signal /= 0_c_int) then
        call terminate_and_reap_workers(pids,reaped,nsim)
        status = status_failed
        write(message,'(a,i0,a,i0,a,i0,a,i0)') 'MTD worker ',i, &
        & ' failed: wait_errno=',int(cstat),', exit=',int(exit_code),', signal=',int(term_signal)
        return
      end if
    end do

    expected_frame_counts = 0
    slice_first = 0
    slice_last = 0
    trajectory_status = status_normal
    trajectory_bytes_by_worker = 0_int64
    trajectory_messages = ''
    trajectory_read_seconds = 0.0_wp
    trajectory_decode_seconds = 0.0_wp
    worker_walls = 0.0_wp
    total_frames64 = 0_int64
    do i = 1,nsim
      workdir = trim(batchdir)//'/worker_'//zero_pad4(i)
      result_file = trim(workdir)//'/result.bin'
      call read_worker_result(result_file,i,mddats(i)%termination_status,worker_wall, &
      & worker_engrad_calls,io,iomessage)
      if (io /= 0 .or. mddats(i)%termination_status /= 0) then
        status = status_failed
        message = 'invalid dynamics result for MTD worker '//integer_string(i)//': '//trim(iomessage)
        return
      end if
      if (worker_engrad_calls /= int(mddats(i)%length_steps,int64)) then
        status = status_failed
        message = 'unexpected energy/gradient call count for MTD worker '//integer_string(i)
        return
      end if
      if (worker_engrad_calls > huge(batch_engrad_calls)-batch_engrad_calls) then
        status = status_failed
        message = 'energy/gradient call accounting overflow for MTD worker '//integer_string(i)
        return
      end if
      batch_engrad_calls = batch_engrad_calls+worker_engrad_calls
      call scan_worker_logs(workdir,io,iomessage)
      if (io /= 0) then
        status = status_failed
        message = 'fatal-output audit failed for MTD worker '//integer_string(i)//': '//trim(iomessage)
        return
      end if
      topology = trim(workdir)//'/gfnff_topo'
      call compare_binary_files_exact(topology_master,topology,identical,io,iomessage)
      if (io /= 0 .or. .not.identical) then
        status = status_failed
        message = 'private topology changed or became unreadable for MTD worker '// &
        & integer_string(i)//': '//trim(iomessage)
        return
      end if
      if (mddats(i)%sdump < 1 .or. mddats(i)%length_steps < 2) then
        status = status_failed
        message = 'invalid trajectory dimensions for MTD worker '//integer_string(i)
        return
      end if
      expected_frames = (mddats(i)%length_steps-1)/mddats(i)%sdump
      if (expected_frames < 1 .or. &
      & int(expected_frames,int64) > huge(total_frames64)-total_frames64) then
        status = status_failed
        message = 'trajectory frame accounting failed for MTD worker '//integer_string(i)
        return
      end if
      expected_frame_counts(i) = expected_frames
      total_frames64 = total_frames64+int(expected_frames,int64)
      worker_walls(i) = worker_wall

      if (.not.present(poststage_out)) then
        call validate_xyz_trajectory(mddats(i)%trajectoryfile,mols(i)%at, &
        & expected_frames,frames,trajectory_bytes,io,iomessage)
        if (io /= 0) then
          status = status_failed
          message = 'trajectory validation failed for MTD worker '// &
          & integer_string(i)//': '//trim(iomessage)
          return
        end if
        trajectory_bytes_by_worker(i) = trajectory_bytes
        call print_process_finish(mddats(i),i,worker_wall,frames,trajectory_bytes)
      end if
    end do

    if (present(poststage_out)) then
      if (total_frames64 > int(huge(1),int64)) then
        status = status_failed
        message = 'poststage trajectory frame total exceeds default-integer indexing'
        return
      end if
      do i = 2,nsim
        if (size(mols(i)%at) /= size(mols(1)%at)) then
          status = status_failed
          message = 'poststage atom count differs for MTD worker '//integer_string(i)
          return
        end if
        if (any(mols(i)%at /= mols(1)%at)) then
          status = status_failed
          message = 'poststage atom order differs for MTD worker '//integer_string(i)
          return
        end if
      end do

      poststage_out%nat = size(mols(1)%at)
      poststage_out%nall = int(total_frames64)
      allocate(poststage_out%at(poststage_out%nat), &
      & poststage_out%xyz(3,poststage_out%nat,poststage_out%nall), &
      & poststage_out%eread(poststage_out%nall), &
      & poststage_out%comments(poststage_out%nall),stat=io,errmsg=iomessage)
      if (io /= 0) then
        call poststage_out%clear()
        status = status_failed
        message = 'poststage trajectory allocation failed: '//trim(iomessage)
        return
      end if
      poststage_out%at = mols(1)%at
      slice_first(1) = 1
      slice_last(1) = expected_frame_counts(1)
      do i = 2,nsim
        slice_first(i) = slice_last(i-1)+1
        slice_last(i) = slice_first(i)+expected_frame_counts(i)-1
      end do
      if (slice_last(nsim) /= poststage_out%nall) then
        call poststage_out%clear()
        status = status_failed
        message = 'poststage trajectory prefix-sum invariant failed'
        return
      end if

      parser_threads = min(nsim,max(1,omp_get_max_threads()))
      parser_thread_setting = ''
      call get_environment_variable('CREST_POSTSTAGE_PARSER_THREADS', &
      & value=parser_thread_setting,status=io)
      if (io == 0 .and. len_trim(parser_thread_setting) > 0) then
        read(parser_thread_setting,*,iostat=io) requested_parser_threads
        if (io /= 0 .or. requested_parser_threads < 1 .or. &
        & requested_parser_threads > nsim) then
          call poststage_out%clear()
          status = status_failed
          write(message,'(a,i0)') 'CREST_POSTSTAGE_PARSER_THREADS must be an integer in 1:',nsim
          return
        end if
        parser_threads = min(parser_threads,requested_parser_threads)
      else if (io < 0) then
        call poststage_out%clear()
        status = status_failed
        message = 'CREST_POSTSTAGE_PARSER_THREADS value is too long'
        return
      end if
      parser_start = omp_get_wtime()
!$omp parallel do default(none) schedule(static,1) num_threads(parser_threads) &
!$omp shared(mddats,mols,expected_frame_counts,slice_first,slice_last,poststage_out) &
!$omp shared(trajectory_bytes_by_worker,trajectory_status,trajectory_messages) &
!$omp shared(trajectory_read_seconds,trajectory_decode_seconds,nsim) &
!$omp private(i)
      do i = 1,nsim
        call read_canonical_trajectory_slice(mddats(i)%trajectoryfile,mols(i)%at, &
        & expected_frame_counts(i), &
        & poststage_out%xyz(:,:,slice_first(i):slice_last(i)), &
        & poststage_out%eread(slice_first(i):slice_last(i)), &
        & poststage_out%comments(slice_first(i):slice_last(i)), &
        & trajectory_bytes_by_worker(i),trajectory_status(i),trajectory_messages(i), &
        & trajectory_read_seconds(i),trajectory_decode_seconds(i))
      end do
!$omp end parallel do
      parser_finish = omp_get_wtime()
      write(stdout,'(1x,a,f12.3,a)') 'Poststage parallel trajectory parse wall time: ', &
      & parser_finish-parser_start,' sec'
      write(stdout,'(1x,a,i0)') 'Poststage parallel trajectory parse threads: ',parser_threads
      write(stdout,'(1x,a)') 'Poststage trajectory parser backend: buffered C strict parser v1'
      write(stdout,'(1x,a,f12.3,a)') 'Poststage maximum worker file-read time: ', &
      & maxval(trajectory_read_seconds),' sec'
      write(stdout,'(1x,a,f12.3,a)') 'Poststage maximum worker numeric-decode time: ', &
      & maxval(trajectory_decode_seconds),' sec'
      write(stdout,'(1x,a,f12.3,a)') 'Poststage aggregate worker file-read time: ', &
      & sum(trajectory_read_seconds),' worker-sec'
      write(stdout,'(1x,a,f12.3,a)') 'Poststage aggregate worker numeric-decode time: ', &
      & sum(trajectory_decode_seconds),' worker-sec'
      if (parser_finish > parser_start) then
        write(stdout,'(1x,a,f12.3,a)') 'Poststage aggregate parser throughput: ', &
        & real(sum(trajectory_bytes_by_worker),wp)/(parser_finish-parser_start)/ &
        & (1024.0_wp*1024.0_wp),' MiB/s'
      end if

      do i = 1,nsim
        if (trajectory_status(i) /= status_normal) then
          status = status_failed
          message = 'trajectory validation failed for MTD worker '// &
          & integer_string(i)//': '//trim(trajectory_messages(i))
          call poststage_out%clear()
          return
        end if
      end do
      if (.not.poststage_out%valid()) then
        status = status_failed
        message = 'parallel trajectory parser produced an invalid poststage ensemble'
        call poststage_out%clear()
        return
      end if
      do i = 1,nsim
        call print_process_finish(mddats(i),i,worker_walls(i), &
        & expected_frame_counts(i),trajectory_bytes_by_worker(i))
      end do
    end if
    if (batch_engrad_calls > huge(engrad_total)-engrad_total) then
      status = status_failed
      message = 'parent energy/gradient call accounting overflow after process MTD'
      if (present(poststage_out)) call poststage_out%clear()
      return
    end if
    ! Do not generalize optimizer behavior as a side effect of this scheduler
    ! redesign.  The optimizer-affinity feature was validated only after the
    ! historical 48 x 4 = 192-core process batch, so preserve that exact arm
    ! condition and leave it disabled after every other valid N x K batch.
    if (nsim == optimizer_reference_workers .and. &
    & inner_threads == optimizer_reference_threads .and. &
    & env%Threads == optimizer_reference_cpus) then
      cstat = c_set_optimizer_affinity_process_batch_ready(1_c_int)
      if (cstat /= 0_c_int) then
        status = status_failed
        message = 'failed to arm historical optimizer affinity after validated 48x4 process MTD batch'
        if (present(poststage_out)) call poststage_out%clear()
        return
      end if
    end if
    engrad_total = engrad_total+batch_engrad_calls
    write(stdout,'(1x,a,i0)') 'Process MTD energy+gradient calls added to parent: ', &
    & batch_engrad_calls
    write(stdout,'(1x,a,i0,a)') 'All ',nsim, &
    & ' process-isolated MTD workers exited 0 and passed trajectory validation.'
  end subroutine run_mtd_process_batch

  subroutine prepare_canonical_capsule_calculator(source_calc,mols,mddats,topology_file, &
  & production_threads,capsule_calc,status,message)
!*******************************************************************************
!* Make the process boundary explicit and prove that dropping process-private
!* GFN-FF runtime storage does not change the prepared calculator result.
!*
!* Native crest_search_multimd2 gives each outer worker an independent deep copy
!* of the untouched post-trial calculator.  Model that semantic boundary with a
!* deterministic T1 bitwise oracle on every first input and on an exactly
!* representable controlled-history displacement.  Do not ask dynamic OpenMP
!* reductions for impossible cross-thread-count byte identity: a separate fresh-T1/fresh-TK
!* production-thread oracle applies the preregistered absolute 1e-12 envelope to
!* energy, components, charges, and gradients only.  Topology, fragments, freeze
!* state, positive-zero frozen gradients, and live HB/XB state remain exact.
!* Neither contract changes a scientific kernel or its result tolerance.
!*******************************************************************************
    type(calcdata),intent(in) :: source_calc
    type(coord),intent(in) :: mols(:)
    type(mddata),intent(in) :: mddats(:)
    character(len=*),intent(in) :: topology_file
    integer,intent(in) :: production_threads
    type(calcdata),intent(out) :: capsule_calc
    integer,intent(out) :: status
    character(len=*),intent(out) :: message
#ifdef WITH_GFNFF
    type(calcdata) :: inherited_calc,fresh_calc,fresh_t1_calc,fresh_tk_calc
    type(coord) :: probe_mol
    real(wp),allocatable :: inherited_gradient(:,:),fresh_gradient(:,:)
    real(wp),allocatable :: fresh_t1_gradient(:,:),fresh_tk_gradient(:,:)
    real(wp) :: inherited_energy,fresh_energy,fresh_t1_energy,fresh_tk_energy
    real(wp) :: differences(4),max_production_differences(4)
    integer(int64) :: engrad_before
    integer :: probe_id,noracle,index,io,ninitialized,ngfnff,saved_threads,saved_levels
    integer :: nat,probe_atom,fragment_io
    integer,allocatable :: oracle_indices(:),expected_fragments(:)
#ifdef WITH_OPENBLAS
    integer :: saved_blas_threads
#endif
    logical :: same,saved_dynamic,run_thread_oracle
    character(len=32) :: thread_oracle_value
    integer :: thread_oracle_length,thread_oracle_status
#endif

    status = status_config
    message = ''
    capsule_calc = source_calc
    if (.not.allocated(capsule_calc%calcs)) then
      message = 'cannot construct canonical calculator from unallocated source settings'
      return
    end if
    if (size(capsule_calc%calcs) /= 1) then
      message = 'cannot construct canonical calculator from malformed source settings'
      return
    end if
    if (source_calc%id < 0 .or. source_calc%id > 1) then
      message = 'process MTD supports only the one-level weighted or direct calculator selector'
      return
    end if
    if (allocated(source_calc%eweight)) then
      if (size(source_calc%eweight) /= 1) then
        message = 'post-trial calculator weight cache has the wrong size'
        return
      end if
      if (source_calc%id == 0) then
        if (.not.same_real_bits_scalar(source_calc%eweight(1), &
        & source_calc%calcs(1)%weight)) then
          message = 'post-trial calculator weight cache differs from the configured level weight'
          return
        end if
      end if
    end if
    capsule_calc%calcs(1)%restart = .true.
    capsule_calc%calcs(1)%restartfile = trim(topology_file)

#ifndef WITH_GFNFF
    message = 'process-isolated MTD canonicalization requires compiled-in GFN-FF'
    return
#else
    if (size(mols) < 1) then
      message = 'calculator oracle received no prepared molecules'
      return
    end if
    nat = mols(1)%nat
    if (nat < 1) then
      message = 'calculator oracle received a nonpositive molecule atom count'
      return
    end if
    if (.not.allocated(source_calc%calcs(1)%ff_dat)) then
      message = 'required post-trial GFN-FF runtime state is unexpectedly unallocated'
      return
    end if
    if (.not.allocated(source_calc%calcs(1)%ff_dat%topo)) then
      message = 'post-trial GFN-FF state has no initialized topology'
      return
    end if
    if (.not.allocated(source_calc%calcs(1)%ff_dat%nlist)) then
      message = 'post-trial GFN-FF state has no initialized neighbour list'
      return
    end if
    if (allocated(source_calc%calcs(1)%gff_fragments)) then
      call gfnff_fragment_ids(source_calc%calcs(1),mols(1),expected_fragments,fragment_io,message)
      if (fragment_io /= 0) then
        message = 'post-trial GFN-FF fragment selections are invalid: '//trim(message)
        return
      end if
      if (.not.allocated(source_calc%calcs(1)%ff_dat%user_fraglist) .or. &
      & size(source_calc%calcs(1)%ff_dat%user_fraglist) /= nat .or. &
      & any(source_calc%calcs(1)%ff_dat%user_fraglist /= expected_fragments)) then
        message = 'post-trial GFN-FF fragment map differs from parsed fragment selections'
        return
      end if
      deallocate(expected_fragments)
    else if (allocated(source_calc%calcs(1)%ff_dat%user_fraglist)) then
      message = 'post-trial GFN-FF state has a fragment map without configured fragment selections'
      return
    end if
    if (allocated(source_calc%freezelist)) then
      if (size(source_calc%freezelist) /= nat .or. count(source_calc%freezelist) /= source_calc%nfreeze) then
        message = 'post-trial CREST freeze mask/count is inconsistent'
        return
      end if
      if (.not.allocated(source_calc%calcs(1)%ff_dat%frozen_mask) .or. &
      & size(source_calc%calcs(1)%ff_dat%frozen_mask) /= nat .or. &
      & .not.all(source_calc%calcs(1)%ff_dat%frozen_mask .eqv. source_calc%freezelist)) then
        message = 'post-trial GFN-FF frozen mask differs from parsed CREST freeze state'
        return
      end if
    else if (source_calc%nfreeze /= 0) then
      message = 'post-trial GFN-FF freeze count is nonzero without a CREST freeze mask'
      return
    else if (allocated(source_calc%calcs(1)%ff_dat%frozen_mask)) then
      if (size(source_calc%calcs(1)%ff_dat%frozen_mask) /= nat .or. &
      & any(source_calc%calcs(1)%ff_dat%frozen_mask)) then
        message = 'post-trial GFN-FF frozen mask is inconsistent with no requested freeze'
        return
      end if
    end if

    if (production_threads < 1) then
      message = 'calculator oracle received a nonpositive production thread count'
      return
    end if
    allocate(oracle_indices(size(mols)))
    call build_oracle_probe_indices(mols,oracle_indices,noracle)
    if (noracle < 1) then
      message = 'calculator oracle could not identify any prepared input geometry'
      return
    end if

    ! These deallocations are on the independent deep copy only.  source_calc
    ! and env%calc remain untouched and available to the inherited side of
    ! every oracle pair.  The clean side now models the actual capsule decoder:
    ! it reconstructs the live one-level weight from level%weight and allocates
    ! calculator scratch on its first energy/gradient call.
    deallocate(capsule_calc%calcs(1)%ff_dat)
    if (allocated(capsule_calc%etmp)) deallocate(capsule_calc%etmp)
    if (allocated(capsule_calc%grdtmp)) deallocate(capsule_calc%grdtmp)
    if (allocated(capsule_calc%eweight)) deallocate(capsule_calc%eweight)
    if (allocated(capsule_calc%weightbackup)) deallocate(capsule_calc%weightbackup)
    if (allocated(capsule_calc%activebackup)) deallocate(capsule_calc%activebackup)
    if (allocated(capsule_calc%etmp2)) deallocate(capsule_calc%etmp2)
    if (allocated(capsule_calc%grdtmp2)) deallocate(capsule_calc%grdtmp2)
    if (allocated(capsule_calc%eweight2)) deallocate(capsule_calc%eweight2)
    if (allocated(capsule_calc%grdfix)) deallocate(capsule_calc%grdfix)
    saved_dynamic = omp_get_dynamic()
    saved_threads = omp_get_max_threads()
    saved_levels = omp_get_max_active_levels()
#ifdef WITH_OPENBLAS
    saved_blas_threads = openblas_get_num_threads()
#endif
    engrad_before = engrad_total
    call omp_set_dynamic(.false.)
    call omp_set_max_active_levels(1)
    call configure_oracle_threading(1,io,message)
    if (io /= 0) goto 900

    ! Contract 1: isolate calculator history from floating reduction order.
    ! Both sides execute with one actual OpenMP thread and one BLAS thread, so
    ! inherited-versus-fresh equality is deliberately bitwise and fail-closed.
    do probe_id = 1,noracle
      index = oracle_indices(probe_id)

      ! Independent clones: no state flows from one input structure to another.
      inherited_calc = source_calc
      fresh_calc = capsule_calc
      probe_mol = mols(index)
      ! Native MTD workers reuse calculators across independent trajectories.
      ! Model the same one-shot geometry-boundary rebuild that the native task
      ! path requests before dynamics; no topology, EEQ, D3, or workspace state
      ! is released.
      call request_gfnff_hbond_update(inherited_calc)
      call prepare_gfnff_topology(probe_mol,fresh_calc,io,ninitialized,ngfnff)
      if (io /= 0 .or. ninitialized /= 1 .or. ngfnff /= 1) then
        write(message,'(a,i0,a,i0,a,i0,a,i0)') 'fresh topology setup failed for input ', &
        & probe_id,': iostat=',io,', initialized=',ninitialized,', GFN-FF levels=',ngfnff
        goto 900
      end if
      allocate(inherited_gradient(3,nat),fresh_gradient(3,nat))
      call engrad(probe_mol,inherited_calc,inherited_energy,inherited_gradient,io)
      if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'inherited calculator probe failed for input ',probe_id, &
        & ' with iostat=',io
        goto 900
      end if
      call engrad(probe_mol,fresh_calc,fresh_energy,fresh_gradient,io)
      if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'fresh calculator probe failed for input ',probe_id, &
        & ' with iostat=',io
        goto 900
      end if
      call compare_gfnff_probe(inherited_calc,fresh_calc,inherited_energy,fresh_energy, &
      & inherited_gradient,fresh_gradient,nat,.true.,same,differences,message)
      if (.not.same) then
        message = 'input '//trim(integer_string(probe_id))//' first-frame '//trim(message)
        goto 900
      end if

      call select_oracle_probe_atom(mddats(index),source_calc,nat,probe_atom,io,message)
      if (io /= 0) goto 900
      if (probe_atom == 0) then
        write(stdout,'(1x,a,i0,a)') 'Process MTD history oracle input ',probe_id, &
        & ': controlled-motion check skipped because every atom is frozen.'
      else
        ! 2^-7 Bohr gives an exact binary displacement.  The selected atom is
        ! resolved from the prepared MTD/freeze state, never from atom ordering.
        write(stdout,'(1x,a,i0,a,i0)') 'Process MTD history oracle input ',probe_id, &
        & ': controlled-motion atom=',probe_atom
        probe_mol%xyz(1,probe_atom) = probe_mol%xyz(1,probe_atom)+scale(1.0_wp,-7)
        call engrad(probe_mol,inherited_calc,inherited_energy,inherited_gradient,io)
        if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'inherited second probe failed for input ',probe_id, &
        & ' with iostat=',io
        goto 900
        end if
        call engrad(probe_mol,fresh_calc,fresh_energy,fresh_gradient,io)
        if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'fresh second probe failed for input ',probe_id, &
        & ' with iostat=',io
        goto 900
        end if
        call compare_gfnff_probe(inherited_calc,fresh_calc,inherited_energy,fresh_energy, &
        & inherited_gradient,fresh_gradient,nat,.true.,same,differences,message)
        if (.not.same) then
        message = 'input '//trim(integer_string(probe_id))//' controlled-history '//trim(message)
        goto 900
        end if
      end if
      deallocate(inherited_gradient,fresh_gradient)
    end do
    write(stdout,'(1x,a,i0,a,i0,a)') 'Process MTD T1 history oracle: ',noracle,'/',noracle, &
    & ' inherited-vs-fresh unique-geometry pairs bitwise exact.'

    ! Contract 2 is a development-time validation of harmless OpenMP reduction
    ! order differences, not a production integrity check.  Keep it available
    ! as an exact opt-in without allowing sub-picohartree variation to block an
    ! otherwise valid process-isolated MTD batch.
    run_thread_oracle = .false.
    thread_oracle_value = ''
    call get_environment_variable(thread_oracle_optin_name,thread_oracle_value, &
    & length=thread_oracle_length,status=thread_oracle_status,trim_name=.true.)
    if (thread_oracle_status == 0) then
      if (thread_oracle_length /= 1 .or. thread_oracle_value(1:1) /= '1') then
        message = thread_oracle_optin_name//' accepts only the exact value 1 when present'
        goto 900
      end if
      run_thread_oracle = .true.
    else if (thread_oracle_status /= 1) then
      message = thread_oracle_optin_name// &
      & ' is present but unreadable or longer than 32 bytes'
      goto 900
    end if
    if (.not.run_thread_oracle) then
      status = status_normal
      message = ''
      write(stdout,'(1x,a)') 'Process MTD production-thread numerical oracle: skipped; set '// &
      & thread_oracle_optin_name//'=1 for strict validation.'
      goto 900
    end if

    ! Contract 2: quantify only the reduction-order effect of the production
    ! production-thread path.  Both calculators start fresh from the same canonical
    ! capsule settings and topology.  Numerical outputs use the registered
    ! absolute envelope; all discrete/static/mutable state remains exact.
    max_production_differences = 0.0_wp
    do probe_id = 1,noracle
      index = oracle_indices(probe_id)

      fresh_t1_calc = capsule_calc
      fresh_tk_calc = capsule_calc
      probe_mol = mols(index)
      call configure_oracle_threading(1,io,message)
      if (io /= 0) goto 900
      call prepare_gfnff_topology(probe_mol,fresh_t1_calc,io,ninitialized,ngfnff)
      if (io /= 0 .or. ninitialized /= 1 .or. ngfnff /= 1) then
        write(message,'(a,i0,a,i0,a,i0,a,i0)') 'fresh T1 topology setup failed for input ', &
        & probe_id,': iostat=',io,', initialized=',ninitialized,', GFN-FF levels=',ngfnff
        goto 900
      end if
      call configure_oracle_threading(production_threads,io,message)
      if (io /= 0) goto 900
      call prepare_gfnff_topology(probe_mol,fresh_tk_calc,io,ninitialized,ngfnff)
      if (io /= 0 .or. ninitialized /= 1 .or. ngfnff /= 1) then
        write(message,'(a,i0,a,i0,a,i0,a,i0)') 'fresh TK topology setup failed for input ', &
        & probe_id,': iostat=',io,', initialized=',ninitialized,', GFN-FF levels=',ngfnff
        goto 900
      end if
      allocate(fresh_t1_gradient(3,nat),fresh_tk_gradient(3,nat))

      call configure_oracle_threading(1,io,message)
      if (io /= 0) goto 900
      call engrad(probe_mol,fresh_t1_calc,fresh_t1_energy,fresh_t1_gradient,io)
      if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'fresh T1 first-frame probe failed for input ', &
        & probe_id,' with iostat=',io
        goto 900
      end if
      call configure_oracle_threading(production_threads,io,message)
      if (io /= 0) goto 900
      call engrad(probe_mol,fresh_tk_calc,fresh_tk_energy,fresh_tk_gradient,io)
      if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'fresh TK first-frame probe failed for input ', &
        & probe_id,' with iostat=',io
        goto 900
      end if
      call compare_gfnff_probe(fresh_t1_calc,fresh_tk_calc,fresh_t1_energy,fresh_tk_energy, &
      & fresh_t1_gradient,fresh_tk_gradient,nat,.false.,same,differences,message)
      max_production_differences = max(max_production_differences,differences)
      if (.not.same) then
        message = 'input '//trim(integer_string(probe_id))// &
        & ' fresh-T1/fresh-TK first-frame '//trim(message)
        goto 900
      end if

      call select_oracle_probe_atom(mddats(index),source_calc,nat,probe_atom,io,message)
      if (io /= 0) goto 900
      if (probe_atom == 0) then
        write(stdout,'(1x,a,i0,a)') 'Process MTD strict thread oracle input ',probe_id, &
        & ': controlled-motion check skipped because every atom is frozen.'
      else
        write(stdout,'(1x,a,i0,a,i0)') 'Process MTD strict thread oracle input ',probe_id, &
        & ': controlled-motion atom=',probe_atom
        probe_mol%xyz(1,probe_atom) = probe_mol%xyz(1,probe_atom)+scale(1.0_wp,-7)
        call configure_oracle_threading(1,io,message)
      if (io /= 0) goto 900
        call engrad(probe_mol,fresh_t1_calc,fresh_t1_energy,fresh_t1_gradient,io)
        if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'fresh T1 controlled-history probe failed for input ', &
        & probe_id,' with iostat=',io
        goto 900
        end if
        call configure_oracle_threading(production_threads,io,message)
      if (io /= 0) goto 900
        call engrad(probe_mol,fresh_tk_calc,fresh_tk_energy,fresh_tk_gradient,io)
        if (io /= 0) then
        write(message,'(a,i0,a,i0)') 'fresh TK controlled-history probe failed for input ', &
        & probe_id,' with iostat=',io
        goto 900
        end if
        call compare_gfnff_probe(fresh_t1_calc,fresh_tk_calc,fresh_t1_energy,fresh_tk_energy, &
        & fresh_t1_gradient,fresh_tk_gradient,nat,.false.,same,differences,message)
        max_production_differences = max(max_production_differences,differences)
        if (.not.same) then
        message = 'input '//trim(integer_string(probe_id))// &
        & ' fresh-T1/fresh-TK controlled-history '//trim(message)
        goto 900
        end if
      end if
      deallocate(fresh_t1_gradient,fresh_tk_gradient)
    end do

    status = status_normal
    message = ''
    write(stdout,'(1x,a,i0,a,i0,a,es10.3,a,4(1x,es12.4))') &
    & 'Process MTD production-thread numerical oracle: ',noracle,'/',noracle, &
    & ' fresh T1/TK pairs within absolute ',production_thread_oracle_atol, &
    & '; max |dE|/component/q/gradient=', &
    & max_production_differences

900 continue
    if (allocated(inherited_gradient)) deallocate(inherited_gradient)
    if (allocated(fresh_gradient)) deallocate(fresh_gradient)
    if (allocated(fresh_t1_gradient)) deallocate(fresh_t1_gradient)
    if (allocated(fresh_tk_gradient)) deallocate(fresh_tk_gradient)
    if (allocated(oracle_indices)) deallocate(oracle_indices)
    if (allocated(expected_fragments)) deallocate(expected_fragments)
    engrad_total = engrad_before
    call omp_set_num_threads(saved_threads)
    call omp_set_max_active_levels(saved_levels)
    call omp_set_dynamic(saved_dynamic)
#ifdef WITH_OPENBLAS
    call openblasset(saved_blas_threads)
#endif
#endif
  end subroutine prepare_canonical_capsule_calculator

  subroutine build_oracle_probe_indices(mols,indices,nunique)
    type(coord),intent(in) :: mols(:)
    integer,intent(out) :: indices(:)
    integer,intent(out) :: nunique
    integer :: i,j
    logical :: duplicate

    nunique = 0
    indices = 0
    do i = 1,size(mols)
      duplicate = .false.
      do j = 1,nunique
        if (mols(i)%nat == mols(indices(j))%nat .and. &
        & all(mols(i)%at == mols(indices(j))%at) .and. &
        & same_real_bits_2d(mols(i)%xyz,mols(indices(j))%xyz)) then
          duplicate = .true.
          exit
        end if
      end do
      if (.not.duplicate) then
        nunique = nunique+1
        indices(nunique) = i
      end if
    end do
  end subroutine build_oracle_probe_indices

  subroutine configure_oracle_threading(expected_threads,io,message)
    integer,intent(in) :: expected_threads
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    integer :: observed_threads

    io = 0
    message = ''
    observed_threads = 0
    call omp_set_num_threads(expected_threads)
    call openblasset(1)
    !$omp parallel default(none) shared(observed_threads)
    !$omp single
    observed_threads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    if (observed_threads /= expected_threads) then
      io = 1
      write(message,'(a,i0,a,i0)') 'calculator oracle requested ',expected_threads, &
      & ' OpenMP threads but formed ',observed_threads
      return
    end if
#ifdef WITH_OPENBLAS
    if (openblas_get_num_threads() /= 1) then
      io = 2
      write(message,'(a,i0)') 'calculator oracle requires one BLAS thread but observed ', &
      & openblas_get_num_threads()
    end if
#endif
  end subroutine configure_oracle_threading

  subroutine select_oracle_probe_atom(md,calc,nat,atom,io,message)
    type(mddata),intent(in) :: md
    type(calcdata),intent(in) :: calc
    integer,intent(in) :: nat
    integer,intent(out) :: atom,io
    character(len=*),intent(out) :: message
    integer :: i

    atom = 0
    io = 0
    message = ''
    if (nat < 1) then
      io = 1
      message = 'oracle probe selection received a nonpositive atom count'
      return
    end if
    ! A CV can include frozen atoms, but the controlled-motion oracle must
    ! never perturb one.  Validate this state before choosing an included atom.
    if (allocated(calc%freezelist)) then
      if (size(calc%freezelist) /= nat) then
        io = 4
        message = 'oracle freeze mask length differs from molecule atom count'
        return
      end if
    else if (calc%nfreeze /= 0) then
      io = 5
      message = 'oracle freeze count is nonzero without a freeze mask'
      return
    end if
    if (md%npot > 0 .and. allocated(md%mtd)) then
      if (size(md%mtd) < 1) then
        io = 2
        message = 'oracle probe selection found an empty MTD potential array'
        return
      end if
      if (allocated(md%mtd(1)%atinclude)) then
        if (size(md%mtd(1)%atinclude) /= nat) then
          io = 3
          message = 'oracle MTD inclusion mask length differs from molecule atom count'
          return
        end if
        do i = 1,nat
          if (md%mtd(1)%atinclude(i) .and. &
          & (.not.allocated(calc%freezelist) .or. .not.calc%freezelist(i))) then
            atom = i
            return
          end if
        end do
      end if
    end if
    if (allocated(calc%freezelist)) then
      do i = 1,nat
        if (.not.calc%freezelist(i)) then
          atom = i
          return
        end if
      end do
    else
      atom = 1
    end if
  end subroutine select_oracle_probe_atom

  subroutine compare_gfnff_probe(calc_a,calc_b,energy_a,energy_b,gradient_a,gradient_b, &
  & expected_atoms,require_bitwise,same,differences,message)
    type(calcdata),intent(in) :: calc_a,calc_b
    real(wp),intent(in) :: energy_a,energy_b
    real(wp),intent(in) :: gradient_a(:,:),gradient_b(:,:)
    integer,intent(in) :: expected_atoms
    logical,intent(in) :: require_bitwise
    logical,intent(out) :: same
    real(wp),intent(out) :: differences(4)
    character(len=*),intent(out) :: message
#ifdef WITH_GFNFF
    real(wp) :: components_a(20),components_b(20)
    integer :: i
#endif
    same = .false.
    differences = 0.0_wp
    message = ''
#ifndef WITH_GFNFF
    message = 'GFN-FF probe comparison is unavailable in this build'
#else
    if (.not.allocated(calc_a%calcs) .or. .not.allocated(calc_b%calcs)) then
      message = 'calculator levels are unallocated after probe'
      return
    end if
    if (size(calc_a%calcs) /= 1 .or. size(calc_b%calcs) /= 1) then
      message = 'calculator level count changed during probe'
      return
    end if
    if (.not.allocated(calc_a%calcs(1)%ff_dat) .or. &
    & .not.allocated(calc_b%calcs(1)%ff_dat)) then
      message = 'calculator state is unallocated after probe'
      return
    end if
    if (.not.allocated(calc_a%calcs(1)%ff_dat%res) .or. &
    & .not.allocated(calc_b%calcs(1)%ff_dat%res)) then
      message = 'GFN-FF result components are unallocated after probe'
      return
    end if
    if (.not.allocated(calc_a%calcs(1)%ff_dat%nlist) .or. &
    & .not.allocated(calc_b%calcs(1)%ff_dat%nlist)) then
      message = 'GFN-FF neighbour list is unallocated after probe'
      return
    end if
    if (.not.allocated(calc_a%calcs(1)%ff_dat%nlist%q) .or. &
    & .not.allocated(calc_b%calcs(1)%ff_dat%nlist%q)) then
      message = 'GFN-FF charges are unallocated after probe'
      return
    end if
    if (expected_atoms < 1 .or. size(calc_a%calcs(1)%ff_dat%nlist%q) /= expected_atoms .or. &
    & size(calc_b%calcs(1)%ff_dat%nlist%q) /= expected_atoms) then
      message = 'GFN-FF charge vector has the wrong atom count'
      return
    end if
    if (size(gradient_a,1) /= 3 .or. size(gradient_b,1) /= 3 .or. &
    & size(gradient_a,2) /= expected_atoms .or. size(gradient_b,2) /= expected_atoms) then
      message = 'probe gradient has the wrong dimensions'
      return
    end if
    if (.not.same_gfnff_static_state(calc_a,calc_b)) then
      message = 'exact topology/fragment/freeze state differs'
      return
    end if
    components_a = result_components(calc_a)
    components_b = result_components(calc_b)
    if (.not.ieee_is_finite(energy_a) .or. .not.ieee_is_finite(energy_b) .or. &
    & .not.all(ieee_is_finite(components_a)) .or. &
    & .not.all(ieee_is_finite(components_b)) .or. &
    & .not.all(ieee_is_finite(calc_a%calcs(1)%ff_dat%nlist%q)) .or. &
    & .not.all(ieee_is_finite(calc_b%calcs(1)%ff_dat%nlist%q)) .or. &
    & .not.all(ieee_is_finite(gradient_a)) .or. &
    & .not.all(ieee_is_finite(gradient_b))) then
      message = 'probe produced nonfinite scientific output'
      return
    end if
    if (allocated(calc_a%freezelist)) then
      if (size(calc_a%freezelist) /= expected_atoms) then
        message = 'probe calculator freeze mask has the wrong atom count'
        return
      end if
      do i = 1,expected_atoms
        if (calc_a%freezelist(i)) then
          if (any(transfer(gradient_a(:,i),0_int64,3) /= 0_int64) .or. &
          & any(transfer(gradient_b(:,i),0_int64,3) /= 0_int64)) then
            message = 'probe did not return exact positive-zero gradients for every frozen atom'
            return
          end if
        end if
      end do
    end if
    differences = [abs(energy_a-energy_b), &
    & maxval(abs(components_a-components_b)), &
    & maxval(abs(calc_a%calcs(1)%ff_dat%nlist%q- &
    &            calc_b%calcs(1)%ff_dat%nlist%q)), &
    & maxval(abs(gradient_a-gradient_b))]
    if (require_bitwise) then
      if (.not.same_real_bits_scalar(energy_a,energy_b) .or. &
      & .not.same_real_bits_1d(components_a,components_b) .or. &
      & .not.same_real_bits_1d(calc_a%calcs(1)%ff_dat%nlist%q, &
      &                        calc_b%calcs(1)%ff_dat%nlist%q) .or. &
      & .not.same_real_bits_2d(gradient_a,gradient_b)) then
        write(message,'(a,4(1x,es12.4))') &
        & 'T1 bitwise scientific outputs differ; |dE|, max component/q/gradient=',differences
        return
      end if
    else
      if (any(differences > production_thread_oracle_atol)) then
        write(message,'(a,es10.3,a,4(1x,es12.4))') &
        & 'production-thread reduction-order difference exceeds absolute ',production_thread_oracle_atol, &
        & '; |dE|, max component/q/gradient=',differences
        return
      end if
    end if
    if (.not.same_gfnff_hbxb_state(calc_a,calc_b)) then
      message = 'exact mutable HB/XB neighbour-list state differs'
      return
    end if
    same = .true.
#endif
  end subroutine compare_gfnff_probe

  logical function same_gfnff_static_state(calc_a,calc_b) result(same)
    type(calcdata),intent(in) :: calc_a,calc_b
    integer :: i
#ifdef WITH_GFNFF
    same = allocated(calc_a%calcs) .and. allocated(calc_b%calcs)
    if (.not.same) return
    same = size(calc_a%calcs) == 1 .and. size(calc_b%calcs) == 1
    if (.not.same) return
    same = calc_a%nfreeze == calc_b%nfreeze .and. &
    & same_logical_alloc_1d(calc_a%freezelist,calc_b%freezelist)
    if (.not.same) return
    same = allocated(calc_a%calcs(1)%gff_fragments) .eqv. &
    & allocated(calc_b%calcs(1)%gff_fragments)
    if (.not.same) return
    if (allocated(calc_a%calcs(1)%gff_fragments)) then
      same = size(calc_a%calcs(1)%gff_fragments) == size(calc_b%calcs(1)%gff_fragments)
      if (.not.same) return
      do i = 1,size(calc_a%calcs(1)%gff_fragments)
        same = allocated(calc_a%calcs(1)%gff_fragments(i)%value) .eqv. &
        & allocated(calc_b%calcs(1)%gff_fragments(i)%value)
        if (.not.same) return
        if (allocated(calc_a%calcs(1)%gff_fragments(i)%value)) then
          same = len(calc_a%calcs(1)%gff_fragments(i)%value) == &
          & len(calc_b%calcs(1)%gff_fragments(i)%value) .and. &
          & calc_a%calcs(1)%gff_fragments(i)%value == &
          & calc_b%calcs(1)%gff_fragments(i)%value
          if (.not.same) return
        end if
      end do
    end if
    same = allocated(calc_a%calcs(1)%ff_dat) .and. &
    & allocated(calc_b%calcs(1)%ff_dat)
    if (.not.same) return
    associate(a => calc_a%calcs(1)%ff_dat,b => calc_b%calcs(1)%ff_dat)
      same = a%ichrg == b%ichrg .and. &
      & same_real_bits_scalar(a%accuracy,b%accuracy) .and. &
      & (a%make_chrg .eqv. b%make_chrg) .and. a%version == b%version .and. &
      & (a%update .eqv. b%update) .and. (a%write_topo .eqv. b%write_topo) .and. &
      & same_int_alloc_1d(a%user_fraglist,b%user_fraglist) .and. &
      & same_logical_alloc_1d(a%frozen_mask,b%frozen_mask)
    end associate
    if (.not.same) return
    same = same_gfnff_topology(calc_a,calc_b)
#else
    same = .false.
#endif
  end function same_gfnff_static_state

  logical function same_gfnff_topology(calc_a,calc_b) result(same)
    type(calcdata),intent(in) :: calc_a,calc_b
#ifdef WITH_GFNFF
    same = allocated(calc_a%calcs(1)%ff_dat%topo) .and. &
    & allocated(calc_b%calcs(1)%ff_dat%topo)
    if (.not.same) return
    associate(a => calc_a%calcs(1)%ff_dat%topo,b => calc_b%calcs(1)%ff_dat%topo)
      same = a%nbond == b%nbond .and. a%nangl == b%nangl .and. &
      & a%ntors == b%ntors .and. a%nathbH == b%nathbH .and. &
      & a%nathbAB == b%nathbAB .and. a%natxbAB == b%natxbAB .and. &
      & a%nbatm == b%nbatm .and. a%nfrag == b%nfrag .and. &
      & a%maxsystem == b%maxsystem .and. a%bond_hb_nr == b%bond_hb_nr .and. &
      & a%b_max == b%b_max .and. a%nbond_blist == b%nbond_blist .and. &
      & a%nbond_vbond == b%nbond_vbond .and. a%nangl_alloc == b%nangl_alloc .and. &
      & a%ntors_alloc == b%ntors_alloc .and. a%read_file_type == b%read_file_type .and. &
      & a%nsystem == b%nsystem
      if (.not.same) return
      ! topo%filename intentionally names different byte-identical original and
      ! private audit copies.  It is provenance, not topology content; the two
      ! files were compared byte-for-byte before this oracle.  Compare every
      ! scientific scalar/array held by the loaded topology instead.
      same = same_character_alloc(a%refcharges,b%refcharges) .and. &
      & same_int_alloc_2d(a%nb,b%nb) .and. same_int_alloc_1d(a%hyb,b%hyb) .and. &
      & same_int_alloc_1d(a%bpair,b%bpair) .and. &
      & same_int_alloc_2d(a%blist,b%blist) .and. &
      & same_int_alloc_2d(a%alist,b%alist) .and. &
      & same_int_alloc_2d(a%tlist,b%tlist) .and. &
      & same_int_alloc_2d(a%b3list,b%b3list) .and. &
      & same_int_alloc_2d(a%sTorsl,b%sTorsl) .and. &
      & same_real_bits_alloc_1d(a%pbo,b%pbo)
      if (.not.same) return
      same = same_int_alloc_1d(a%nr_hb,b%nr_hb) .and. &
      & same_int_alloc_2d(a%bond_hb_AH,b%bond_hb_AH) .and. &
      & same_int_alloc_2d(a%bond_hb_B,b%bond_hb_B) .and. &
      & same_int_alloc_1d(a%bond_hb_Bn,b%bond_hb_Bn) .and. &
      & same_int_alloc_2d(a%hbatABl,b%hbatABl) .and. &
      & same_int_alloc_2d(a%xbatABl,b%xbatABl) .and. &
      & same_int_alloc_1d(a%hbatHl,b%hbatHl) .and. &
      & same_int_alloc_1d(a%fraglist,b%fraglist) .and. &
      & same_int_alloc_1d(a%qpdb,b%qpdb)
      if (.not.same) return
      same = same_real_bits_alloc_2d(a%vbond,b%vbond) .and. &
      & same_real_bits_alloc_2d(a%vangl,b%vangl) .and. &
      & same_real_bits_alloc_2d(a%vtors,b%vtors) .and. &
      & same_real_bits_alloc_1d(a%chieeq,b%chieeq) .and. &
      & same_real_bits_alloc_1d(a%gameeq,b%gameeq) .and. &
      & same_real_bits_alloc_1d(a%alpeeq,b%alpeeq) .and. &
      & same_real_bits_alloc_1d(a%alphanb,b%alphanb) .and. &
      & same_real_bits_alloc_1d(a%qa,b%qa) .and. &
      & same_real_bits_alloc_2d(a%xyze0,b%xyze0) .and. &
      & same_real_bits_alloc_1d(a%zetac6,b%zetac6) .and. &
      & same_real_bits_alloc_1d(a%qfrag,b%qfrag) .and. &
      & same_real_bits_alloc_1d(a%hbbas,b%hbbas) .and. &
      & same_real_bits_alloc_1d(a%hbaci,b%hbaci) .and. &
      & same_int_alloc_2d(a%ispinsyst,b%ispinsyst) .and. &
      & same_int_alloc_1d(a%nspinsyst,b%nspinsyst)
      if (.not.same) return
      same = same_real_bits_scalar(a%dispm%s6,b%dispm%s6) .and. &
      & same_real_bits_scalar(a%dispm%s8,b%dispm%s8) .and. &
      & same_real_bits_scalar(a%dispm%s10,b%dispm%s10) .and. &
      & same_real_bits_scalar(a%dispm%a1,b%dispm%a1) .and. &
      & same_real_bits_scalar(a%dispm%a2,b%dispm%a2) .and. &
      & same_real_bits_scalar(a%dispm%s9,b%dispm%s9) .and. &
      & a%dispm%alp == b%dispm%alp .and. &
      & same_real_bits_scalar(a%dispm%wf,b%dispm%wf) .and. &
      & same_real_bits_scalar(a%dispm%g_a,b%dispm%g_a) .and. &
      & same_real_bits_scalar(a%dispm%g_c,b%dispm%g_c) .and. &
      & same_int_alloc_1d(a%dispm%atoms,b%dispm%atoms) .and. &
      & same_int_alloc_1d(a%dispm%nref,b%dispm%nref) .and. &
      & same_int_alloc_2d(a%dispm%ncount,b%dispm%ncount) .and. &
      & same_real_bits_alloc_2d(a%dispm%cn,b%dispm%cn) .and. &
      & same_real_bits_alloc_2d(a%dispm%q,b%dispm%q) .and. &
      & same_real_bits_alloc_3d(a%dispm%alpha,b%dispm%alpha) .and. &
      & same_real_bits_alloc_4d(a%dispm%c6,b%dispm%c6)
    end associate
#else
    same = .false.
#endif
  end function same_gfnff_topology

  function result_components(calc) result(values)
    type(calcdata),intent(in) :: calc
    real(wp) :: values(20)
#ifdef WITH_GFNFF
    associate(res => calc%calcs(1)%ff_dat%res)
      values = [res%e_total,res%e_rep,res%e_es,res%e_disp,res%e_xb, &
      & res%g_born,res%g_sasa,res%g_hb,res%g_shift,res%dipole, &
      & res%g_solv,res%gnorm,res%e_bond,res%e_angl,res%e_tors, &
      & res%e_hb,res%e_batm,res%e_ext]
    end associate
#else
    values = 0.0_wp
#endif
  end function result_components

  logical function same_gfnff_hbxb_state(calc_a,calc_b) result(same)
    type(calcdata),intent(in) :: calc_a,calc_b
#ifdef WITH_GFNFF
    associate(a => calc_a%calcs(1)%ff_dat%nlist, &
    &         b => calc_b%calcs(1)%ff_dat%nlist)
      same = (a%initialized .eqv. b%initialized) .and. &
      & (a%force_hbond_update .eqv. b%force_hbond_update) .and. &
      & a%nhb1 == b%nhb1 .and. a%nhb2 == b%nhb2 .and. a%nxb == b%nxb
      if (.not.same) return
      ! Charges are scientific numerical output: bitwise at T1 and bounded at
      ! the production thread count by compare_gfnff_probe.  Keep the actual HB/XB history exact here.
      same = same_real_bits_alloc_2d(a%hbrefgeo,b%hbrefgeo) .and. &
      & same_active_hb_list(a%hblist1,b%hblist1,a%nhb1) .and. &
      & same_active_hb_list(a%hblist2,b%hblist2,a%nhb2) .and. &
      & same_active_hb_list(a%hblist3,b%hblist3,a%nxb)
    end associate
#else
    same = .false.
#endif
  end function same_gfnff_hbxb_state

  logical function same_real_bits_scalar(a,b) result(same)
    real(wp),intent(in) :: a,b
    same = transfer(a,0_int64) == transfer(b,0_int64)
  end function same_real_bits_scalar

  logical function same_character_alloc(a,b) result(same)
    character(len=:),allocatable,intent(in) :: a,b
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = len(a) == len(b) .and. a == b
  end function same_character_alloc

  logical function same_int_alloc_1d(a,b) result(same)
    integer,allocatable,intent(in) :: a(:),b(:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = size(a) == size(b)
    if (same) same = all(a == b)
  end function same_int_alloc_1d

  logical function same_int_alloc_2d(a,b) result(same)
    integer,allocatable,intent(in) :: a(:,:),b(:,:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = all(shape(a) == shape(b))
    if (same) same = all(a == b)
  end function same_int_alloc_2d

  logical function same_logical_alloc_1d(a,b) result(same)
    logical,allocatable,intent(in) :: a(:),b(:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = size(a) == size(b)
    if (same) same = all(a .eqv. b)
  end function same_logical_alloc_1d

  logical function same_real_bits_1d(a,b) result(same)
    real(wp),intent(in) :: a(:),b(:)
    same = size(a) == size(b)
    if (.not.same) return
    same = all(transfer(a,0_int64,size(a)) == transfer(b,0_int64,size(b)))
  end function same_real_bits_1d

  logical function same_real_bits_2d(a,b) result(same)
    real(wp),intent(in) :: a(:,:),b(:,:)
    same = all(shape(a) == shape(b))
    if (.not.same) return
    same = all(transfer(a,0_int64,size(a)) == transfer(b,0_int64,size(b)))
  end function same_real_bits_2d

  logical function same_real_bits_alloc_1d(a,b) result(same)
    real(wp),allocatable,intent(in) :: a(:),b(:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = same_real_bits_1d(a,b)
  end function same_real_bits_alloc_1d

  logical function same_real_bits_alloc_2d(a,b) result(same)
    real(wp),allocatable,intent(in) :: a(:,:),b(:,:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = same_real_bits_2d(a,b)
  end function same_real_bits_alloc_2d

  logical function same_real_bits_alloc_3d(a,b) result(same)
    real(wp),allocatable,intent(in) :: a(:,:,:),b(:,:,:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = all(shape(a) == shape(b))
    if (same) same = all(transfer(a,0_int64,size(a)) == transfer(b,0_int64,size(b)))
  end function same_real_bits_alloc_3d

  logical function same_real_bits_alloc_4d(a,b) result(same)
    real(wp),allocatable,intent(in) :: a(:,:,:,:),b(:,:,:,:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = all(shape(a) == shape(b))
    if (same) same = all(transfer(a,0_int64,size(a)) == transfer(b,0_int64,size(b)))
  end function same_real_bits_alloc_4d

  logical function same_active_hb_list(a,b,nactive) result(same)
    integer,allocatable,intent(in) :: a(:,:),b(:,:)
    integer,intent(in) :: nactive
    ! GFN-FF reads these arrays only through 1:nhb1, 1:nhb2, or
    ! 1:nxb.  The allocated second dimension is spare capacity (new() uses
    ! 5x/5x/3x over-allocation), and a rebuilt list can leave dead tail data.
    ! Require a valid capacity on both sides and compare every live entry,
    ! without treating allocator history as scientific state.
    same = allocated(a) .and. allocated(b)
    if (.not.same) return
    same = nactive >= 0 .and. size(a,1) == 3 .and. size(b,1) == 3
    if (.not.same) return
    same = size(a,2) >= nactive .and. size(b,2) >= nactive
    if (.not.same .or. nactive == 0) return
    same = all(a(:,1:nactive) == b(:,1:nactive))
  end function same_active_hb_list

  subroutine mtd_process_worker_dispatch(handled,status)
    logical,intent(out) :: handled
    integer,intent(out) :: status
    type(mtd_capsule_meta) :: meta
    type(coord) :: mol
    type(mddata) :: md
    type(calcdata) :: calc
    character(len=4096) :: flag,capsule_arg,workdir_arg,index_arg,workers_arg,threads_arg, &
    & parent_arg,cpu_arg
    character(len=1024) :: message
    character(kind=c_char),allocatable :: c_workdir(:),c_cpu_list(:)
    integer :: argc,io,worker_index,worker_count,thread_count,parent_pid,term
    integer(int64) :: clock_start,clock_end,clock_rate,clock_max
    integer(c_int) :: cstat
    real(wp) :: elapsed
    logical :: exists

    handled = .false.
    status = status_normal
    argc = command_argument_count()
    if (argc < 1) return
    call get_command_argument(1,flag)
    if (trim(flag) /= worker_flag) return
    handled = .true.
    if (argc /= 8) then
      write(stdout,'(1x,a)') 'Invalid internal process MTD worker argument count.'
      status = status_config
      return
    end if
    call get_command_argument(2,capsule_arg)
    call get_command_argument(3,workdir_arg)
    call get_command_argument(4,index_arg)
    call get_command_argument(5,workers_arg)
    call get_command_argument(6,threads_arg)
    call get_command_argument(7,parent_arg)
    call get_command_argument(8,cpu_arg)
    read(index_arg,*,iostat=io) worker_index
    if (io == 0) read(workers_arg,*,iostat=io) worker_count
    if (io == 0) read(threads_arg,*,iostat=io) thread_count
    if (io == 0) read(parent_arg,*,iostat=io) parent_pid
    if (io /= 0) then
      status = status_config
      return
    end if
    ! Establish orphan handling, the private cwd, CPU mask, and NUMA policy
    ! before allocating/first-touching any prepared scientific capsule payload.
    call to_c_string(trim(workdir_arg),c_workdir)
    call to_c_string(trim(cpu_arg),c_cpu_list)
    cstat = c_prepare_worker(c_workdir,int(worker_index,c_int),int(worker_count,c_int), &
    & int(thread_count,c_int),int(parent_pid,c_int),c_cpu_list)
    if (cstat /= 0_c_int) then
      write(stdout,'(1x,a,i0)') 'Worker affinity/orphan preparation failed; errno=',int(cstat)
      status = status_config
      return
    end if
    ! This must remain the first application OpenMP parallel region in an
    ! exec'd worker.  It proves every assigned CPU before capsule allocation
    ! or scientific initialization can perform useful work.
    call omp_set_dynamic(.false.)
    call omp_set_max_active_levels(1)
    call omp_set_num_threads(thread_count)
    call verify_worker_openmp_binding(trim(cpu_arg),thread_count,io,message)
    if (io /= 0) then
      write(stdout,'(1x,a)') 'Worker OpenMP binding verification failed: '//trim(message)
      status = status_config
      return
    end if
    call openblasset(1)
    call read_mtd_capsule(trim(capsule_arg),meta,mol,md,calc,io,message)
    if (io /= 0) then
      write(stdout,'(1x,a)') 'Worker capsule read failed: '//trim(message)
      status = status_config
      return
    end if
    call verify_mtd_capsule_roundtrip(trim(capsule_arg),meta,mol,md,calc,io,message)
    if (io /= 0) then
      write(stdout,'(1x,a)') 'Worker capsule canonical-byte verification failed: '//trim(message)
      status = status_config
      return
    end if
    if (meta%worker_index /= worker_index .or. worker_index < 1 .or. &
    & worker_index > worker_count .or. meta%parent_pid /= parent_pid .or. &
    & meta%inner_threads /= thread_count .or. meta%expected_atoms < 1 .or. &
    & trim(meta%workdir) /= trim(workdir_arg) .or. mol%nat /= meta%expected_atoms .or. &
    & md%md_index /= worker_index .or. md%samerand) then
      write(stdout,'(1x,a)') 'Worker capsule/argv invariant mismatch.'
      status = status_config
      return
    end if
    inquire(file=meta%topology_file,exist=exists,iostat=io)
    if (io /= 0 .or. .not.exists) then
      write(stdout,'(1x,a)') 'Private worker topology is absent.'
      status = status_config
      return
    end if
    call system_clock(clock_start,clock_rate,clock_max)
    if (clock_rate <= 0_int64) then
      write(stdout,'(1x,a)') 'Worker system clock is unavailable.'
      status = status_config
      return
    end if
    call dynamics(mol,md,calc,.false.,term)
    call system_clock(clock_end)
    if (clock_end >= clock_start) then
      elapsed = real(clock_end-clock_start,wp)/real(clock_rate,wp)
    else
      elapsed = real((clock_max-clock_start)+clock_end+1_int64,wp)/real(clock_rate,wp)
    end if
    call write_worker_result(meta%result_file,worker_index,term,elapsed,engrad_total,io,message)
    if (io /= 0) then
      write(stdout,'(1x,a)') 'Worker result write failed: '//trim(message)
      status = status_failed
      return
    end if
    if (term /= 0) then
      status = status_failed
    else
      status = status_normal
    end if
  end subroutine mtd_process_worker_dispatch

  subroutine validate_process_configuration(env,mols,mddats,nsim,outer_threads, &
  & inner_threads,status,message)
    type(systemdata),intent(in) :: env
    integer,intent(in) :: nsim,outer_threads,inner_threads
    type(coord),intent(in) :: mols(nsim)
    type(mddata),intent(in) :: mddats(nsim)
    integer,intent(out) :: status
    character(len=*),intent(out) :: message
    integer :: i,nat,fragment_io
    integer,allocatable :: fragment_ids(:)
    integer(c_int) :: layout_status
    integer(int64) :: required_threads

    status = status_config
    message = ''
#ifndef WITH_OPENBLAS
    message = 'process MTD requires a compiled OpenBLAS build'
    return
#endif
    if (nsim < 2) then
      message = 'process MTD batch scheduler requires at least two trajectories'
      return
    end if
    if (outer_threads /= nsim) then
      write(message,'(a,i0,a,i0,a)') 'process MTD requires one concurrent worker per trajectory: ', &
      & nsim,' trajectories but outer scheduler resolved ',outer_threads, &
      & '; increase -TMD/-T or reduce the batch size'
      return
    end if
    if (inner_threads < 1) then
      message = 'process MTD resolved a nonpositive calculator thread count'
      return
    end if
    if (.not.allocated(mols(1)%at) .or. .not.allocated(mols(1)%xyz) .or. &
    & mols(1)%nat < 1) then
      message = 'process MTD first prepared molecule lacks a positive atom count, atoms, or coordinates'
      return
    end if
    nat = mols(1)%nat
    required_threads = int(nsim,int64)*int(inner_threads,int64)
    if (required_threads > int(env%Threads,int64)) then
      write(message,'(a,i0,a,i0,a,i0,a)') 'process MTD layout needs ',required_threads, &
      & ' CPUs (workers=',nsim,', threads/worker=',inner_threads, &
      & ') but the MTD thread budget is smaller'
      return
    end if
    layout_status = c_validate_layout(int(nsim,c_int),int(inner_threads,c_int))
    if (layout_status /= 0_c_int) then
      write(message,'(a,i0,a,i0,a,i0)') 'process MTD L3/NUMA layout is invalid: errno=', &
      & int(layout_status),', workers=',nsim,', threads/worker=',inner_threads
      return
    end if
    if (env%calc%ncalculations /= 1 .or. .not.allocated(env%calc%calcs)) then
      message = 'process MTD requires exactly one calculator level'
      return
    end if
    if (size(env%calc%calcs) /= 1) then
      message = 'process MTD calculator allocation/count is inconsistent'
      return
    end if
    if (env%calc%calcs(1)%id /= jobtype%gfnff .or. .not.env%calc%calcs(1)%active) then
      message = 'process MTD requires one active built-in GFN-FF calculator'
      return
    end if
#ifndef WITH_GFNFF
    message = 'process MTD requires a compiled GFN-FF implementation'
    return
#else
    if (.not.allocated(env%calc%calcs(1)%ff_dat)) then
      message = 'process MTD requires initialized parent GFN-FF runtime state'
      return
    end if
    if (.not.allocated(env%calc%calcs(1)%ff_dat%topo)) then
      message = 'process MTD requires initialized parent GFN-FF topology state'
      return
    end if
#endif
    if (env%calc%nscans /= 0 .or. allocated(env%calc%ONIOM)) then
      message = 'process MTD rejects scans and ONIOM because their prepared state is not serialized'
      return
    end if
    if (.not.supported_process_nci_walls(env%calc,nat)) then
      message = 'process MTD requires exactly one canonical automatic all-atom NCI log-Fermi wall'
      return
    end if
    if (env%calc%nfreeze < 0) then
      message = 'process MTD freeze count is negative'
      return
    end if
    if (allocated(env%calc%freezelist)) then
      if (size(env%calc%freezelist) /= nat .or. &
      & count(env%calc%freezelist) /= env%calc%nfreeze) then
        message = 'process MTD freeze mask/count is inconsistent with molecule atom count'
        return
      end if
    else if (env%calc%nfreeze /= 0) then
      message = 'process MTD freeze count is nonzero without a freeze mask'
      return
    end if
    if (allocated(env%calc%calcs(1)%gff_fragments)) then
      call gfnff_fragment_ids(env%calc%calcs(1),mols(1),fragment_ids,fragment_io,message)
      if (fragment_io /= 0) then
        message = 'process MTD GFN-FF fragment selections are invalid: '//trim(message)
        return
      end if
      deallocate(fragment_ids)
    end if
    do i = 1,nsim
      if (mols(i)%nat /= nat .or. .not.allocated(mols(i)%at) .or. &
      & .not.allocated(mols(i)%xyz)) then
        message = 'process MTD prepared molecule atom count differs across trajectories'
        return
      end if
      if (size(mols(i)%at) /= nat .or. size(mols(i)%xyz,1) /= 3 .or. &
      & size(mols(i)%xyz,2) /= nat .or. .not.all(ieee_is_finite(mols(i)%xyz))) then
        message = 'prepared process MTD molecule has invalid dimensions or coordinates'
        return
      end if
      if (mols(i)%pdb%nat /= 0 .or. mols(i)%pdb%frag /= 0 .or. &
      & allocated(mols(i)%pdb%athet) .or. allocated(mols(i)%pdb%pdbat) .or. &
      & allocated(mols(i)%pdb%pdbas) .or. allocated(mols(i)%pdb%pdbfrag) .or. &
      & allocated(mols(i)%pdb%pdbgrp) .or. allocated(mols(i)%pdb%pdbocc) .or. &
      & allocated(mols(i)%pdb%pdbtf)) then
        message = 'process MTD prototype rejects nonempty PDB metadata'
        return
      end if
      if (i > 1 .and. any(mols(i)%at /= mols(1)%at)) then
        message = 'process MTD molecule atom ordering differs across jobs'
        return
      end if
      if (.not.mddats(i)%requested .or. mddats(i)%simtype /= type_mtd .or. &
      & mddats(i)%md_index /= i .or. mddats(i)%samerand) then
        message = 'resolved MTD job metadata is invalid or requests samerand'
        return
      end if
      if (mddats(i)%restart) then
        message = 'process MTD v1 rejects MD restart input; only fresh stochastic MTD is supported'
        return
      end if
      if (mddats(i)%shake .or. mddats(i)%nshake /= 0 .or. &
      & mddats(i)%shk%initialized .or. mddats(i)%shk%shake_mode /= 0 .or. &
      & mddats(i)%shk%nusr /= 0 .or. mddats(i)%shk%ncons /= 0 .or. &
      & allocated(mddats(i)%shk%conslistu) .or. allocated(mddats(i)%shk%wbo) .or. &
      & allocated(mddats(i)%shk%conslist) .or. allocated(mddats(i)%shk%distcons) .or. &
      & allocated(mddats(i)%shk%dro) .or. allocated(mddats(i)%shk%dr) .or. &
      & allocated(mddats(i)%shk%xyzt) .or. associated(mddats(i)%shk%freezeptr)) then
        message = 'process MTD does not yet serialize active SHAKE state'
        return
      end if
      if (mddats(i)%length_steps <= 1 .or. mddats(i)%sdump <= 0) then
        message = 'process MTD requires positive fresh step and dump counts'
        return
      end if
      if ((mddats(i)%length_steps-1)/mddats(i)%sdump < 1 .or. &
      & mddats(i)%dumped /= 0) then
        message = 'process MTD requires a fresh run with a positive exact trajectory frame count'
        return
      end if
      if (allocated(mddats(i)%blockrege) .or. allocated(mddats(i)%blocke) .or. &
      & allocated(mddats(i)%blockt)) then
        message = 'process MTD rejects preallocated dynamics block scratch arrays'
        return
      end if
      if (mddats(i)%termination_status /= -1 .or. mddats(i)%iblock /= 0 .or. &
      & mddats(i)%nblock /= 0 .or. mddats(i)%blocknreg /= 0) then
        message = 'process MTD requires canonical fresh termination and block counters'
        return
      end if
      if (mddats(i)%npot /= 1 .or. .not.allocated(mddats(i)%mtd) .or. &
      & .not.allocated(mddats(i)%cvtype)) then
        message = 'process MTD requires one already-resolved RMSD bias per job'
        return
      end if
      if (size(mddats(i)%mtd) /= 1 .or. size(mddats(i)%cvtype) /= 1) then
        message = 'process MTD potential/CV allocation count is not one'
        return
      end if
      if (mddats(i)%cvtype(1) /= cv_rmsd .or. mddats(i)%mtd(1)%mtdtype /= cv_rmsd) then
        message = 'process MTD supports resolved RMSD metadynamics only'
        return
      end if
      if (.not.ieee_is_finite(mddats(i)%mtd(1)%com_factor) .or. &
      & .not.ieee_is_finite(mddats(i)%mtd(1)%com_width)) then
        message = 'process MTD COM bias settings contain a nonfinite value'
        return
      end if
      if (mddats(i)%mtd(1)%com_bias .and. &
      & (mddats(i)%mtd(1)%com_factor <= 0.0_wp .or. &
      &  mddats(i)%mtd(1)%com_width <= 0.0_wp)) then
        message = 'enabled process MTD COM bias has a nonpositive factor or width'
        return
      end if
      if (allocated(mddats(i)%mtd(1)%atinclude)) then
        if (size(mddats(i)%mtd(1)%atinclude) /= nat) then
          message = 'process MTD inclusion mask length differs from molecule atom count'
          return
        end if
      end if
      if (mddats(i)%mtd(1)%nmax /= 0 .or. mddats(i)%mtd(1)%ncur /= 0 .or. &
      & mddats(i)%mtd(1)%cvdump /= 0 .or. mddats(i)%mtd(1)%cvdumpstep /= 0 .or. &
      & mddats(i)%mtd(1)%maxsave /= 0 .or. allocated(mddats(i)%mtd(1)%cv) .or. &
      & allocated(mddats(i)%mtd(1)%cvgrd) .or. allocated(mddats(i)%mtd(1)%cvxyz) .or. &
      & allocated(mddats(i)%mtd(1)%biasfile) .or. allocated(mddats(i)%mtd(1)%damping)) then
        message = 'process MTD requires fresh, already-resolved RMSD bias settings without history'
        return
      end if
      if (mddats(i)%input_structure_id < 0 .or. mddats(i)%bias_configuration_id < 0) then
        message = 'process MTD input/bias metadata contains a negative identifier'
        return
      end if
    end do
    status = status_normal
  end subroutine validate_process_configuration

  subroutine verify_worker_openmp_binding(cpu_list,expected_threads,io,message)
    character(len=*),intent(in) :: cpu_list
    integer,intent(in) :: expected_threads
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    integer :: team_size,thread_id,i,j,read_io,declared_cpu
    integer,allocatable :: thread_places(:),expected_cpus(:),place_cpus(:)
    integer(c_int),allocatable :: thread_cpus(:),affinity_status(:), &
    & affinity_counts(:),affinity_cpu_ids(:,:)
    character(len=256) :: expected_set,observed_set

    io = 0
    message = ''
    if (expected_threads < 1) then
      io = 1
      message = 'worker OpenMP verification received a nonpositive thread count'
      return
    end if
    allocate(thread_places(expected_threads),thread_cpus(expected_threads), &
    & expected_cpus(expected_threads),place_cpus(expected_threads), &
    & affinity_status(expected_threads),affinity_counts(expected_threads), &
    & affinity_cpu_ids(affinity_report_limit,expected_threads))
    team_size = 0
    thread_places = -1
    thread_cpus = -1
    place_cpus = -1
    affinity_status = -1
    affinity_counts = -1
    affinity_cpu_ids = -1
    read(cpu_list,*,iostat=read_io) expected_cpus
    if (read_io /= 0) then
      io = 1
      message = 'cannot parse the expected worker CPU IDs'
      return
    end if
    if (omp_get_num_places() /= expected_threads) then
      io = 2
      call format_cpu_set(expected_cpus,expected_threads,expected_set)
      write(message,'(a,i0,a,i0,a,a,a)') 'libgomp exposed ',omp_get_num_places(), &
      & ' places but worker requires ',expected_threads,'; expected CPU set=', &
      & trim(expected_set),'; observed place/CPU set=unavailable'
      return
    end if
    do i = 1,expected_threads
      if (omp_get_place_num_procs(i-1) /= 1) then
        io = 3
        message = 'an OpenMP place is not exactly one logical CPU'
        return
      end if
      call omp_get_place_proc_ids(i-1,place_cpus(i:i))
    end do
    !$omp parallel default(none) shared(team_size,thread_places,thread_cpus,expected_threads, &
    !$omp& place_cpus,affinity_status,affinity_counts,affinity_cpu_ids) private(thread_id)
    thread_id = omp_get_thread_num()
    !$omp single
    team_size = omp_get_num_threads()
    !$omp end single
    if (thread_id >= 0 .and. thread_id < expected_threads) then
      thread_places(thread_id+1) = omp_get_place_num()
      affinity_status(thread_id+1) = c_get_thread_affinity( &
      & int(affinity_report_limit,c_int),affinity_counts(thread_id+1), &
      & affinity_cpu_ids(:,thread_id+1))
      thread_cpus(thread_id+1) = c_current_cpu()
    end if
    !$omp end parallel
    if (team_size /= expected_threads .or. any(thread_places < 0)) then
      io = 4
      call format_cpu_set(expected_cpus,expected_threads,expected_set)
      call format_nonnegative_cpu_set(int(thread_cpus),observed_set)
      write(message,'(a,i0,a,i0,a,a,a,a)') 'worker formed ',team_size, &
      & ' OpenMP threads but expected ',expected_threads,'; expected CPU set=', &
      & trim(expected_set),'; observed CPU samples=',trim(observed_set)
      return
    end if
    do i = 1,expected_threads
      if (thread_places(i) < 0 .or. thread_places(i) >= expected_threads) then
        io = 5
        call worker_affinity_failure_message(message,'out-of-range OpenMP place',i-1, &
        & thread_places(i),-1,thread_cpus(i),affinity_status(i),affinity_counts(i), &
        & affinity_cpu_ids(:,i))
        return
      end if
      declared_cpu = place_cpus(thread_places(i)+1)
      if (affinity_status(i) /= 0_c_int) then
        io = 5
        call worker_affinity_failure_message(message,'sched_getaffinity failed',i-1, &
        & thread_places(i),declared_cpu,thread_cpus(i),affinity_status(i), &
        & affinity_counts(i),affinity_cpu_ids(:,i))
        return
      end if
      if (affinity_counts(i) /= 1_c_int) then
        io = 5
        call worker_affinity_failure_message(message,'kernel affinity mask is not singleton', &
        & i-1,thread_places(i),declared_cpu,thread_cpus(i),affinity_status(i), &
        & affinity_counts(i),affinity_cpu_ids(:,i))
        return
      end if
      if (affinity_cpu_ids(1,i) /= int(declared_cpu,c_int)) then
        io = 5
        call worker_affinity_failure_message(message, &
        & 'kernel affinity mask differs from declared singleton place',i-1, &
        & thread_places(i),declared_cpu,thread_cpus(i),affinity_status(i), &
        & affinity_counts(i),affinity_cpu_ids(:,i))
        return
      end if
      if (thread_cpus(i) < 0_c_int) then
        io = 5
        call worker_affinity_failure_message(message,'sched_getcpu failed',i-1, &
        & thread_places(i),declared_cpu,thread_cpus(i),affinity_status(i), &
        & affinity_counts(i),affinity_cpu_ids(:,i))
        return
      end if
      if (thread_cpus(i) /= int(declared_cpu,c_int)) then
        io = 5
        call worker_affinity_failure_message(message, &
        & 'worker thread is not executing on its declared singleton OpenMP place',i-1, &
        & thread_places(i),declared_cpu,thread_cpus(i),affinity_status(i), &
        & affinity_counts(i),affinity_cpu_ids(:,i))
        return
      end if
      if (.not.any(declared_cpu == expected_cpus)) then
        io = 6
        call worker_affinity_failure_message(message, &
        & 'declared singleton CPU is outside the scheduler assignment',i-1, &
        & thread_places(i),declared_cpu,thread_cpus(i),affinity_status(i), &
        & affinity_counts(i),affinity_cpu_ids(:,i))
        return
      end if
      do j = i+1,expected_threads
        if (thread_places(i) == thread_places(j)) then
          io = 7
          call format_cpu_set(expected_cpus,expected_threads,expected_set)
          call format_nonnegative_cpu_set(int(thread_cpus),observed_set)
          write(message,'(a,i0,a,i0,a,i0,a,a,a,a)') 'worker OpenMP threads ',i-1, &
          & ' and ',j-1,' share place ',thread_places(i),'; expected CPU set=', &
          & trim(expected_set),'; observed CPU set=',trim(observed_set)
          return
        end if
        if (thread_cpus(i) == thread_cpus(j)) then
          io = 8
          call format_cpu_set(expected_cpus,expected_threads,expected_set)
          call format_nonnegative_cpu_set(int(thread_cpus),observed_set)
          write(message,'(a,i0,a,i0,a,i0,a,a,a,a)') 'worker OpenMP threads ',i-1, &
          & ' and ',j-1,' sampled CPU ',int(thread_cpus(i)),'; expected CPU set=', &
          & trim(expected_set),'; observed CPU set=',trim(observed_set)
          return
        end if
      end do
    end do
    do i = 1,expected_threads
      if (.not.any(place_cpus == expected_cpus(i))) then
        io = 9
        call format_cpu_set(expected_cpus,expected_threads,expected_set)
        call format_cpu_set(place_cpus,expected_threads,observed_set)
        message = 'OpenMP place/CPU set mismatch: expected='//trim(expected_set)// &
        & '; observed='//trim(observed_set)
        return
      end if
    end do
  end subroutine verify_worker_openmp_binding

  subroutine worker_affinity_failure_message(message,reason,thread_id,place,declared_cpu, &
  & current_cpu,affinity_errno,affinity_count,affinity_cpu_ids)
    character(len=*),intent(out) :: message
    character(len=*),intent(in) :: reason
    integer,intent(in) :: thread_id,place,declared_cpu
    integer(c_int),intent(in) :: current_cpu,affinity_errno,affinity_count
    integer(c_int),intent(in) :: affinity_cpu_ids(:)
    character(len=256) :: observed_mask,expected_mask
    character(len=64) :: value

    message = 'worker affinity verification failed: '//trim(reason)
    write(value,'(i0)') thread_id
    call append_bounded(message,'; thread='//trim(value))
    write(value,'(i0)') place
    call append_bounded(message,'; place='//trim(value))
    if (declared_cpu >= 0) then
      write(value,'(i0)') declared_cpu
      call append_bounded(message,'; declared CPU='//trim(value))
      write(expected_mask,'(a,i0,a)') '{',declared_cpu,'}'
    else
      call append_bounded(message,'; declared CPU=unavailable')
      expected_mask = '{unavailable}'
    end if
    if (current_cpu >= 0_c_int) then
      write(value,'(i0)') int(current_cpu)
      call append_bounded(message,'; current CPU='//trim(value))
    else
      call append_bounded(message,'; current CPU=unavailable')
      write(value,'(i0)') int(-current_cpu)
      call append_bounded(message,'; sched_getcpu errno='//trim(value))
    end if
    if (affinity_errno == 0_c_int) then
      call format_cpu_set(int(affinity_cpu_ids),int(affinity_count),observed_mask)
    else
      observed_mask = '{unavailable}'
    end if
    call append_bounded(message,'; observed mask='//trim(observed_mask))
    call append_bounded(message,'; expected mask='//trim(expected_mask))
    if (affinity_errno > 0_c_int) then
      write(value,'(i0)') int(affinity_errno)
      call append_bounded(message,'; sched_getaffinity errno='//trim(value))
    end if
  end subroutine worker_affinity_failure_message

  subroutine format_nonnegative_cpu_set(cpus,text)
    integer,intent(in) :: cpus(:)
    character(len=*),intent(out) :: text
    integer :: compact(size(cpus)),i,count

    compact = -1
    count = 0
    do i = 1,size(cpus)
      if (cpus(i) >= 0) then
        count = count+1
        compact(count) = cpus(i)
      end if
    end do
    call format_cpu_set(compact,count,text)
  end subroutine format_nonnegative_cpu_set

  subroutine format_cpu_set(cpus,count,text)
    integer,intent(in) :: cpus(:),count
    character(len=*),intent(out) :: text
    character(len=32) :: value
    integer :: i,shown

    text = ''
    if (count < 0) then
      text = '{unavailable}'
      return
    end if
    call append_bounded(text,'{')
    shown = min(count,size(cpus),affinity_report_limit)
    do i = 1,shown
      if (i > 1) call append_bounded(text,',')
      write(value,'(i0)') cpus(i)
      call append_bounded(text,trim(value))
    end do
    if (count > shown) call append_bounded(text,',...')
    call append_bounded(text,'}')
  end subroutine format_cpu_set

  subroutine append_bounded(text,piece)
    character(len=*),intent(inout) :: text
    character(len=*),intent(in) :: piece
    integer :: available,ncopy,start

    start = len_trim(text)+1
    available = len(text)-start+1
    ncopy = min(available,len_trim(piece))
    if (ncopy > 0) text(start:start+ncopy-1) = piece(1:ncopy)
  end subroutine append_bounded

  subroutine copy_file_exact(source,destination,io,message)
    character(len=*),intent(in) :: source,destination
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    integer :: input_unit,output_unit,n
    integer(int64) :: source_size,position,remaining
    integer(int8),allocatable :: buffer(:)
    integer,parameter :: chunk = 1048576
    message = ''
    io = 0
    source_size = 0_int64
    inquire(file=source,size=source_size,iostat=io,iomsg=message)
    if (io /= 0 .or. source_size < 1) then
      if (io == 0) io = 1
      return
    end if
    open(newunit=input_unit,file=source,access='stream',form='unformatted', &
    & status='old',action='read',iostat=io,iomsg=message)
    if (io /= 0) return
    open(newunit=output_unit,file=destination,access='stream',form='unformatted', &
    & status='new',action='write',iostat=io,iomsg=message)
    if (io /= 0) then
      close(input_unit)
      return
    end if
    allocate(buffer(chunk))
    position = 1_int64
    remaining = source_size
    do while (remaining > 0)
      n = int(min(remaining,int(chunk,int64)))
      read(input_unit,pos=position,iostat=io,iomsg=message) buffer(1:n)
      if (io /= 0) exit
      write(output_unit,iostat=io,iomsg=message) buffer(1:n)
      if (io /= 0) exit
      position = position+int(n,int64)
      remaining = remaining-int(n,int64)
    end do
    close(input_unit)
    close(output_unit)
  end subroutine copy_file_exact

  subroutine write_worker_result(path,worker_index,termination,wall,engrad_calls,io,message)
    character(len=*),intent(in) :: path
    integer,intent(in) :: worker_index,termination
    real(wp),intent(in) :: wall
    integer(int64),intent(in) :: engrad_calls
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    integer :: unit,close_io
    message = ''
    open(newunit=unit,file=path,access='stream',form='unformatted',status='new', &
    & action='write',iostat=io,iomsg=message)
    if (io /= 0) return
    write(unit,iostat=io,iomsg=message) result_magic,result_version, &
    & int(worker_index,int32),int(termination,int32),wall,engrad_calls
    close(unit,iostat=close_io)
    if (io == 0) io = close_io
  end subroutine write_worker_result

  subroutine read_worker_result(path,expected_worker,termination,wall,engrad_calls,io,message)
    character(len=*),intent(in) :: path
    integer,intent(in) :: expected_worker
    integer,intent(out) :: termination
    real(wp),intent(out) :: wall
    integer(int64),intent(out) :: engrad_calls
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    character(len=32) :: magic
    integer(int32) :: version,worker,term
    integer :: unit
    integer(int64) :: file_size,expected_size
    message = ''
    termination = -1
    wall = -1.0_wp
    engrad_calls = -1_int64
    expected_size = int(len(result_magic),int64)+ &
    & 3_int64*int(storage_size(result_version)/8,int64)+ &
    & int(storage_size(wall)/8,int64)+int(storage_size(engrad_calls)/8,int64)
    inquire(file=path,size=file_size,iostat=io,iomsg=message)
    if (io /= 0) return
    if (file_size /= expected_size) then
      io = 1
      message = 'worker result byte size is not canonical'
      return
    end if
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old', &
    & action='read',iostat=io,iomsg=message)
    if (io /= 0) return
    read(unit,iostat=io,iomsg=message) magic,version,worker,term,wall,engrad_calls
    close(unit)
    if (io /= 0) return
    if (magic /= result_magic .or. version /= result_version .or. &
    & worker /= expected_worker .or. .not.ieee_is_finite(wall) .or. wall < 0.0_wp .or. &
    & engrad_calls < 0_int64) then
      io = 2
      message = 'worker result invariant mismatch'
      return
    end if
    termination = int(term)
  end subroutine read_worker_result

  subroutine validate_xyz_trajectory(path,expected_types,expected_frames,frames, &
  & file_bytes,io,message)
    character(len=*),intent(in) :: path
    integer,intent(in) :: expected_types(:),expected_frames
    integer,intent(out) :: frames
    integer(int64),intent(out) :: file_bytes
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    integer :: unit,natoms,atom,comment_io
    real(wp) :: x,y,z,energy
    character(len=32) :: symbol,energy_label,equals
    character(len=4096) :: line
    frames = 0
    file_bytes = 0_int64
    io = 0
    message = ''
    inquire(file=path,size=file_bytes,iostat=io,iomsg=message)
    if (io /= 0 .or. file_bytes < 1) then
      if (io == 0) then
        io = 1
        message = 'trajectory is empty'
      end if
      return
    end if
    open(newunit=unit,file=path,status='old',action='read',form='formatted', &
    & iostat=io,iomsg=message)
    if (io /= 0) return
    do
      read(unit,*,iostat=io,iomsg=message) natoms
      if (io == iostat_end) then
        io = 0
        exit
      else if (io /= 0) then
        exit
      end if
      if (natoms /= size(expected_types)) then
        io = 2
        message = 'trajectory frame atom count differs from prepared molecule'
        exit
      end if
      read(unit,'(a)',iostat=io,iomsg=message) line
      if (io /= 0) exit
      energy_label = ''
      equals = ''
      energy = 0.0_wp
      read(line,*,iostat=comment_io) energy_label,equals,energy
      if (comment_io /= 0) then
        io = 3
        message = 'trajectory frame has a malformed Epot comment'
        exit
      end if
      if (trim(energy_label) /= 'Epot' .or. trim(equals) /= '=' .or. &
      & .not.ieee_is_finite(energy)) then
        io = 3
        message = 'trajectory frame lacks a finite canonical Epot comment'
        exit
      end if
      do atom = 1,natoms
        read(unit,'(a)',iostat=io,iomsg=message) line
        if (io /= 0) exit
        symbol = ''
        x = 0.0_wp
        y = 0.0_wp
        z = 0.0_wp
        read(line,*,iostat=io) symbol,x,y,z
        if (io /= 0) then
          io = 4
          message = 'trajectory contains a malformed atom record'
          exit
        end if
        if (trim(symbol) /= trim(i2e(expected_types(atom),'nc')) .or. &
        & .not.all([ieee_is_finite(x),ieee_is_finite(y),ieee_is_finite(z)])) then
          io = 4
          message = 'trajectory contains malformed or nonfinite coordinates'
          exit
        end if
      end do
      if (io /= 0) exit
      frames = frames+1
    end do
    close(unit)
    if (io == 0 .and. frames /= expected_frames) then
      io = 5
      write(message,'(a,i0,a,i0)') 'trajectory frame count ',frames, &
      & ' differs from exact expected count ',expected_frames
    end if
  end subroutine validate_xyz_trajectory

  subroutine scan_worker_logs(workdir,io,message)
    character(len=*),intent(in) :: workdir
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    integer :: i,normal_count,total_normal
    character(len=:),allocatable :: path
    character(len=13),parameter :: names(2) = [character(len=13) :: &
    & 'worker.stdout','worker.stderr']
    io = 0
    message = ''
    total_normal = 0
    do i = 1,size(names)
      path = trim(workdir)//'/'//trim(names(i))
      call scan_one_worker_log(path,normal_count,io,message)
      if (io /= 0) return
      total_normal = total_normal+normal_count
    end do
    if (total_normal /= 1) then
      io = 2
      write(message,'(a,i0)') 'worker logs contain unexpected normal-termination marker count ', &
      & total_normal
    end if
  end subroutine scan_worker_logs

  subroutine scan_one_worker_log(path,normal_count,io,message)
    character(len=*),intent(in) :: path
    integer,intent(out) :: normal_count
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    character(len=48),parameter :: fatal(42) = [character(len=48) :: &
    & 'program received signal','sigsegv','segmentation fault', &
    & 'invalid memory reference','backtrace for this error', &
    & 'bad memory unallocation','openblas warning','precompiled num_threads exceeded', &
    & 'fortran runtime error','error stop','**error**','crest terminated abnormally', &
    & 'crest terminated with','crest terminated due to','sigbus','bus error','sigabrt', &
    & 'double free','invalid pointer','corrupted size','stack smashing','out of memory', &
    & 'cannot allocate memory','floating point exception','error in md calculation', &
    & 'command terminated by signal','command exited with non-zero status', &
    & 'error while loading shared libraries','symbol lookup error','factorisation failed', &
    & 'solving linear equations failed','solving linear system failed', &
    & 'failed to generate topology','failed to generate charges','oom-kill','core dumped', &
    & 'munmap_chunk','ieee_invalid_flag','ieee_divide_by_zero','ieee_overflow_flag', &
    & 'fatal runtime','fatal error']
    character(len=4096) :: line,lower
    integer :: unit,line_number,i,read_io
    io = 0
    message = ''
    normal_count = 0
    open(newunit=unit,file=path,status='old',action='read',form='formatted', &
    & iostat=io,iomsg=message)
    if (io /= 0) return
    line_number = 0
    do
      read(unit,'(a)',iostat=read_io,iomsg=message) line
      if (read_io == iostat_end) exit
      if (read_io /= 0) then
        io = read_io
        close(unit)
        return
      end if
      line_number = line_number+1
      lower = ascii_lower(line)
      if (index(lower,'crest terminated normally') > 0) normal_count = normal_count+1
      if (contains_nan_token(lower)) then
        io = 1
        write(message,'(a,i0,a,a)') 'NaN token at line ',line_number,' of ',trim(path)
        close(unit)
        return
      end if
      do i = 1,size(fatal)
        if (index(lower,trim(fatal(i))) > 0) then
          io = 1
          write(message,'(a,a,a,i0,a,a)') 'fatal signature "',trim(fatal(i)), &
          & '" at line ',line_number,' of ',trim(path)
          close(unit)
          return
        end if
      end do
    end do
    close(unit)
  end subroutine scan_one_worker_log

  pure function ascii_lower(value) result(lower)
    character(len=*),intent(in) :: value
    character(len=len(value)) :: lower
    integer :: i,code
    lower = value
    do i = 1,len(value)
      code = iachar(lower(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code+32)
    end do
  end function ascii_lower

  pure logical function contains_nan_token(value) result(found)
    character(len=*),intent(in) :: value
    integer :: position,before,after,relative
    found = .false.
    position = index(value,'nan')
    do while (position > 0)
      before = 0
      after = 0
      if (position > 1) before = iachar(value(position-1:position-1))
      if (position+3 <= len(value)) after = iachar(value(position+3:position+3))
      if ((position == 1 .or. .not.is_ascii_letter(before)) .and. &
      & (position+2 == len(value) .or. .not.is_ascii_letter(after))) then
        found = .true.
        return
      end if
      if (position+3 > len(value)) exit
      relative = index(value(position+3:),'nan')
      if (relative == 0) exit
      position = position+2+relative
    end do
  end function contains_nan_token

  pure logical function is_ascii_letter(code) result(letter)
    integer,intent(in) :: code
    letter = (code >= iachar('a') .and. code <= iachar('z')) .or. &
    & (code >= iachar('A') .and. code <= iachar('Z'))
  end function is_ascii_letter

  subroutine terminate_and_reap_workers(pids,reaped,count)
    integer(c_int),intent(in) :: pids(:)
    logical,intent(inout) :: reaped(:)
    integer,intent(in) :: count
    integer(c_int) :: ignored,exit_code,term_signal
    integer :: i
    do i = 1,min(count,size(pids))
      if (pids(i) > 0_c_int .and. .not.reaped(i)) ignored = c_kill_worker(pids(i))
    end do
    do i = 1,min(count,size(pids))
      if (pids(i) > 0_c_int .and. .not.reaped(i)) then
        ignored = c_reap_worker_bounded(pids(i),100_c_int,exit_code,term_signal)
        reaped(i) = ignored == 0_c_int
      end if
    end do
  end subroutine terminate_and_reap_workers

  subroutine print_process_finish(md,index,wall,frames,file_bytes)
    type(mddata),intent(in) :: md
    integer,intent(in) :: index,frames
    real(wp),intent(in) :: wall
    integer(int64),intent(in) :: file_bytes
    integer(int64) :: minutes
    real(wp) :: seconds
    minutes = int(wall/60.0_wp,int64)
    seconds = wall-real(minutes,wp)*60.0_wp
    write(stdout,'(a,i3,a,i9,a,f6.3,a)') '*MTD ',index, &
    & ' completed successfully ...',minutes,' min, ',seconds,' sec'
    write(stdout,'(2x,a,i0,a,i0,a,i0,a,i0)') 'process worker metadata: input=', &
    & md%input_structure_id,', bias=',md%bias_configuration_id,', frames=',frames, &
    & ', bytes=',file_bytes
    flush(stdout)
  end subroutine print_process_finish

  subroutine to_c_string(value,cvalue)
    character(len=*),intent(in) :: value
    character(kind=c_char),allocatable,intent(out) :: cvalue(:)
    integer :: i,n
    n = len_trim(value)
    allocate(cvalue(n+1))
    do i = 1,n
      cvalue(i) = value(i:i)
    end do
    cvalue(n+1) = c_null_char
  end subroutine to_c_string

  subroutine secure_directory(path,io,message)
    character(len=*),intent(in) :: path
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    character(kind=c_char),allocatable :: c_path(:)
    integer(c_int) :: cstat
    call to_c_string(path,c_path)
    cstat = c_secure_directory(c_path)
    io = int(cstat)
    if (io == 0) then
      message = ''
    else
      write(message,'(a,i0)') 'secure-directory errno=',io
    end if
  end subroutine secure_directory

  function integer_string(value) result(text)
    integer,intent(in) :: value
    character(len=:),allocatable :: text
    character(len=64) :: buffer
    write(buffer,'(i0)') value
    text = trim(buffer)
  end function integer_string

  function zero_pad4(value) result(text)
    integer,intent(in) :: value
    character(len=4) :: text
    write(text,'(i4.4)') value
  end function zero_pad4

end module mtd_process_scheduler
