!> mlip_client.f90
!> ─────────────────────────────────────────────────────────────────────────────
!> Fortran module for managing one or more persistent Python MLIP server
!> processes.  Each "instance" maps to one Python process + one socket
!> connection.  Up to MLIP_MAX_INSTANCES can run in parallel.
!>
!> Public API
!> ----------
!>   mlip_init      (iid, port, server_cmd, timeout_sec, ierr)
!>       Spawn Python server, wait for "READY", open socket.
!>       server_cmd is the full shell command, e.g.:
!>         "fmlip-relay-server --port 54321 --backend mace --model /path/to/model.model"
!>
!>   mlip_compute   (iid, natoms, positions, cell, pbc, compute_stress, charge, spin,
!>                   energy, forces, stress, ierr)
!>       Send one COMPUTE request, receive energy/forces/stress.
!>       charge (total charge) and spin (multiplicity 2S+1) are optional and
!>       default to a neutral, ground-state singlet (charge=0, spin=1); they
!>       are only used by charge/spin-aware backends (e.g. UMA omol).
!>
!>   mlip_ping      (iid, ierr)
!>       Lightweight round-trip check.
!>
!>   mlip_finalize  (iid, ierr)
!>       Send QUIT, close socket, wait for child process to exit.
!>
!>   mlip_finalize_all (ierr)
!>       Finalise every active instance.
!>
!> Notes
!> -----
!> • All integers are default INTEGER (32-bit).
!> • Positions in Angstrom, energy in eV, forces in eV/Ang, stress in eV/Ang^3.
!> • Thread-safety: each instance owns independent socket fd; the module uses
!>   only instance-local state, so distinct instances may be used from
!>   distinct OpenMP threads without locking.
!> ─────────────────────────────────────────────────────────────────────────────

module fmlip_relay_client
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use iso_c_binding,only:c_int,c_double,c_char,c_null_char,c_loc
  implicit none
  private

  ! ── public symbols ──────────────────────────────────────────────────────────
  public :: mlip_init,mlip_compute,mlip_ping,mlip_finalize,mlip_finalize_all
  public :: MLIP_MAX_INSTANCES,MLIP_OK,MLIP_ERR,GET_PID

  ! ── constants ────────────────────────────────────────────────────────────────
  integer,parameter :: MLIP_MAX_INSTANCES = 16
  integer,parameter :: MLIP_OK = 0
  integer,parameter :: MLIP_ERR = 1

  ! wire protocol
  integer,parameter :: MSG_COMPUTE = 1
  integer,parameter :: MSG_QUIT = 2
  integer,parameter :: MSG_PING = 3
  integer,parameter :: STATUS_OK = 0

  ! ── C interfaces ─────────────────────────────────────────────────────────────
  interface
    function sock_connect_c(host,port) bind(C,name='sock_connect_c')
      import c_int,c_char
      character(kind=c_char),intent(in) :: host(*)
      integer(c_int),value      :: port
      integer(c_int)                     :: sock_connect_c
    end function

    function sock_send_c(fd,buf,nbytes) bind(C,name='sock_send_c')
      import c_int
      integer(c_int),value         :: fd,nbytes
      integer(c_int),intent(in)    :: buf(*)   ! raw bytes
      integer(c_int)                :: sock_send_c
    end function

    function sock_recv_c(fd,buf,nbytes) bind(C,name='sock_recv_c')
      import c_int
      integer(c_int),value         :: fd,nbytes
      integer(c_int),intent(inout) :: buf(*)
      integer(c_int)                :: sock_recv_c
    end function

    subroutine sock_close_c(fd) bind(C,name='sock_close_c')
      import c_int
      integer(c_int),value :: fd
    end subroutine

    function spawn_process_c(cmd,stdout_fd) bind(C,name='spawn_process_c')
      import c_int,c_char
      character(kind=c_char),intent(in) :: cmd(*)
      integer(c_int),intent(out)        :: stdout_fd
      integer(c_int)                     :: spawn_process_c
    end function

    function wait_for_line_c(fd,prefix,timeout_sec) bind(C,name='wait_for_line_c')
      import c_int,c_char
      integer(c_int),value              :: fd,timeout_sec
      character(kind=c_char),intent(in) :: prefix(*)
      integer(c_int)                     :: wait_for_line_c
    end function

    subroutine kill_process_c(pid) bind(C,name='kill_process_c')
      import c_int
      integer(c_int),value :: pid
    end subroutine

    subroutine usleep(useconds) bind(C,name='usleep')
      import c_int
      integer(c_int),value :: useconds
    end subroutine

    function get_pid_shim() bind(C,name="get_pid_shim") result(pid)
      import c_int
      integer(c_int) :: pid
    end function get_pid_shim
  end interface

  ! ── instance state ────────────────────────────────────────────────────────────
  type :: mlip_instance_t
    logical  :: active = .false.
    integer  :: sock_fd = -1      ! TCP socket file descriptor
    integer  :: proc_pid = -1     ! child PID
    integer  :: stdout_fd = -1    ! pipe to child stdout (used at startup)
    integer  :: port = -1
  end type

  type(mlip_instance_t),save :: instances(MLIP_MAX_INSTANCES)

contains

  ! ════════════════════════════════════════════════════════════════════════════
  !> mlip_init – spawn server and connect
  ! ════════════════════════════════════════════════════════════════════════════
  subroutine mlip_init(iid,port,server_cmd,timeout_sec,ierr)
    integer,intent(in)  :: iid          !< instance id  [1..MLIP_MAX_INSTANCES]
    integer,intent(in)  :: port         !< TCP port for this instance
    character(len=*),intent(in)  :: server_cmd   !< full shell command to launch the server
    !<   e.g. "fmlip-relay-server --port 54321 --backend mace --model /path/model.model"
    integer,intent(in)  :: timeout_sec  !< seconds to wait for READY
    integer,intent(out) :: ierr

    character(kind=c_char,len=1024) :: cmd_c
    integer(c_int) :: pid_c,stdout_fd_c,rc,retries
    integer :: retry_delay_us = 100000  ! 0.1 s

    ierr = MLIP_ERR

    if (iid < 1.or.iid > MLIP_MAX_INSTANCES) then
      write (stdout,'(A,I0)') "mlip_init: invalid iid=",iid; return
    end if
    if (instances(iid)%active) then
      write (stdout,'(A,I0)') "mlip_init: instance already active, iid=",iid; return
    end if

    ! ── convert server_cmd to null-terminated C string ──────────────────────
    call f_to_c_str(server_cmd,cmd_c)

    ! ── spawn process ────────────────────────────────────────────────────────
    pid_c = spawn_process_c(cmd_c,stdout_fd_c)
    if (pid_c < 0) then
      write (stdout,'(A,I0,A)') "[mlip_init] instance ",iid,": failed to spawn process"; return
    end if
    instances(iid)%proc_pid = int(pid_c)
    instances(iid)%stdout_fd = int(stdout_fd_c)
    instances(iid)%port = port

    ! ── wait for "READY" on stdout ────────────────────────────────────────────
    rc = wait_for_line_c(stdout_fd_c,c_str("READY"),int(timeout_sec,c_int))
    if (rc /= 0) then
      write (stdout,'(A,I0,A)') "[mlip_init] instance ",iid,": timeout waiting for READY signal"
      call kill_process_c(pid_c)
      return
    end if

    ! ── connect socket (retry briefly in case server isn't listening yet) ────
    retries = 0
    instances(iid)%sock_fd = -1
    do while (instances(iid)%sock_fd < 0.and.retries < 10)
      instances(iid)%sock_fd = int( &
                               sock_connect_c(c_str("127.0.0.1"),int(port,c_int)))
      if (instances(iid)%sock_fd < 0) then
        call usleep_f(retry_delay_us)
        retries = retries+1
      end if
    end do
    if (instances(iid)%sock_fd < 0) then
      write (stdout,'(A,I0,A)') "[mlip_init] instance ",iid,": failed to connect to server"; return
    end if

    instances(iid)%active = .true.
    ierr = MLIP_OK
    write (stdout,'(A,I0,A,I0)') "[mlip_init] instance ",iid,": ready on port ",port
  end subroutine mlip_init

  ! ════════════════════════════════════════════════════════════════════════════
  !> mlip_compute – send coordinates, receive energy/forces/stress
  ! ════════════════════════════════════════════════════════════════════════════
  subroutine mlip_compute(iid,natoms,atomic_numbers,positions,cell,pbc,compute_stress, &
     &                    charge,spin, &
     &                    energy,forces,stress,ierr)
    integer,intent(in)   :: iid
    integer,intent(in)   :: natoms
    integer,intent(in)   :: atomic_numbers(natoms)  !< Z for each atom (e.g. 13=Al, 8=O)
    real(wp),intent(in)  :: positions(3,natoms)     !< (xyz, atom), Angstrom
    real(wp),intent(in)  :: cell(3,3)               !< row vectors
    integer,intent(in)   :: pbc(3)                  !< 1=periodic, 0=not
    integer,intent(in)   :: compute_stress          !< 1 to request stress
    integer,intent(in),optional :: charge           !< total molecular charge (default 0)
    integer,intent(in),optional :: spin             !< spin multiplicity 2S+1 (default 1, ground state)
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: forces(3,natoms)
    real(wp),intent(out) :: stress(3,3)
    integer,intent(out)  :: ierr

    integer(c_int) :: hdr(2),pbc_c(3),stress_flag_c(1),status(1),rc
    integer(c_int) :: atomic_numbers_c(natoms),charge_c(1),spin_c(1)
    real(wp)       :: pos_buf(3*natoms),cell_buf(9),stress_buf(9),energy_buf(1)
    integer        :: charge_l,spin_l

    ierr = MLIP_ERR
    energy = 0.0_wp
    forces = 0.0_wp
    stress = 0.0_wp

    ! ── resolve optional charge/spin: neutral, ground-state singlet default ──
    charge_l = 0
    if (present(charge)) charge_l = charge
    spin_l = 1
    if (present(spin)) spin_l = spin
    if (spin_l < 1) spin_l = 1   ! multiplicity 2S+1 is >= 1; clamp stray 0

    if (.not.instances(iid)%active) then
      write (stdout,'(A,I0)') "mlip_compute: instance not active, iid=",iid; return
    end if

    ! ── pack positions (C-order: atom-major) ────────────────────────────────
    pos_buf = reshape(positions, [3*natoms])
    cell_buf = reshape(cell, [9])
    pbc_c = int(pbc,c_int)
    atomic_numbers_c = int(atomic_numbers,c_int)
    stress_flag_c(1) = int(compute_stress,c_int)
    charge_c(1) = int(charge_l,c_int)
    spin_c(1) = int(spin_l,c_int)

    ! ── send header: msg_type, natoms ────────────────────────────────────────
    hdr(1) = int(MSG_COMPUTE,c_int)
    hdr(2) = int(natoms,c_int)
    rc = sock_send_bytes(iid,hdr,2*4); if (rc /= 0) return

    ! ── send atomic numbers (natoms int32) ───────────────────────────────────
    rc = sock_send_bytes(iid,atomic_numbers_c,natoms*4); if (rc /= 0) return

    ! ── send positions (3*natoms float64) ───────────────────────────────────
    rc = sock_send_r8(iid,pos_buf,3*natoms); if (rc /= 0) return

    ! ── send cell (9 float64) ────────────────────────────────────────────────
    rc = sock_send_r8(iid,cell_buf,9); if (rc /= 0) return

    ! ── send pbc (3 int32) ───────────────────────────────────────────────────
    rc = sock_send_bytes(iid,pbc_c,3*4); if (rc /= 0) return

    ! ── send stress flag (1 int32) ───────────────────────────────────────────
    rc = sock_send_bytes(iid,stress_flag_c,4); if (rc /= 0) return

    ! ── send molecular charge (1 int32) ──────────────────────────────────────
    rc = sock_send_bytes(iid,charge_c,4); if (rc /= 0) return

    ! ── send spin multiplicity (1 int32) ─────────────────────────────────────
    rc = sock_send_bytes(iid,spin_c,4); if (rc /= 0) return

    ! ── receive status (int32) ───────────────────────────────────────────────
    rc = sock_recv_bytes(iid,status,4); if (rc /= 0) return
    if (status(1) /= STATUS_OK) then
      write (stdout,*) "mlip_compute: server returned error status"; return
    end if

    ! ── receive energy (float64) ─────────────────────────────────────────────
    rc = sock_recv_r8(iid,energy_buf,1); if (rc /= 0) return
    energy = energy_buf(1)

    ! ── receive forces (3*natoms float64) ───────────────────────────────────
    rc = sock_recv_r8(iid,forces,3*natoms); if (rc /= 0) return

    ! ── receive stress (9 float64) ───────────────────────────────────────────
    rc = sock_recv_r8(iid,stress_buf,9); if (rc /= 0) return
    stress = reshape(stress_buf, [3,3])

    ierr = MLIP_OK
  end subroutine mlip_compute

  ! ════════════════════════════════════════════════════════════════════════════
  !> mlip_ping – check connectivity
  ! ════════════════════════════════════════════════════════════════════════════
  subroutine mlip_ping(iid,ierr)
    integer,intent(in)  :: iid
    integer,intent(out) :: ierr
    integer(c_int) :: msg(1),status(1)
    integer        :: rc

    ierr = MLIP_ERR
    if (.not.instances(iid)%active) return

    msg(1) = int(MSG_PING,c_int)
    rc = sock_send_bytes(iid,msg,4); if (rc /= 0) return
    rc = sock_recv_bytes(iid,status,4); if (rc /= 0) return
    if (status(1) == STATUS_OK) ierr = MLIP_OK
  end subroutine mlip_ping

  ! ════════════════════════════════════════════════════════════════════════════
  !> mlip_finalize – graceful shutdown of one instance
  ! ════════════════════════════════════════════════════════════════════════════
  subroutine mlip_finalize(iid,ierr)
    integer,intent(in)  :: iid
    integer,intent(out) :: ierr
    integer(c_int) :: msg(1)
    integer        :: rc

    ierr = MLIP_ERR
    if (.not.instances(iid)%active) then
      ierr = MLIP_OK; return    ! idempotent
    end if

    ! send QUIT
    msg(1) = int(MSG_QUIT,c_int)
    rc = sock_send_bytes(iid,msg,4)

    ! close socket
    call sock_close_c(int(instances(iid)%sock_fd,c_int))

    ! close stdout pipe
    if (instances(iid)%stdout_fd >= 0) &
      call sock_close_c(int(instances(iid)%stdout_fd,c_int))

    ! send SIGTERM just in case
    if (instances(iid)%proc_pid > 0) &
      call kill_process_c(int(instances(iid)%proc_pid,c_int))

    instances(iid) = mlip_instance_t()   ! reset to default
    write (stdout,'(A,I0,A)') "[mlip_finalize] instance ",iid,": shut down"
    ierr = MLIP_OK
  end subroutine mlip_finalize

  ! ════════════════════════════════════════════════════════════════════════════
  !> mlip_finalize_all
  ! ════════════════════════════════════════════════════════════════════════════
  subroutine mlip_finalize_all(ierr)
    integer,intent(out) :: ierr
    integer :: i,local_err
    ierr = MLIP_OK
    do i = 1,MLIP_MAX_INSTANCES
      if (instances(i)%active) then
        call mlip_finalize(i,local_err)
        if (local_err /= MLIP_OK) ierr = local_err
      end if
    end do
  end subroutine mlip_finalize_all

  ! ════════════════════════════════════════════════════════════════════════════
  !  Private helpers
  ! ════════════════════════════════════════════════════════════════════════════

  !> Send raw bytes via instance iid socket
  function sock_send_bytes(iid,buf,nbytes) result(rc)
    integer,intent(in) :: iid,nbytes
    integer(c_int),intent(in) :: buf(*)
    integer :: rc
    rc = int(sock_send_c(int(instances(iid)%sock_fd,c_int),buf, &
                         int(nbytes,c_int)))
  end function sock_send_bytes

  !> Receive raw bytes via instance iid socket
  function sock_recv_bytes(iid,buf,nbytes) result(rc)
    integer,intent(in)    :: iid,nbytes
    integer(c_int),intent(inout) :: buf(*)
    integer :: rc
    rc = int(sock_recv_c(int(instances(iid)%sock_fd,c_int),buf, &
                         int(nbytes,c_int)))
  end function sock_recv_bytes

  !> Send n float64 values
  function sock_send_r8(iid,arr,n) result(rc)
    integer,intent(in) :: iid,n
    real(wp),intent(in) :: arr(*)
    integer :: rc
    integer(c_int) :: buf(2*n)   ! 2 × int32 per float64
    buf = transfer(arr(1:n),buf)
    rc = int(sock_send_c(int(instances(iid)%sock_fd,c_int),buf,int(n*8,c_int)))
  end function sock_send_r8

  !> Receive n float64 values into arr(1:n)
  function sock_recv_r8(iid,arr,n) result(rc)
    integer,intent(in)    :: iid,n
    real(wp),intent(inout) :: arr(*)
    integer :: rc
    integer(c_int) :: buf(2*n)
    rc = int(sock_recv_c(int(instances(iid)%sock_fd,c_int),buf,int(n*8,c_int)))
    if (rc == 0) arr(1:n) = transfer(buf,arr(1:n))
  end function sock_recv_r8

  !> Convert Fortran string to null-terminated C string
  function c_str(s) result(cs)
    character(len=*),intent(in) :: s
    character(kind=c_char,len=len_trim(s)+1) :: cs
    cs = trim(s)//c_null_char
  end function c_str

  !> Copy Fortran string to fixed-length C char array
  subroutine f_to_c_str(s,cs)
    character(len=*),intent(in)  :: s
    character(kind=c_char,len=*),intent(out) :: cs
    integer :: i,n
    n = min(len_trim(s),len(cs)-1)
    do i = 1,n
      cs(i:i) = s(i:i)
    end do
    cs(n+1:n+1) = c_null_char
  end subroutine f_to_c_str

  !> Portable microsecond sleep via C
  subroutine usleep_f(us)
    integer,intent(in) :: us
    call usleep(us)   ! POSIX usleep – available on Linux/macOS
  end subroutine usleep_f

  !> Returns a default Fortran integer for the current process PID
  function get_pid() result(pid)
    integer :: pid
    pid = int(get_pid_shim())
  end function get_pid

end module fmlip_relay_client
