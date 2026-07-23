!> example/example_usage.f90
!> ─────────────────────────────────────────────────────────────────────────────
!> Demonstrates single-instance and parallel multi-instance usage of
!> the fmlip_relay client module against an installed fmlip-relay-server.
!>
!> Prerequisites:
!>   pip install -e .[mace]          (from the fmlip_relay project root)
!>
!> Build:
!>   cd example && make
!>
!> Run (dummy backend – no model required):
!>   ./build/run_example dummy
!>
!> Run (Lennard-Jones backend):
!>   ./build/run_example lj
!>
!> Run (MACE-MP foundation model – downloads automatically):
!>   ./build/run_example mace_mp
!>   ./build/run_example mace_mp medium-mpa-0
!>
!> Run (MACE-OFF foundation model – organic elements only):
!>   ./build/run_example mace_off
!>
!> Run (custom-trained MACE model):
!>   ./build/run_example mace /path/to/model.model
!>
!> Geometry: FCC carbon supercell (64 atoms, a = 3.57 Å).
!> Carbon (Z=6) is supported by all MACE model families.
!> ─────────────────────────────────────────────────────────────────────────────

program example_usage
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use fmlip_relay_client
  use omp_lib,only:omp_get_thread_num
  implicit none

  ! ── runtime arguments ────────────────────────────────────────────────────────
  character(len=64)  :: backend_name
  character(len=256) :: model_arg             ! model path (mace) or size tag (mace_mp/mace_off)
  character(len=512) :: server_cmd

  ! ── parameters ────────────────────────────────────────────────────────────────
  integer,parameter :: TIMEOUT_SEC = 120
  integer,parameter :: BASE_PORT = 54320
  integer,parameter :: NINSTANCES = 2
  integer,parameter :: NATOMS = 64
  !> FCC lattice parameter for carbon (Å).
  !> Nearest-neighbour distance = a/sqrt(2) ≈ 2.52 Å – reasonable for a demo.
  real(wp),parameter :: ALAT_DEFAULT = 3.57_wp

  ! ── local variables ──────────────────────────────────────────────────────────
  real(wp) :: positions(3,NATOMS),cell(3,3)
  integer :: pbc(3),atomic_numbers(NATOMS)
  real(wp) :: energy,forces(3,NATOMS),stress(3,3)
  real(wp) :: lj_epsilon,lj_sigma,lj_cutoff
  integer :: ierr,i,iid,pid

  ! ── parse arguments ──────────────────────────────────────────────────────────
  ! Usage:
  !   run_example dummy
  !   run_example lj [epsilon] [sigma] [cutoff]
  !   run_example mace_mp [model_size]         e.g. medium-mpa-0
  !   run_example mace_off [model_size]        e.g. small | medium | large
  !   run_example mace <model_path>
  if (command_argument_count() < 1) then
    write (stdout,*) "Usage: run_example <backend> [options]"
    write (stdout,*) "  dummy"
    write (stdout,*) "  lj  [epsilon_eV]  [sigma_ang]  [cutoff_ang]"
    write (stdout,*) "  mace_mp  [model_size]"
    write (stdout,*) "  mace_off [model_size]"
    write (stdout,*) "  mace  model_path"
    stop 1
  end if
  call get_command_argument(1,backend_name)

  ! Defaults
  lj_epsilon = 0.0104_wp
  lj_sigma = 3.40_wp
  lj_cutoff = 8.50_wp
  model_arg = ""

  select case (trim(backend_name))
  case ("lj")
    if (command_argument_count() >= 2) then
      call get_command_argument(2,model_arg); read (model_arg,*) lj_epsilon
    end if
    if (command_argument_count() >= 3) then
      call get_command_argument(3,model_arg); read (model_arg,*) lj_sigma
    end if
    if (command_argument_count() >= 4) then
      call get_command_argument(4,model_arg); read (model_arg,*) lj_cutoff
    end if
    model_arg = ""
  case ("mace","mace_mp","mace_off")
    ! optional: model path (mace) or size tag (mace_mp / mace_off)
    if (command_argument_count() >= 2) call get_command_argument(2,model_arg)
  end select

  ! ════════════════════════════════════════════════════════════════════════════
  ! Demo 1: single instance
  ! ════════════════════════════════════════════════════════════════════════════
  write (stdout,*) "=== Demo 1: single instance ==="

  pid = get_pid()
  write (stdout,*) "Process PID:",pid

  call make_fcc_positions(NATOMS,ALAT_DEFAULT,positions,cell)
  pbc = [1,1,1]
  atomic_numbers = 6    ! carbon – supported by all MACE model families

  call build_server_cmd(backend_name,model_arg,BASE_PORT+1, &  !> BASE_PORT could be replaced with pid
                        lj_epsilon,lj_sigma,lj_cutoff,server_cmd)
  call mlip_init(1,BASE_PORT+1,trim(server_cmd),TIMEOUT_SEC,ierr)
  if (ierr /= MLIP_OK) stop "mlip_init failed"

  call mlip_ping(1,ierr)
  if (ierr == MLIP_OK) write (stdout,*) " PING OK"

  ! charge=0, spin=1 (neutral, ground-state singlet). Both are optional and
  ! default to exactly these values, so they may be omitted via keyword args:
  !   call mlip_compute(1,NATOMS,...,pbc,1,energy=energy,forces=forces, &
  !                     stress=stress,ierr=ierr)
  call mlip_compute(1,NATOMS,atomic_numbers,positions,cell,pbc,1,0,1, &
                    energy,forces,stress,ierr)
  if (ierr /= MLIP_OK) stop "mlip_compute failed"

  write (stdout,'(A,F16.8,A)') "  Energy          = ",energy," eV"
  write (stdout,'(A,3F12.6)') "  Forces[atom 1]  = ",forces(:,1)
  write (stdout,'(A,3F12.6)') "  Stress diagonal = ", &
    stress(1,1),stress(2,2),stress(3,3)

  call mlip_finalize(1,ierr)
  write (stdout,*) "Instance 1 finalized."

  ! ════════════════════════════════════════════════════════════════════════════
  ! Demo 2: multiple parallel instances (OpenMP)
  ! ════════════════════════════════════════════════════════════════════════════
  write (stdout,*)
  write (stdout,*) "=== Demo 2: parallel instances (OpenMP) ==="

  pid = get_pid()
  write (stdout,*) "Process PID:",pid

  !> NOTE we are using the PID to assign communication channels in this example
  do i = 1,NINSTANCES

    call build_server_cmd(backend_name,model_arg,pid+i, &
                          lj_epsilon,lj_sigma,lj_cutoff,server_cmd)
    call mlip_init(i,pid+i,trim(server_cmd),TIMEOUT_SEC,ierr)

    if (ierr /= MLIP_OK) then
      write (stdout,'(A,I0)') "mlip_init failed for instance ",i; stop
    end if
  end do

  !$omp parallel do num_threads(NINSTANCES) &
  !$omp& private(i, iid, positions, cell, pbc, atomic_numbers, &
  !$omp&         energy, forces, stress, ierr) &
  !$omp& schedule(static,1)
  do i = 1,NINSTANCES
    iid = omp_get_thread_num()+1

    call make_fcc_positions(NATOMS,ALAT_DEFAULT+0.01_wp*i,positions,cell)
    pbc = [1,1,1]
    atomic_numbers = 6

    call mlip_compute(iid,NATOMS,atomic_numbers,positions,cell,pbc,0,0,1, &
                      energy,forces,stress,ierr)
    if (ierr == MLIP_OK) then
      !$omp critical
      write (stdout,'(A,I0,A,F14.6,A)') "  Thread ",iid,"  energy = ",energy," eV"
      !$omp end critical
    end if
  end do
  !$omp end parallel do

  call mlip_finalize_all(ierr)
  write (stdout,*) "All instances finalized."

contains

  !> Assemble the full shell command for fmlip-relay-server.
  subroutine build_server_cmd(backend,model,port,eps,sig,rc,cmd)
    character(len=*),intent(in)  :: backend,model
    integer,intent(in)  :: port
    real(wp),intent(in)  :: eps,sig,rc
    character(len=*),intent(out) :: cmd

    select case (trim(backend))
    case ("mace")
      write (cmd,'("fmlip-relay-server --port ",I0, &
                  &" --backend mace --model ",A)') port,trim(model)
    case ("mace_mp")
      if (len_trim(model) > 0) then
        write (cmd,'("fmlip-relay-server --port ",I0, &
                    &" --backend mace_mp --mace-model ",A)') port,trim(model)
      else
        write (cmd,'("fmlip-relay-server --port ",I0, &
                    &" --backend mace_mp")') port
      end if
    case ("mace_off")
      if (len_trim(model) > 0) then
        write (cmd,'("fmlip-relay-server --port ",I0, &
                    &" --backend mace_off --mace-model ",A)') port,trim(model)
      else
        write (cmd,'("fmlip-relay-server --port ",I0, &
                    &" --backend mace_off")') port
      end if
    case ("lj")
      write (cmd,'("fmlip-relay-server --port ",I0, &
                  &" --backend lj --lj-epsilon ",F10.6, &
                  &" --lj-sigma ",F10.6," --lj-cutoff ",F10.6)') &
            port,eps,sig,rc
    case default   ! dummy or anything unrecognised
      write (cmd,'("fmlip-relay-server --port ",I0, &
                  &" --backend dummy")') port
    end select
  end subroutine build_server_cmd

  !> Build a simple FCC supercell (placeholder geometry).
  subroutine make_fcc_positions(natoms,alat,pos,cellout)
    integer,intent(in)  :: natoms
    real(wp),intent(in)  :: alat
    real(wp),intent(out) :: pos(3,natoms),cellout(3,3)
    real(wp) :: dx(3,4)
    integer :: ia,ix,iy,iz,ncell,j

    dx(:,1) = [0.0_wp,0.0_wp,0.0_wp]
    dx(:,2) = [0.5_wp,0.5_wp,0.0_wp]
    dx(:,3) = [0.5_wp,0.0_wp,0.5_wp]
    dx(:,4) = [0.0_wp,0.5_wp,0.5_wp]

    ncell = nint((natoms/4.0_wp)**(1.0_wp/3.0_wp))
    cellout = 0.0_wp
    cellout(1,1) = ncell*alat
    cellout(2,2) = ncell*alat
    cellout(3,3) = ncell*alat

    ia = 0
    outer: do iz = 0,ncell-1
      do iy = 0,ncell-1
        do ix = 0,ncell-1
          do j = 1,4
            ia = ia+1
            if (ia > natoms) exit outer
            pos(:,ia) = (dx(:,j)+[dble(ix),dble(iy),dble(iz)])*alat
          end do
        end do
      end do
    end do outer
  end subroutine make_fcc_positions

end program example_usage
