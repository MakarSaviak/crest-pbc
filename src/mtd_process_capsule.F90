!================================================================================!
! Process-isolated MTD capsule for prepared GFN-FF state.
!
! The format is intentionally native-endian and tied to this executable.  It is
! versioned and round-tripped before a child is launched.  The parent and child
! are the same executable, so this preserves every represented real64 bit rather
! than routing prepared MTD state through text or through the original inputs.
!================================================================================!
module mtd_process_capsule
  use iso_fortran_env,only:int8,int32,int64,real64
  use ieee_arithmetic,only:ieee_is_finite
  use crest_parameters,only:wp
  use strucrd,only:coord
  use crest_calculator,only:calcdata,calculation_settings,jobtype
  use calc_type,only:alloc_string
  use dynamics_module,only:mddata,mtdpot,type_mtd,cv_rmsd
  use shake_module,only:shakedata
  implicit none
  private

  integer(int32),parameter,public :: mtd_capsule_version = 5_int32
  integer(int64),parameter :: max_capsule_elements = 10000000_int64
  character(len=32),parameter :: capsule_magic = 'CREST_MTD_PROCESS_CAPSULE_V5    '
  character(len=32),parameter :: capsule_trailer = 'END_CREST_MTD_CAPSULE_V5       '

  type,public :: mtd_capsule_meta
    integer :: worker_index = 0
    integer :: parent_pid = 0
    integer :: inner_threads = 0
    integer :: expected_atoms = 0
    character(len=:),allocatable :: workdir
    character(len=:),allocatable :: topology_file
    character(len=:),allocatable :: result_file
  end type mtd_capsule_meta

  public :: write_mtd_capsule,read_mtd_capsule
  public :: verify_mtd_capsule_roundtrip,compare_binary_files_exact
  public :: supported_process_nci_walls

contains

  subroutine write_mtd_capsule(path,meta,mol,md,calc,iostat_out,message)
    character(len=*),intent(in) :: path
    type(mtd_capsule_meta),intent(in) :: meta
    type(coord),intent(in) :: mol
    type(mddata),intent(in) :: md
    type(calcdata),intent(in) :: calc
    integer,intent(out) :: iostat_out
    character(len=*),intent(out) :: message
    integer :: unit,io

    call validate_capsule_source_state(meta,mol,md,calc,iostat_out,message)
    if (iostat_out /= 0) return
    open(newunit=unit,file=trim(path),access='stream',form='unformatted', &
    & status='replace',action='write',iostat=io,iomsg=message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    write(unit,iostat=io,iomsg=message) capsule_magic
    if (io == 0) write(unit,iostat=io,iomsg=message) mtd_capsule_version
    if (io == 0) write(unit,iostat=io,iomsg=message) int(z'01020304',int32)
    if (io == 0) call write_meta(unit,meta,io,message)
    if (io == 0) call write_coord(unit,mol,io,message)
    if (io == 0) call write_mddata(unit,md,io,message)
    if (io == 0) call write_restricted_calc(unit,mol%nat,calc,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) capsule_trailer
    close(unit,iostat=iostat_out)
    if (io /= 0) iostat_out = io
    if (iostat_out /= 0 .and. len_trim(message) == 0) message = 'capsule write failed'
  end subroutine write_mtd_capsule

  subroutine read_mtd_capsule(path,meta,mol,md,calc,iostat_out,message)
    character(len=*),intent(in) :: path
    type(mtd_capsule_meta),intent(out) :: meta
    type(coord),intent(out) :: mol
    type(mddata),intent(out) :: md
    type(calcdata),intent(out) :: calc
    integer,intent(out) :: iostat_out
    character(len=*),intent(out) :: message
    character(len=32) :: magic,trailer
    integer(int32) :: version,endian
    integer :: unit,io,close_io
    integer(int64) :: stream_position,file_size

    iostat_out = 0
    message = ''
    magic = ''
    trailer = ''
    version = 0_int32
    endian = 0_int32
    stream_position = 0_int64
    file_size = 0_int64
    open(newunit=unit,file=trim(path),access='stream',form='unformatted', &
    & status='old',action='read',iostat=io,iomsg=message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    read(unit,iostat=io,iomsg=message) magic
    if (io == 0 .and. magic /= capsule_magic) then
      io = 1
      message = 'invalid MTD capsule magic'
    end if
    if (io == 0) read(unit,iostat=io,iomsg=message) version
    if (io == 0 .and. version /= mtd_capsule_version) then
      io = 2
      message = 'unsupported MTD capsule version'
    end if
    if (io == 0) read(unit,iostat=io,iomsg=message) endian
    if (io == 0 .and. endian /= int(z'01020304',int32)) then
      io = 3
      message = 'MTD capsule endian mismatch'
    end if
    if (io == 0) call read_meta(unit,meta,io,message)
    if (io == 0) call read_coord(unit,meta%expected_atoms,mol,io,message)
    if (io == 0) call read_mddata(unit,meta%expected_atoms,md,io,message)
    if (io == 0) call read_restricted_calc(unit,meta,calc,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) trailer
    if (io == 0 .and. trailer /= capsule_trailer) then
      io = 4
      message = 'invalid MTD capsule trailer'
    end if
    if (io == 0) then
      inquire(unit=unit,pos=stream_position,iostat=close_io,iomsg=message)
      if (close_io == 0) inquire(file=trim(path),size=file_size,iostat=close_io,iomsg=message)
      if (close_io /= 0) then
        io = close_io
      else if (stream_position /= file_size+1_int64) then
        io = 5
        message = 'trailing data after MTD capsule trailer'
      end if
    end if
    if (io == 0) call validate_capsule_source_state(meta,mol,md,calc,io,message)
    close(unit,iostat=close_io)
    iostat_out = io
    if (iostat_out /= 0 .and. len_trim(message) == 0) message = 'capsule read failed'
  end subroutine read_mtd_capsule

  subroutine verify_mtd_capsule_roundtrip(path,expected_meta,expected_mol, &
  & expected_md,expected_calc,iostat_out,message)
    character(len=*),intent(in) :: path
    type(mtd_capsule_meta),intent(in) :: expected_meta
    type(coord),intent(in) :: expected_mol
    type(mddata),intent(in) :: expected_md
    type(calcdata),intent(in) :: expected_calc
    integer,intent(out) :: iostat_out
    character(len=*),intent(out) :: message
    type(mtd_capsule_meta) :: meta
    type(coord) :: mol
    type(mddata) :: md
    type(calcdata) :: calc
    character(len=:),allocatable :: roundtrip
    integer :: io,compare_io,unit
    logical :: identical

    iostat_out = 0
    message = ''
    roundtrip = trim(path)//'.roundtrip'
    call read_mtd_capsule(path,meta,mol,md,calc,io,message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    call compare_supported_state(expected_meta,expected_mol,expected_md, &
    & expected_calc,meta,mol,md,calc,io,message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    call write_mtd_capsule(roundtrip,meta,mol,md,calc,io,message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    call compare_binary_files_exact(path,roundtrip,identical,compare_io,message)
    if (compare_io /= 0) then
      iostat_out = compare_io
      return
    end if
    if (.not.identical) then
      iostat_out = 6
      message = 'MTD capsule failed byte-exact roundtrip'
      return
    end if
    open(newunit=unit,file=roundtrip,status='old',iostat=io)
    if (io == 0) close(unit,status='delete')
  end subroutine verify_mtd_capsule_roundtrip

  subroutine compare_binary_files_exact(path_a,path_b,identical,iostat_out,message)
    character(len=*),intent(in) :: path_a,path_b
    logical,intent(out) :: identical
    integer,intent(out) :: iostat_out
    character(len=*),intent(out) :: message
    integer(int64) :: size_a,size_b,position,remaining,nread
    integer :: ua,ub,io
    integer(int8),allocatable :: a(:),b(:)
    integer,parameter :: chunk_size = 1048576

    identical = .false.
    iostat_out = 0
    message = ''
    inquire(file=trim(path_a),size=size_a,iostat=io,iomsg=message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    inquire(file=trim(path_b),size=size_b,iostat=io,iomsg=message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    if (size_a /= size_b) return
    open(newunit=ua,file=trim(path_a),access='stream',form='unformatted', &
    & status='old',action='read',iostat=io,iomsg=message)
    if (io /= 0) then
      iostat_out = io
      return
    end if
    open(newunit=ub,file=trim(path_b),access='stream',form='unformatted', &
    & status='old',action='read',iostat=io,iomsg=message)
    if (io /= 0) then
      close(ua)
      iostat_out = io
      return
    end if
    allocate(a(chunk_size),b(chunk_size))
    position = 1_int64
    remaining = size_a
    identical = .true.
    do while (remaining > 0_int64)
      nread = min(remaining,int(chunk_size,int64))
      read(ua,pos=position,iostat=io,iomsg=message) a(1:int(nread))
      if (io /= 0) exit
      read(ub,pos=position,iostat=io,iomsg=message) b(1:int(nread))
      if (io /= 0) exit
      if (any(a(1:int(nread)) /= b(1:int(nread)))) then
        identical = .false.
        exit
      end if
      position = position+nread
      remaining = remaining-nread
    end do
    close(ua)
    close(ub)
    if (io /= 0) then
      identical = .false.
      iostat_out = io
    end if
  end subroutine compare_binary_files_exact

  subroutine validate_capsule_source_state(meta,mol,md,calc,io,message)
    type(mtd_capsule_meta),intent(in) :: meta
    type(coord),intent(in) :: mol
    type(mddata),intent(in) :: md
    type(calcdata),intent(in) :: calc
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    integer :: i,j
    character(len=32) :: worker_text
    io = 0
    message = ''
    if (storage_size(0.0_wp) /= storage_size(0.0_real64)) then
      io = 60
      message = 'capsule requires a real64 CREST working precision'
      return
    end if
    if (.not.all_int32([meta%worker_index,meta%parent_pid,meta%inner_threads, &
    & meta%expected_atoms]) .or. meta%worker_index < 1 .or. meta%parent_pid < 2 .or. &
    & meta%inner_threads < 1 .or. meta%expected_atoms < 1 .or. &
    & .not.nonempty_alloc_string(meta%workdir) .or. &
    & .not.nonempty_alloc_string(meta%topology_file) .or. &
    & .not.nonempty_alloc_string(meta%result_file)) then
      io = 61
      message = 'capsule metadata is incomplete or outside int32 range'
      return
    end if
    if (meta%workdir(1:1) /= '/' .or. meta%topology_file(1:1) /= '/' .or. &
    & meta%result_file(1:1) /= '/') then
      io = 61
      message = 'capsule metadata paths must be absolute'
      return
    end if
    if (trim(meta%topology_file) /= trim(meta%workdir)//'/gfnff_topo' .or. &
    & trim(meta%result_file) /= trim(meta%workdir)//'/result.bin') then
      io = 61
      message = 'capsule topology/result paths escape the private worker contract'
      return
    end if
    write(worker_text,'(i0)') meta%worker_index
    if (.not.all_int32([mol%nat,mol%chrg,mol%uhf,mol%nbd]) .or. mol%nat < 1 .or. &
    & mol%nat /= meta%expected_atoms .or. .not.ieee_is_finite(mol%energy)) then
      io = 62
      message = 'capsule molecule is malformed, nonfinite, or outside int32 range'
      return
    end if
    if (.not.allocated(mol%at) .or. .not.allocated(mol%xyz)) then
      io = 62
      message = 'capsule molecule lacks atoms or coordinates'
      return
    end if
    if (size(mol%at) /= mol%nat .or. size(mol%xyz,1) /= 3 .or. &
    & size(mol%xyz,2) /= mol%nat) then
      io = 62
      message = 'capsule molecule atom/coordinate dimensions are inconsistent'
      return
    end if
    if (.not.all_int_array32(mol%at) .or. .not.all(ieee_is_finite(mol%xyz))) then
      io = 62
      message = 'capsule molecule atom/coordinate payload is invalid'
      return
    end if
    if (allocated(mol%bond) .or. allocated(mol%lat) .or. allocated(mol%qat)) then
      io = 63
      message = 'capsule rejects unused molecule bond, lattice, and charge payloads'
      return
    end if
    if (allocated(mol%bond)) then
      if (.not.all_int_array32_2d(mol%bond)) then
        io = 63
        message = 'capsule bond array is outside int32 range'
        return
      end if
    end if
    if (allocated(mol%lat)) then
      if (.not.all(ieee_is_finite(mol%lat))) then
        io = 64
        message = 'capsule lattice contains a nonfinite value'
        return
      end if
    end if
    if (allocated(mol%qat)) then
      if (.not.all(ieee_is_finite(mol%qat))) then
        io = 65
        message = 'capsule atomic charges contain a nonfinite value'
        return
      end if
    end if
    if (mol%pdb%nat /= 0 .or. mol%pdb%frag /= 0 .or. allocated(mol%pdb%athet) .or. &
    & allocated(mol%pdb%pdbat) .or. allocated(mol%pdb%pdbas) .or. &
    & allocated(mol%pdb%pdbfrag) .or. allocated(mol%pdb%pdbgrp) .or. &
    & allocated(mol%pdb%pdbocc) .or. allocated(mol%pdb%pdbtf)) then
      io = 66
      message = 'capsule does not support PDB metadata'
      return
    end if
    if (.not.all_int32([md%md_index,md%input_structure_id,md%bias_configuration_id, &
    & md%termination_status,md%simtype,md%length_steps,md%sdump,md%dumped, &
    & md%printstep,md%nshake,md%blockl,md%iblock,md%nblock,md%blocknreg, &
    & md%maxblock,md%npot]) .or. .not.all(ieee_is_finite([md%length_ps,md%tstep, &
    & md%dumpstep,md%md_hmass,md%tsoll,md%thermo_damp]))) then
      io = 67
      message = 'capsule MD scalar is nonfinite or outside int32 range'
      return
    end if
    if (md%md_index /= meta%worker_index .or. md%input_structure_id < 0 .or. &
    & md%bias_configuration_id < 0 .or. md%length_steps <= 1 .or. md%tstep <= 0.0_wp .or. &
    & md%sdump <= 0 .or. md%dumped /= 0) then
      io = 67
      message = 'capsule MD index/mapping/step state is invalid'
      return
    end if
    if (.not.nonempty_alloc_string(md%trajectoryfile) .or. &
    & .not.nonempty_alloc_string(md%restartfile)) then
      io = 67
      message = 'capsule fresh dynamics output paths are absent'
      return
    end if
    if (trim(md%trajectoryfile) /= trim(meta%workdir)//'/crest_'// &
    & trim(worker_text)//'.trj' .or. trim(md%restartfile) /= trim(meta%workdir)// &
    & '/crest_'//trim(worker_text)//'.mdrestart') then
      io = 67
      message = 'capsule dynamics outputs escape the private worker contract'
      return
    end if
    if (allocated(md%blockrege) .or. allocated(md%blocke) .or. allocated(md%blockt)) then
      io = 68
      message = 'capsule rejects MD block scratch arrays that dynamics reallocates'
      return
    end if
    if (md%termination_status /= -1 .or. md%iblock /= 0 .or. md%nblock /= 0 .or. &
    & md%blocknreg /= 0) then
      io = 68
      message = 'capsule requires canonical fresh termination and block counters'
      return
    end if
    if (md%shake .or. md%nshake /= 0 .or. md%shk%initialized .or. &
    & md%shk%shake_mode /= 0 .or. md%shk%nusr /= 0 .or. md%shk%ncons /= 0 .or. &
    & allocated(md%shk%conslistu) .or. allocated(md%shk%wbo) .or. &
    & allocated(md%shk%conslist) .or. allocated(md%shk%distcons) .or. &
    & allocated(md%shk%dro) .or. allocated(md%shk%dr) .or. allocated(md%shk%xyzt) .or. &
    & associated(md%shk%freezeptr)) then
      io = 68
      message = 'capsule does not serialize active SHAKE state'
      return
    end if
    if (md%npot /= 1 .or. .not.allocated(md%mtd) .or. .not.allocated(md%cvtype) .or. &
    & allocated(md%active_potentials)) then
      io = 69
      message = 'capsule requires exactly one MTD/CV and no active-potential override'
      return
    end if
    if (size(md%mtd) /= 1 .or. size(md%cvtype) /= 1) then
      io = 69
      message = 'capsule MTD/CV allocation has the wrong size'
      return
    end if
    if (.not.md%requested .or. md%simtype /= type_mtd .or. md%restart .or. md%samerand .or. &
    & md%cvtype(1) /= cv_rmsd .or. md%mtd(1)%mtdtype /= cv_rmsd) then
      io = 69
      message = 'capsule requires fresh stochastic resolved RMSD metadynamics'
      return
    end if
    if (allocated(md%cvtype)) then
      if (.not.all_int_array32(md%cvtype)) then
        io = 69
        message = 'capsule CV type array is outside int32 range'
        return
      end if
    end if
    if (allocated(md%active_potentials)) then
      if (.not.all_int_array32(md%active_potentials)) then
        io = 70
        message = 'capsule active-potential array is outside int32 range'
        return
      end if
    end if
    call validate_shake_state(md%shk,mol%nat,io,message)
    if (io /= 0) return
    if (allocated(md%mtd)) then
      if (size(md%mtd) /= md%npot) then
        io = 71
        message = 'capsule MTD allocation/count mismatch'
        return
      end if
      do i = 1,size(md%mtd)
        if (.not.all_int32([md%mtd(i)%mtdtype,md%mtd(i)%nmax,md%mtd(i)%ncur, &
        & md%mtd(i)%cvdump,md%mtd(i)%cvdumpstep,md%mtd(i)%maxsave, &
        & md%mtd(i)%damptype]) .or. .not.all(ieee_is_finite([md%mtd(i)%kpush, &
        & md%mtd(i)%alpha,md%mtd(i)%com_factor,md%mtd(i)%com_width, &
        & md%mtd(i)%cvdump_fs,md%mtd(i)%ramp,md%mtd(i)%damp]))) then
          io = 72
          message = 'capsule MTD potential is nonfinite or outside int32 range'
          return
        end if
        if (allocated(md%mtd(i)%cv)) then
          if (.not.all(ieee_is_finite(md%mtd(i)%cv))) then
            io = 73; message = 'capsule MTD CV contains a nonfinite value'; return
          end if
        end if
        if (allocated(md%mtd(i)%cvgrd)) then
          if (.not.all(ieee_is_finite(md%mtd(i)%cvgrd))) then
            io = 74; message = 'capsule MTD CV gradient contains a nonfinite value'; return
          end if
        end if
        if (allocated(md%mtd(i)%cvxyz)) then
          if (.not.all(ieee_is_finite(md%mtd(i)%cvxyz))) then
            io = 75; message = 'capsule MTD CV coordinates contain a nonfinite value'; return
          end if
        end if
        if (allocated(md%mtd(i)%damping)) then
          io = 76; message = 'capsule rejects preloaded MTD damping history'; return
        end if
        if (md%mtd(i)%nmax /= 0 .or. md%mtd(i)%ncur /= 0 .or. &
        & md%mtd(i)%cvdump /= 0 .or. md%mtd(i)%cvdumpstep /= 0 .or. &
        & md%mtd(i)%maxsave /= 0) then
          io = 76; message = 'capsule requires canonical fresh MTD counters'; return
        end if
        if (allocated(md%mtd(i)%cv) .or. allocated(md%mtd(i)%cvgrd) .or. &
        & allocated(md%mtd(i)%cvxyz) .or. allocated(md%mtd(i)%biasfile)) then
          io = 76; message = 'capsule rejects preloaded MTD CV/history state'; return
        end if
        if (md%mtd(i)%com_bias .and. &
        & (md%mtd(i)%com_factor <= 0.0_wp .or. md%mtd(i)%com_width <= 0.0_wp)) then
          io = 76; message = 'enabled capsule MTD COM bias has a nonpositive factor or width'; return
        end if
        if (allocated(md%mtd(i)%atinclude)) then
          if (size(md%mtd(i)%atinclude) /= mol%nat) then
            io = 76; message = 'capsule MTD inclusion mask length differs from molecule atom count'; return
          end if
        end if
      end do
    else if (md%npot /= 0) then
      io = 77
      message = 'capsule MTD potentials are missing'
      return
    end if
    if (calc%ncalculations /= 1 .or. .not.allocated(calc%calcs) .or. &
    & .not.supported_process_nci_walls(calc,mol%nat) .or. calc%nscans /= 0 .or. &
    & calc%pr_energies .or. allocated(calc%elog) .or. &
    & allocated(calc%ONIOM) .or. allocated(calc%ONIOMmols) .or. &
    & allocated(calc%ONIOMmap) .or. allocated(calc%ONIOMrevmap)) then
      io = 78
      message = 'capsule supports one automatic NCI wall and no scans or ONIOM state'
      return
    end if
    if (size(calc%calcs) /= 1) then
      io = 78
      message = 'capsule calculator allocation/count is inconsistent'
      return
    end if
    if (calc%id < 0 .or. calc%id > 1) then
      io = 78
      message = 'capsule supports only the one-level weighted or direct calculator selector'
      return
    end if
    if (allocated(calc%etmp) .or. allocated(calc%grdtmp) .or. &
    & allocated(calc%eweight) .or. allocated(calc%weightbackup) .or. &
    & allocated(calc%activebackup) .or. allocated(calc%etmp2) .or. &
    & allocated(calc%grdtmp2) .or. allocated(calc%eweight2) .or. &
    & allocated(calc%grdfix)) then
      io = 78
      message = 'capsule calculator contains process-local runtime scratch state'
      return
    end if
    if (.not.all_int32([calc%id,calc%refine_stage,calc%nfreeze]) .or. &
    & calc%nfreeze < 0) then
      io = 79
      message = 'capsule calculator freeze mask/count is invalid'
      return
    end if
    if (allocated(calc%freezelist)) then
      if (size(calc%freezelist) /= mol%nat .or. count(calc%freezelist) /= calc%nfreeze) then
        io = 79
        message = 'capsule calculator freeze mask/count is inconsistent'
        return
      end if
    else if (calc%nfreeze /= 0) then
      io = 79
      message = 'capsule freeze count is nonzero without a freeze mask'
      return
    end if
    associate(level => calc%calcs(1))
      if (level%id /= jobtype%gfnff .or. .not.level%active .or. &
      & .not.all_int32([level%id,level%refine_lvl,level%chrg,level%uhf]) .or. &
      & .not.ieee_is_finite(level%weight) .or. level%apiclean .or. &
      & .not.level%restart .or. .not.nonempty_alloc_string(level%restartfile)) then
        io = 80
        message = 'capsule requires one persistent active GFN-FF topology restart'
        return
      end if
      if (level%pr .or. level%prappend .or. level%prstdout .or. level%numgrad .or. &
      & .not.level%rdgrad .or. level%rdwbo .or. level%rdqat .or. level%dumpq .or. &
      & level%rddip .or. level%rddipgrad .or. allocated(level%getsasa) .or. &
      & level%getlmocent .or. allocated(level%parametrisation) .or. &
      & allocated(level%refgeo) .or. allocated(level%refcharges) .or. &
      & allocated(level%solvmodel) .or. allocated(level%solvent) .or. &
      & allocated(level%tblite) .or. allocated(level%g0calc) .or. &
      & allocated(level%ff_dat) .or. allocated(level%libpvol) .or. &
      & level%ONIOM_highlowroot /= 0 .or. &
      & level%ONIOM_id /= 0) then
        io = 81
        message = 'capsule calculator contains unsupported output, runtime, path, model, or ONIOM state'
        return
      end if
      if (allocated(level%gff_fragments)) then
        if (size(level%gff_fragments) < 1) then
          io = 82
          message = 'capsule GFN-FF fragment array is allocated but empty'
          return
        end if
        do i = 1,size(level%gff_fragments)
          if (.not.allocated(level%gff_fragments(i)%value)) then
            io = 82
            message = 'capsule GFN-FF fragment string is unallocated'
            return
          end if
          if (len_trim(level%gff_fragments(i)%value) == 0) then
            io = 82
            message = 'capsule GFN-FF fragment string is empty'
            return
          end if
          do j = 1,len_trim(level%gff_fragments(i)%value)
            if (iachar(level%gff_fragments(i)%value(j:j)) < 32 .or. &
            & iachar(level%gff_fragments(i)%value(j:j)) > 126) then
              io = 83
              message = 'capsule GFN-FF fragment string contains non-printable bytes'
              return
            end if
          end do
        end do
      end if
    end associate
  end subroutine validate_capsule_source_state

  subroutine validate_shake_state(shk,nat,io,message)
    type(shakedata),intent(in) :: shk
    integer,intent(in) :: nat
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    io = 0
    message = ''
    if (.not.all_int32([shk%shake_mode,shk%nusr,shk%ncons,shk%maxcyc]) .or. &
    & .not.ieee_is_finite(shk%tolshake) .or. nat < 1 .or. shk%nusr < 0 .or. shk%ncons < 0 .or. &
    & int(shk%nusr,int64) > 8_int64*int(nat,int64) .or. &
    & int(shk%ncons,int64) > 8_int64*int(nat,int64)) then
      io = 84; message = 'capsule SHAKE scalar is invalid'; return
    end if
    if (allocated(shk%conslistu)) then
      if (size(shk%conslistu,1) /= 2 .or. size(shk%conslistu,2) /= shk%nusr .or. &
      & .not.all_int_array32_2d(shk%conslistu)) then
        io = 85; message = 'capsule SHAKE user list is outside int32 range'; return
      end if
    end if
    if (shk%nusr > 0 .and. .not.allocated(shk%conslistu)) then
      io = 85; message = 'capsule SHAKE user constraints are missing'; return
    end if
    if (allocated(shk%conslist)) then
      if (size(shk%conslist,1) /= 2 .or. size(shk%conslist,2) /= shk%ncons .or. &
      & .not.all_int_array32_2d(shk%conslist)) then
        io = 86; message = 'capsule SHAKE list is outside int32 range'; return
      end if
    end if
    if (shk%ncons > 0 .and. .not.allocated(shk%conslist)) then
      io = 86; message = 'capsule SHAKE constraints are missing'; return
    end if
    if (allocated(shk%wbo)) then
      if (size(shk%wbo,1) /= nat .or. size(shk%wbo,2) /= nat .or. &
      & .not.all(ieee_is_finite(shk%wbo))) then
        io = 87; message = 'capsule SHAKE WBO is nonfinite'; return
      end if
    end if
    if (allocated(shk%distcons)) then
      if (size(shk%distcons) /= shk%ncons .or. .not.all(ieee_is_finite(shk%distcons))) then
        io = 88; message = 'capsule SHAKE distances are nonfinite'; return
      end if
    end if
    if (shk%ncons > 0 .and. .not.allocated(shk%distcons)) then
      io = 88; message = 'capsule SHAKE constraint distances are missing'; return
    end if
    if (allocated(shk%dro)) then
      if (size(shk%dro,1) /= 3 .or. size(shk%dro,2) /= shk%ncons) then
        io = 89; message = 'capsule SHAKE dro scratch has the wrong shape'; return
      end if
    end if
    if (shk%ncons > 0 .and. .not.allocated(shk%dro)) then
      io = 89; message = 'capsule SHAKE dro state is missing'; return
    end if
    if (allocated(shk%dr)) then
      if (size(shk%dr,1) /= 4 .or. size(shk%dr,2) /= shk%ncons) then
        io = 90; message = 'capsule SHAKE dr scratch has the wrong shape'; return
      end if
    end if
    if (shk%ncons > 0 .and. .not.allocated(shk%dr)) then
      io = 90; message = 'capsule SHAKE dr state is missing'; return
    end if
    if (allocated(shk%xyzt)) then
      io = 91
      message = 'capsule requires fresh SHAKE coordinate scratch'
      return
    end if
  end subroutine validate_shake_state

  subroutine compare_supported_state(expected_meta,expected_mol,expected_md, &
  & expected_calc,actual_meta,actual_mol,actual_md,actual_calc,io,message)
    type(mtd_capsule_meta),intent(in) :: expected_meta,actual_meta
    type(coord),intent(in) :: expected_mol,actual_mol
    type(mddata),intent(in) :: expected_md,actual_md
    type(calcdata),intent(in) :: expected_calc,actual_calc
    integer,intent(out) :: io
    character(len=*),intent(out) :: message
    io = 0
    message = ''
    if (expected_meta%worker_index /= actual_meta%worker_index .or. &
    & expected_meta%parent_pid /= actual_meta%parent_pid .or. &
    & expected_meta%inner_threads /= actual_meta%inner_threads .or. &
    & expected_meta%expected_atoms /= actual_meta%expected_atoms .or. &
    & .not.same_alloc_string(expected_meta%workdir,actual_meta%workdir) .or. &
    & .not.same_alloc_string(expected_meta%topology_file,actual_meta%topology_file) .or. &
    & .not.same_alloc_string(expected_meta%result_file,actual_meta%result_file)) then
      io = 92
      message = 'decoded capsule metadata differs from parent state'
      return
    end if
    if (.not.same_coord(expected_mol,actual_mol)) then
      io = 93
      message = 'decoded capsule molecule differs from parent state'
      return
    end if
    if (.not.same_mddata(expected_md,actual_md)) then
      io = 94
      message = 'decoded capsule MD state differs from parent state'
      return
    end if
    if (.not.same_restricted_calc(expected_calc,actual_calc)) then
      io = 95
      message = 'decoded capsule calculator state differs from parent state'
      return
    end if
  end subroutine compare_supported_state

  logical function same_coord(a,b) result(same)
    type(coord),intent(in) :: a,b
    same = a%nat == b%nat .and. a%energy == b%energy .and. a%chrg == b%chrg .and. &
    & a%uhf == b%uhf .and. a%nbd == b%nbd .and. same_int_1d(a%at,b%at) .and. &
    & same_real_2d(a%xyz,b%xyz) .and. same_alloc_string(a%comment,b%comment) .and. &
    & same_int_2d(a%bond,b%bond) .and. same_real_2d(a%lat,b%lat) .and. &
    & same_real_1d(a%qat,b%qat)
  end function same_coord

  logical function same_mddata(a,b) result(same)
    type(mddata),intent(in) :: a,b
    integer :: i
    same = a%requested .eqv. b%requested
    same = same .and. a%md_index == b%md_index .and. &
    & a%input_structure_id == b%input_structure_id .and. &
    & a%bias_configuration_id == b%bias_configuration_id .and. &
    & a%termination_status == b%termination_status .and. a%simtype == b%simtype .and. &
    & (a%restart .eqv. b%restart) .and. same_alloc_string(a%restartfile,b%restartfile) .and. &
    & same_alloc_string(a%trajectoryfile,b%trajectoryfile) .and. &
    & a%length_ps == b%length_ps .and. a%length_steps == b%length_steps .and. &
    & a%tstep == b%tstep .and. a%dumpstep == b%dumpstep .and. a%sdump == b%sdump .and. &
    & a%dumped == b%dumped .and. a%printstep == b%printstep .and. &
    & a%md_hmass == b%md_hmass .and. (a%shake .eqv. b%shake) .and. &
    & a%nshake == b%nshake .and. same_shake(a%shk,b%shk) .and. a%tsoll == b%tsoll .and. &
    & (a%thermostat .eqv. b%thermostat) .and. a%thermotype == b%thermotype .and. &
    & a%thermo_damp == b%thermo_damp .and. (a%samerand .eqv. b%samerand) .and. &
    & a%blockl == b%blockl .and. a%iblock == b%iblock .and. a%nblock == b%nblock .and. &
    & a%blocknreg == b%blocknreg .and. a%maxblock == b%maxblock .and. &
    & a%npot == b%npot .and. same_int_1d(a%cvtype,b%cvtype) .and. &
    & same_int_1d(a%active_potentials,b%active_potentials)
    if (.not.same) return
    if (allocated(a%mtd) .neqv. allocated(b%mtd)) then
      same = .false.
      return
    end if
    if (allocated(a%mtd)) then
      if (size(a%mtd) /= size(b%mtd)) then
        same = .false.
        return
      end if
      do i = 1,size(a%mtd)
        if (.not.same_mtdpot(a%mtd(i),b%mtd(i))) then
          same = .false.
          return
        end if
      end do
    end if
  end function same_mddata

  logical function same_mtdpot(a,b) result(same)
    type(mtdpot),intent(in) :: a,b
    same = a%mtdtype == b%mtdtype .and. a%nmax == b%nmax .and. a%ncur == b%ncur .and. &
    & a%kpush == b%kpush .and. a%alpha == b%alpha .and. &
    & (a%com_bias .eqv. b%com_bias) .and. a%com_factor == b%com_factor .and. &
    & a%com_width == b%com_width .and. &
    & (a%com_mass_weighted .eqv. b%com_mass_weighted) .and. &
    & same_real_1d(a%cv,b%cv) .and. same_real_2d(a%cvgrd,b%cvgrd) .and. &
    & a%cvdump == b%cvdump .and. a%cvdump_fs == b%cvdump_fs .and. &
    & a%cvdumpstep == b%cvdumpstep .and. a%maxsave == b%maxsave .and. &
    & same_alloc_string(a%biasfile,b%biasfile) .and. &
    & same_logical_1d(a%atinclude,b%atinclude) .and. &
    & same_real_3d(a%cvxyz,b%cvxyz) .and. a%damptype == b%damptype .and. &
    & a%ramp == b%ramp .and. a%damp == b%damp .and. same_real_1d(a%damping,b%damping)
  end function same_mtdpot

  logical function same_shake(a,b) result(same)
    type(shakedata),intent(in) :: a,b
    same = (a%initialized .eqv. b%initialized) .and. a%shake_mode == b%shake_mode .and. &
    & a%nusr == b%nusr .and. same_int_2d(a%conslistu,b%conslistu) .and. &
    & same_real_2d(a%wbo,b%wbo) .and. a%ncons == b%ncons .and. &
    & same_int_2d(a%conslist,b%conslist) .and. same_real_1d(a%distcons,b%distcons) .and. &
    & a%maxcyc == b%maxcyc .and. a%tolshake == b%tolshake
  end function same_shake

  logical function same_restricted_calc(a,b) result(same)
    type(calcdata),intent(in) :: a,b
    same = a%id == b%id .and. a%refine_stage == b%refine_stage .and. &
    & a%ncalculations == b%ncalculations .and. a%nconstraints == b%nconstraints .and. &
    & a%nfreeze == b%nfreeze .and. same_logical_1d(a%freezelist,b%freezelist) .and. &
    & same_process_nci_walls(a,b)
    if (.not.same .or. a%ncalculations /= 1 .or. .not.allocated(a%calcs) .or. &
    & .not.allocated(b%calcs)) then
      same = .false.
      return
    end if
    if (size(a%calcs) /= 1 .or. size(b%calcs) /= 1) then
      same = .false.
      return
    end if
    same = a%calcs(1)%id == b%calcs(1)%id .and. &
    & a%calcs(1)%refine_lvl == b%calcs(1)%refine_lvl .and. &
    & a%calcs(1)%chrg == b%calcs(1)%chrg .and. a%calcs(1)%uhf == b%calcs(1)%uhf .and. &
    & a%calcs(1)%weight == b%calcs(1)%weight .and. &
    & (a%calcs(1)%active .eqv. b%calcs(1)%active) .and. &
    & (a%calcs(1)%apiclean .eqv. b%calcs(1)%apiclean) .and. &
    & (a%calcs(1)%restart .eqv. b%calcs(1)%restart) .and. &
    & same_alloc_string(a%calcs(1)%parametrisation,b%calcs(1)%parametrisation) .and. &
    & same_alloc_string(a%calcs(1)%refgeo,b%calcs(1)%refgeo) .and. &
    & same_alloc_string(a%calcs(1)%refcharges,b%calcs(1)%refcharges) .and. &
    & same_alloc_string(a%calcs(1)%solvmodel,b%calcs(1)%solvmodel) .and. &
    & same_alloc_string(a%calcs(1)%solvent,b%calcs(1)%solvent) .and. &
    & same_fragment_string_array(a%calcs(1)%gff_fragments,b%calcs(1)%gff_fragments)
    if (allocated(a%calcs(1)%ff_dat) .or. allocated(b%calcs(1)%ff_dat) .or. &
    & allocated(a%etmp) .or. allocated(b%etmp) .or. &
    & allocated(a%grdtmp) .or. allocated(b%grdtmp) .or. &
    & allocated(a%eweight) .or. allocated(b%eweight) .or. &
    & allocated(a%weightbackup) .or. allocated(b%weightbackup) .or. &
    & allocated(a%activebackup) .or. allocated(b%activebackup) .or. &
    & allocated(a%etmp2) .or. allocated(b%etmp2) .or. &
    & allocated(a%grdtmp2) .or. allocated(b%grdtmp2) .or. &
    & allocated(a%eweight2) .or. allocated(b%eweight2) .or. &
    & allocated(a%grdfix) .or. allocated(b%grdfix)) same = .false.
  end function same_restricted_calc

  logical function supported_process_nci_walls(calc,nat) result(supported)
    type(calcdata),intent(in) :: calc
    integer,intent(in) :: nat

    supported = .false.
    if (nat < 1 .or. calc%nconstraints /= 1) return
    if (.not.allocated(calc%cons)) return
    if (size(calc%cons) /= 1) return
    if (.not.calc%cons(1)%is_exact_auto_nci_wall(nat)) return
    if (.not.all_int_array32(calc%cons(1)%atms)) return
    supported = .true.
  end function supported_process_nci_walls

  logical function same_process_nci_walls(a,b) result(same)
    type(calcdata),intent(in) :: a,b
    integer :: i

    if (a%nconstraints /= 1 .or. b%nconstraints /= 1 .or. .not.allocated(a%cons) .or. &
    & .not.allocated(b%cons) .or. size(a%cons) /= 1 .or. size(b%cons) /= 1) then
      same = .false.
      return
    end if
    same = supported_process_nci_walls(a,a%cons(1)%n) .and. &
    & supported_process_nci_walls(b,b%cons(1)%n)
    if (.not.same) return
    do i = 1,a%nconstraints
      same = (a%cons(i)%active .eqv. b%cons(i)%active) .and. &
      & (a%cons(i)%auto_nci_wall .eqv. b%cons(i)%auto_nci_wall) .and. &
      & a%cons(i)%type == b%cons(i)%type .and. a%cons(i)%n == b%cons(i)%n .and. &
      & a%cons(i)%subtype == b%cons(i)%subtype .and. &
      & a%cons(i)%wscal == b%cons(i)%wscal .and. &
      & same_int_1d(a%cons(i)%atms,b%cons(i)%atms) .and. &
      & same_real_1d(a%cons(i)%ref,b%cons(i)%ref) .and. &
      & same_real_1d(a%cons(i)%fc,b%cons(i)%fc)
      if (.not.same) return
    end do
  end function same_process_nci_walls

  logical function same_alloc_string(a,b) result(same)
    character(len=:),allocatable,intent(in) :: a,b
    same = allocated(a) .eqv. allocated(b)
    if (same .and. allocated(a)) same = len(a) == len(b) .and. a == b
  end function same_alloc_string

  logical function same_fragment_string_array(a,b) result(same)
    type(alloc_string),allocatable,intent(in) :: a(:),b(:)
    integer :: i
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = size(a) == size(b)
    if (.not.same) return
    do i = 1,size(a)
      same = allocated(a(i)%value) .eqv. allocated(b(i)%value)
      if (.not.same) return
      if (allocated(a(i)%value)) then
        same = len(a(i)%value) == len(b(i)%value)
        if (.not.same) return
        same = a(i)%value == b(i)%value
        if (.not.same) return
      end if
    end do
  end function same_fragment_string_array

  logical function same_int_1d(a,b) result(same)
    integer,allocatable,intent(in) :: a(:),b(:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = size(a) == size(b)
    if (same) same = all(a == b)
  end function same_int_1d

  logical function same_int_2d(a,b) result(same)
    integer,allocatable,intent(in) :: a(:,:),b(:,:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = all(shape(a) == shape(b))
    if (same) same = all(a == b)
  end function same_int_2d

  logical function same_logical_1d(a,b) result(same)
    logical,allocatable,intent(in) :: a(:),b(:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = size(a) == size(b)
    if (same) same = all(a .eqv. b)
  end function same_logical_1d

  logical function same_real_1d(a,b) result(same)
    real(wp),allocatable,intent(in) :: a(:),b(:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = size(a) == size(b)
    if (same) same = all(a == b)
  end function same_real_1d

  logical function same_real_2d(a,b) result(same)
    real(wp),allocatable,intent(in) :: a(:,:),b(:,:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = all(shape(a) == shape(b))
    if (same) same = all(a == b)
  end function same_real_2d

  logical function same_real_3d(a,b) result(same)
    real(wp),allocatable,intent(in) :: a(:,:,:),b(:,:,:)
    same = allocated(a) .eqv. allocated(b)
    if (.not.same .or. .not.allocated(a)) return
    same = all(shape(a) == shape(b))
    if (same) same = all(a == b)
  end function same_real_3d

  logical function all_int32(values) result(valid)
    integer,intent(in) :: values(:)
    integer(int64),parameter :: lower = -2147483648_int64,upper = 2147483647_int64
    valid = all(int(values,int64) >= lower .and. int(values,int64) <= upper)
  end function all_int32

  logical function all_int_array32(values) result(valid)
    integer,intent(in) :: values(:)
    valid = all_int32(values)
  end function all_int_array32

  logical function all_int_array32_2d(values) result(valid)
    integer,intent(in) :: values(:,:)
    integer(int64),parameter :: lower = -2147483648_int64,upper = 2147483647_int64
    valid = all(int(values,int64) >= lower .and. int(values,int64) <= upper)
  end function all_int_array32_2d

  logical function nonempty_alloc_string(value) result(valid)
    character(len=:),allocatable,intent(in) :: value
    valid = allocated(value)
    if (valid) valid = len_trim(value) > 0
  end function nonempty_alloc_string

  subroutine write_meta(unit,meta,io,message)
    integer,intent(in) :: unit
    type(mtd_capsule_meta),intent(in) :: meta
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    if (io /= 0) return
    write(unit,iostat=io,iomsg=message) int(meta%worker_index,int32), &
    & int(meta%parent_pid,int32),int(meta%inner_threads,int32), &
    & int(meta%expected_atoms,int32)
    if (io == 0) call write_alloc_string(unit,meta%workdir,io,message)
    if (io == 0) call write_alloc_string(unit,meta%topology_file,io,message)
    if (io == 0) call write_alloc_string(unit,meta%result_file,io,message)
  end subroutine write_meta

  subroutine read_meta(unit,meta,io,message)
    integer,intent(in) :: unit
    type(mtd_capsule_meta),intent(out) :: meta
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32) :: worker,parent,threads,natoms
    if (io /= 0) return
    worker = 0_int32
    parent = 0_int32
    threads = 0_int32
    natoms = 0_int32
    read(unit,iostat=io,iomsg=message) worker,parent,threads,natoms
    if (io /= 0) return
    meta%worker_index = int(worker)
    meta%parent_pid = int(parent)
    meta%inner_threads = int(threads)
    meta%expected_atoms = int(natoms)
    call read_alloc_string(unit,meta%workdir,io,message)
    if (io == 0) call read_alloc_string(unit,meta%topology_file,io,message)
    if (io == 0) call read_alloc_string(unit,meta%result_file,io,message)
  end subroutine read_meta

  subroutine write_coord(unit,mol,io,message)
    integer,intent(in) :: unit
    type(coord),intent(in) :: mol
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    if (io /= 0) return
    write(unit,iostat=io,iomsg=message) int(mol%nat,int32),mol%energy, &
    & int(mol%chrg,int32),int(mol%uhf,int32),int(mol%nbd,int32)
    if (io == 0) call write_int_1d(unit,mol%at,io,message)
    if (io == 0) call write_real_2d(unit,mol%xyz,io,message)
    if (io == 0) call write_alloc_string(unit,mol%comment,io,message)
    if (io == 0) call write_int_2d(unit,mol%bond,io,message)
    if (io == 0) call write_real_2d(unit,mol%lat,io,message)
    if (io == 0) call write_real_1d(unit,mol%qat,io,message)
    if (io == 0) then
      ! PDB payload is deliberately unsupported by the process capsule.
      call write_bool(unit,.false.,io,message)
    end if
  end subroutine write_coord

  subroutine read_coord(unit,expected_atoms,mol,io,message)
    integer,intent(in) :: unit
    integer,intent(in) :: expected_atoms
    type(coord),intent(out) :: mol
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32) :: nat,chrg,uhf,nbd
    logical :: pdb_present
    if (io /= 0) return
    nat = 0_int32
    chrg = 0_int32
    uhf = 0_int32
    nbd = 0_int32
    pdb_present = .false.
    read(unit,iostat=io,iomsg=message) nat,mol%energy,chrg,uhf,nbd
    if (io /= 0) return
    if (expected_atoms < 1 .or. nat /= int(expected_atoms,int32)) then
      io = 20
      message = 'process MTD capsule atom count differs from metadata'
      return
    end if
    mol%nat = int(nat)
    mol%chrg = int(chrg)
    mol%uhf = int(uhf)
    mol%nbd = int(nbd)
    call read_int_1d(unit,mol%at,io,message,expected_atoms)
    if (io == 0) call read_real_2d(unit,mol%xyz,io,message,3,expected_atoms)
    if (io == 0) call read_alloc_string(unit,mol%comment,io,message)
    if (io == 0) call read_required_absent(unit,'molecule bond array',io,message)
    if (io == 0) call read_required_absent(unit,'molecule lattice',io,message)
    if (io == 0) call read_required_absent(unit,'molecule charges',io,message)
    if (io == 0) call read_bool(unit,pdb_present,io,message)
    if (io == 0 .and. pdb_present) then
      io = 21
      message = 'PDB payload is unsupported in process MTD capsule'
    end if
    if (io == 0) then
      if (.not.allocated(mol%at) .or. .not.allocated(mol%xyz)) then
        io = 22
        message = 'capsule molecule lacks atoms or coordinates'
      else if (size(mol%at) /= mol%nat .or. size(mol%xyz,1) /= 3 .or. &
      & size(mol%xyz,2) /= mol%nat) then
        io = 23
        message = 'capsule molecule dimensions are inconsistent'
      end if
    end if
  end subroutine read_coord

  subroutine write_mddata(unit,md,io,message)
    integer,intent(in) :: unit
    type(mddata),intent(in) :: md
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer :: i
    if (io /= 0) return
    call write_bool(unit,md%requested,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(md%md_index,int32), &
    & int(md%input_structure_id,int32),int(md%bias_configuration_id,int32), &
    & int(md%termination_status,int32),int(md%simtype,int32)
    if (io == 0) call write_bool(unit,md%restart,io,message)
    if (io == 0) call write_alloc_string(unit,md%restartfile,io,message)
    if (io == 0) call write_alloc_string(unit,md%trajectoryfile,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) md%length_ps, &
    & int(md%length_steps,int32),md%tstep,md%dumpstep,int(md%sdump,int32), &
    & int(md%dumped,int32),int(md%printstep,int32),md%md_hmass
    if (io == 0) call write_bool(unit,md%shake,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(md%nshake,int32)
    if (io == 0) call write_shakedata(unit,md%shk,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) md%tsoll
    if (io == 0) call write_bool(unit,md%thermostat,io,message)
    if (io == 0) call write_fixed_string(unit,md%thermotype,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) md%thermo_damp
    if (io == 0) call write_bool(unit,md%samerand,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(md%blockl,int32), &
    & int(md%iblock,int32),int(md%nblock,int32),int(md%blocknreg,int32), &
    & int(md%maxblock,int32)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(md%npot,int32)
    if (io == 0) then
      call write_bool(unit,allocated(md%mtd),io,message)
      if (io == 0 .and. allocated(md%mtd)) then
        write(unit,iostat=io,iomsg=message) int(size(md%mtd),int32)
        do i = 1,size(md%mtd)
          if (io == 0) call write_mtdpot(unit,md%mtd(i),io,message)
        end do
      end if
    end if
    if (io == 0) call write_int_1d(unit,md%cvtype,io,message)
    if (io == 0) call write_int_1d(unit,md%active_potentials,io,message)
  end subroutine write_mddata

  subroutine read_mddata(unit,expected_atoms,md,io,message)
    integer,intent(in) :: unit
    integer,intent(in) :: expected_atoms
    type(mddata),intent(out) :: md
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32) :: md_index,input_id,bias_id,termination,simtype
    integer(int32) :: length_steps,sdump,dumped,printstep,nshake
    integer(int32) :: blockl,iblock,nblock,blocknreg,maxblock,npot,nmtd
    logical :: mtd_allocated
    character(len=:),allocatable :: thermotype
    if (io /= 0) return
    md_index = 0_int32
    input_id = 0_int32
    bias_id = 0_int32
    termination = 0_int32
    simtype = 0_int32
    length_steps = 0_int32
    sdump = 0_int32
    dumped = 0_int32
    printstep = 0_int32
    nshake = 0_int32
    blockl = 0_int32
    iblock = 0_int32
    nblock = 0_int32
    blocknreg = 0_int32
    maxblock = 0_int32
    npot = 0_int32
    nmtd = 0_int32
    mtd_allocated = .false.
    call read_bool(unit,md%requested,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) md_index,input_id,bias_id, &
    & termination,simtype
    if (io == 0) call read_bool(unit,md%restart,io,message)
    if (io == 0) call read_alloc_string(unit,md%restartfile,io,message)
    if (io == 0) call read_alloc_string(unit,md%trajectoryfile,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) md%length_ps,length_steps, &
    & md%tstep,md%dumpstep,sdump,dumped,printstep,md%md_hmass
    if (io == 0) call read_bool(unit,md%shake,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) nshake
    if (io == 0) call read_shakedata(unit,md%shk,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) md%tsoll
    if (io == 0) call read_bool(unit,md%thermostat,io,message)
    if (io == 0) call read_fixed_string(unit,thermotype,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) md%thermo_damp
    if (io == 0) call read_bool(unit,md%samerand,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) blockl,iblock,nblock, &
    & blocknreg,maxblock
    if (io == 0) read(unit,iostat=io,iomsg=message) npot
    if (io /= 0) return
    if (npot /= 1_int32) then
      io = 30
      message = 'process MTD capsule requires exactly one MTD potential'
      return
    end if
    call read_bool(unit,mtd_allocated,io,message)
    if (io /= 0) return
    if (.not.mtd_allocated) then
      io = 30
      message = 'process MTD capsule MTD potential is absent'
      return
    end if
    read(unit,iostat=io,iomsg=message) nmtd
    if (io /= 0) return
    if (nmtd /= 1_int32) then
      io = 30
      message = 'process MTD capsule MTD allocation count is not one'
      return
    end if
    allocate(md%mtd(1))
    call read_mtdpot(unit,expected_atoms,md%mtd(1),io,message)
    if (io == 0) call read_int_1d(unit,md%cvtype,io,message,1)
    if (io == 0 .and. .not.allocated(md%cvtype)) then
      io = 30
      message = 'process MTD capsule CV type is absent'
    end if
    if (io == 0) call read_required_absent(unit,'active potential override',io,message)
    if (io /= 0) return
    md%md_index = int(md_index)
    md%input_structure_id = int(input_id)
    md%bias_configuration_id = int(bias_id)
    md%termination_status = int(termination)
    md%simtype = int(simtype)
    md%length_steps = int(length_steps)
    md%sdump = int(sdump)
    md%dumped = int(dumped)
    md%printstep = int(printstep)
    md%nshake = int(nshake)
    md%thermotype = thermotype
    md%blockl = int(blockl)
    md%iblock = int(iblock)
    md%nblock = int(nblock)
    md%blocknreg = int(blocknreg)
    md%maxblock = int(maxblock)
    md%npot = int(npot)
    if (allocated(md%mtd)) then
      if (size(md%mtd) /= md%npot) then
        io = 31
        message = 'MTD potential count mismatch in capsule'
      end if
    else if (md%npot /= 0) then
      io = 32
      message = 'missing MTD potentials in capsule'
    end if
  end subroutine read_mddata

  subroutine write_mtdpot(unit,pot,io,message)
    integer,intent(in) :: unit
    type(mtdpot),intent(in) :: pot
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    if (io /= 0) return
    write(unit,iostat=io,iomsg=message) int(pot%mtdtype,int32), &
    & int(pot%nmax,int32),int(pot%ncur,int32),pot%kpush,pot%alpha
    if (io == 0) call write_bool(unit,pot%com_bias,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) pot%com_factor,pot%com_width
    if (io == 0) call write_bool(unit,pot%com_mass_weighted,io,message)
    if (io == 0) call write_real_1d(unit,pot%cv,io,message)
    if (io == 0) call write_real_2d(unit,pot%cvgrd,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(pot%cvdump,int32), &
    & pot%cvdump_fs,int(pot%cvdumpstep,int32),int(pot%maxsave,int32)
    if (io == 0) call write_alloc_string(unit,pot%biasfile,io,message)
    if (io == 0) call write_logical_1d(unit,pot%atinclude,io,message)
    if (io == 0) call write_real_3d(unit,pot%cvxyz,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(pot%damptype,int32), &
    & pot%ramp,pot%damp
    if (io == 0) call write_real_1d(unit,pot%damping,io,message)
  end subroutine write_mtdpot

  subroutine read_mtdpot(unit,expected_atoms,pot,io,message)
    integer,intent(in) :: unit
    integer,intent(in) :: expected_atoms
    type(mtdpot),intent(out) :: pot
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32) :: mtdtype,nmax,ncur,cvdump,cvdumpstep,maxsave,damptype
    if (io /= 0) return
    mtdtype = 0_int32
    nmax = 0_int32
    ncur = 0_int32
    cvdump = 0_int32
    cvdumpstep = 0_int32
    maxsave = 0_int32
    damptype = 0_int32
    read(unit,iostat=io,iomsg=message) mtdtype,nmax,ncur,pot%kpush,pot%alpha
    if (io == 0) call read_bool(unit,pot%com_bias,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) pot%com_factor,pot%com_width
    if (io == 0) call read_bool(unit,pot%com_mass_weighted,io,message)
    if (io == 0) call read_required_absent(unit,'MTD CV history',io,message)
    if (io == 0) call read_required_absent(unit,'MTD CV gradient scratch',io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) cvdump,pot%cvdump_fs, &
    & cvdumpstep,maxsave
    if (io == 0) call read_required_absent(unit,'MTD bias file',io,message)
    if (io == 0) call read_logical_1d(unit,pot%atinclude,io,message,expected_atoms)
    if (io == 0) call read_required_absent(unit,'MTD CV coordinate history',io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) damptype,pot%ramp,pot%damp
    if (io == 0) call read_required_absent(unit,'MTD damping history',io,message)
    if (io /= 0) return
    pot%mtdtype = int(mtdtype)
    pot%nmax = int(nmax)
    pot%ncur = int(ncur)
    pot%cvdump = int(cvdump)
    pot%cvdumpstep = int(cvdumpstep)
    pot%maxsave = int(maxsave)
    pot%damptype = int(damptype)
  end subroutine read_mtdpot

  subroutine write_shakedata(unit,shk,io,message)
    integer,intent(in) :: unit
    type(shakedata),intent(in) :: shk
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    if (io /= 0) return
    call write_bool(unit,shk%initialized,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(shk%shake_mode,int32), &
    & int(shk%nusr,int32)
    if (io == 0) call write_int_2d(unit,shk%conslistu,io,message)
    if (io == 0) call write_real_2d(unit,shk%wbo,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(shk%ncons,int32)
    if (io == 0) call write_int_2d(unit,shk%conslist,io,message)
    if (io == 0) call write_real_1d(unit,shk%distcons,io,message)
    ! init_shake() can allocate dro/dr without defining their contents.  This
    ! exact shake=false target rejects those arrays, and no scratch payload is
    ! serialized or reconstructed.
    if (io == 0) call write_bool(unit,.false.,io,message)
    if (io == 0) call write_bool(unit,.false.,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(shk%maxcyc,int32),shk%tolshake
    ! xyzt is also per-step scratch and must be absent in the fresh target run.
    if (io == 0) call write_bool(unit,.false.,io,message)
    ! A Fortran pointer address cannot cross exec.  The exact shake=false target
    ! requires this semantically-dead pointer to be explicitly null.
    if (io == 0) call write_bool(unit,.false.,io,message)
  end subroutine write_shakedata

  subroutine read_shakedata(unit,shk,io,message)
    integer,intent(in) :: unit
    type(shakedata),intent(out) :: shk
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32) :: shake_mode,nusr,ncons,maxcyc
    logical :: had_freezeptr
    if (io /= 0) return
    shake_mode = 0_int32
    nusr = 0_int32
    ncons = 0_int32
    maxcyc = 0_int32
    had_freezeptr = .false.
    call read_bool(unit,shk%initialized,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) shake_mode,nusr
    if (io /= 0) return
    if (shk%initialized .or. shake_mode /= 0_int32 .or. nusr /= 0_int32) then
      io = 33
      message = 'process MTD capsule rejects initialized or active SHAKE state'
      return
    end if
    call read_required_absent(unit,'SHAKE user constraints',io,message)
    if (io == 0) call read_required_absent(unit,'SHAKE WBO state',io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) ncons
    if (io /= 0) return
    if (ncons /= 0_int32) then
      io = 34
      message = 'process MTD capsule rejects SHAKE constraints'
      return
    end if
    call read_required_absent(unit,'SHAKE constraints',io,message)
    if (io == 0) call read_required_absent(unit,'SHAKE constraint distances',io,message)
    if (io == 0) call read_required_absent(unit,'SHAKE dro scratch',io,message)
    if (io == 0) call read_required_absent(unit,'SHAKE dr scratch',io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) maxcyc,shk%tolshake
    if (io == 0) call read_required_absent(unit,'SHAKE coordinate scratch',io,message)
    if (io == 0) call read_bool(unit,had_freezeptr,io,message)
    if (io /= 0) return
    if (had_freezeptr) then
      io = 35
      message = 'capsule cannot transfer a process-local SHAKE freeze pointer'
      return
    end if
    shk%shake_mode = int(shake_mode)
    shk%nusr = int(nusr)
    shk%ncons = int(ncons)
    shk%maxcyc = int(maxcyc)
    ! The exact shake=false target keeps this process-local pointer null.
    nullify(shk%freezeptr)
  end subroutine read_shakedata

  subroutine write_restricted_calc(unit,nat,calc,io,message)
    integer,intent(in) :: unit
    integer,intent(in) :: nat
    type(calcdata),intent(in) :: calc
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    if (io /= 0) return
    if (calc%ncalculations /= 1 .or. .not.allocated(calc%calcs)) then
      io = 40
      message = 'process MTD capsule supports exactly one calculator level'
      return
    end if
    associate(level => calc%calcs(1))
    if (level%id /= jobtype%gfnff .or. .not.level%active) then
      io = 41
      message = 'process MTD capsule supports one active GFN-FF level only'
      return
    end if
    if (.not.supported_process_nci_walls(calc,nat) .or. calc%nscans /= 0 .or. &
    & allocated(calc%ONIOM)) then
      io = 42
      message = 'process MTD supports one automatic all-atom NCI wall; scans and ONIOM are unsupported'
      return
    end if
    if (level%numgrad .or. level%rdwbo .or. level%rdqat .or. level%dumpq .or. &
    & level%rddip .or. level%rddipgrad .or. &
    & allocated(level%getsasa) .or. level%getlmocent) then
      io = 43
      message = 'requested calculator properties are unsupported in process MTD prototype'
      return
    end if
    if (level%apiclean .or. .not.level%restart .or. .not.allocated(level%restartfile)) then
      io = 44
      message = 'process MTD requires persistent GFN-FF state from an explicit topology restart'
      return
    end if
    if (allocated(level%parametrisation) .or. allocated(level%refgeo) .or. &
    & allocated(level%refcharges) .or. allocated(level%solvmodel) .or. &
    & allocated(level%solvent)) then
      io = 45
      message = 'custom GFN-FF parameter/reference/solvent inputs are unsupported in the process MTD capsule'
      return
    end if
    write(unit,iostat=io,iomsg=message) int(calc%id,int32), &
    & int(calc%refine_stage,int32),int(calc%nfreeze,int32)
    if (io == 0) call write_logical_1d(unit,calc%freezelist,io,message)
    if (io == 0) call write_process_nci_walls(unit,nat,calc,io,message)
    if (io == 0) write(unit,iostat=io,iomsg=message) int(level%id,int32), &
    & int(level%refine_lvl,int32),int(level%chrg,int32),int(level%uhf,int32), &
    & level%weight
    if (io == 0) call write_bool(unit,level%active,io,message)
    if (io == 0) call write_bool(unit,level%apiclean,io,message)
    if (io == 0) call write_bool(unit,level%restart,io,message)
    if (io == 0) call write_alloc_string(unit,level%parametrisation,io,message)
    if (io == 0) call write_alloc_string(unit,level%refgeo,io,message)
    if (io == 0) call write_alloc_string(unit,level%refcharges,io,message)
    if (io == 0) call write_alloc_string(unit,level%solvmodel,io,message)
    if (io == 0) call write_alloc_string(unit,level%solvent,io,message)
    if (io == 0) call write_fragment_string_array(unit,level%gff_fragments,io,message)
    end associate
  end subroutine write_restricted_calc

  subroutine read_restricted_calc(unit,meta,calc,io,message)
    integer,intent(in) :: unit
    type(mtd_capsule_meta),intent(in) :: meta
    type(calcdata),intent(out) :: calc
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    type(calculation_settings) :: level
    integer(int32) :: calc_id,refine_stage,nfreeze,level_id,refine_lvl,chrg,uhf
    if (io /= 0) return
    calc_id = 0_int32
    refine_stage = 0_int32
    nfreeze = 0_int32
    level_id = 0_int32
    refine_lvl = 0_int32
    chrg = 0_int32
    uhf = 0_int32
    read(unit,iostat=io,iomsg=message) calc_id,refine_stage,nfreeze
    if (io == 0) call read_logical_1d(unit,calc%freezelist,io,message,meta%expected_atoms)
    if (io == 0) call read_process_nci_walls(unit,meta%expected_atoms,calc,io,message)
    if (io == 0) read(unit,iostat=io,iomsg=message) level_id,refine_lvl,chrg,uhf,level%weight
    if (io == 0) call read_bool(unit,level%active,io,message)
    if (io == 0) call read_bool(unit,level%apiclean,io,message)
    if (io == 0) call read_bool(unit,level%restart,io,message)
    if (io == 0) call read_alloc_string(unit,level%parametrisation,io,message)
    if (io == 0) call read_alloc_string(unit,level%refgeo,io,message)
    if (io == 0) call read_alloc_string(unit,level%refcharges,io,message)
    if (io == 0) call read_alloc_string(unit,level%solvmodel,io,message)
    if (io == 0) call read_alloc_string(unit,level%solvent,io,message)
    if (io == 0) call read_fragment_string_array(unit,level%gff_fragments,io,message)
    if (io /= 0) return
    if (level_id /= jobtype%gfnff .or. .not.level%active) then
      io = 46
      message = 'capsule calculator is not active GFN-FF'
      return
    end if
    if (.not.allocated(meta%topology_file) .or. .not.allocated(meta%workdir)) then
      io = 47
      message = 'capsule lacks private topology/work directory'
      return
    end if
    calc%id = int(calc_id)
    calc%refine_stage = int(refine_stage)
    calc%nfreeze = int(nfreeze)
    if (calc%nfreeze > 0) then
      if (.not.allocated(calc%freezelist)) then
        io = 48
        message = 'capsule freeze count has no mask'
        return
      end if
    end if
    level%id = int(level_id)
    level%refine_lvl = int(refine_lvl)
    level%chrg = int(chrg)
    level%uhf = int(uhf)
    calc%ncalculations = 1
    allocate(calc%calcs(1))
    calc%calcs(1)%id = level%id
    calc%calcs(1)%refine_lvl = level%refine_lvl
    calc%calcs(1)%chrg = level%chrg
    calc%calcs(1)%uhf = level%uhf
    calc%calcs(1)%weight = level%weight
    calc%calcs(1)%active = level%active
    calc%calcs(1)%apiclean = .false.
    calc%calcs(1)%restart = .true.
    calc%calcs(1)%restartfile = meta%topology_file
    calc%calcs(1)%calcspace = trim(meta%workdir)//'/CALC'
    calc%calcs(1)%pr = .false.
    calc%calcs(1)%prstdout = .false.
    calc%calcs(1)%prappend = .false.
    calc%calcs(1)%rdgrad = .true.
    calc%calcs(1)%rdwbo = .false.
    calc%calcs(1)%rdqat = .false.
    calc%calcs(1)%dumpq = .false.
    calc%calcs(1)%rddip = .false.
    calc%calcs(1)%rddipgrad = .false.
    calc%calcs(1)%numgrad = .false.
    if (allocated(level%gff_fragments)) &
    & call move_alloc(level%gff_fragments,calc%calcs(1)%gff_fragments)
    calc%pr_energies = .false.
  end subroutine read_restricted_calc

  subroutine write_process_nci_walls(unit,nat,calc,io,message)
    integer,intent(in) :: unit
    integer,intent(in) :: nat
    type(calcdata),intent(in) :: calc
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer :: i

    if (io /= 0) return
    if (.not.supported_process_nci_walls(calc,nat)) then
      io = 49
      message = 'calculator does not contain the supported automatic NCI wall'
      return
    end if
    write(unit,iostat=io,iomsg=message) int(calc%nconstraints,int32)
    do i = 1,calc%nconstraints
      associate(wall => calc%cons(i))
        if (io == 0) write(unit,iostat=io,iomsg=message) int(wall%type,int32), &
        & int(wall%n,int32),int(wall%subtype,int32),wall%wscal
        if (io == 0) call write_bool(unit,wall%active,io,message)
        if (io == 0) call write_bool(unit,wall%auto_nci_wall,io,message)
        if (io == 0) call write_int_1d(unit,wall%atms,io,message)
        if (io == 0) call write_real_1d(unit,wall%ref,io,message)
        if (io == 0) call write_real_1d(unit,wall%fc,io,message)
      end associate
    end do
  end subroutine write_process_nci_walls

  subroutine read_process_nci_walls(unit,expected_atoms,calc,io,message)
    integer,intent(in) :: unit
    integer,intent(in) :: expected_atoms
    type(calcdata),intent(inout) :: calc
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32) :: nconstraints,wall_type,wall_n,wall_subtype
    integer :: i

    if (io /= 0) return
    nconstraints = 0_int32
    wall_type = 0_int32
    wall_n = 0_int32
    wall_subtype = 0_int32
    read(unit,iostat=io,iomsg=message) nconstraints
    if (io /= 0) return
    if (expected_atoms < 1 .or. nconstraints /= 1_int32) then
      io = 49
      message = 'capsule calculator NCI wall count is unsupported'
      return
    end if
    calc%nconstraints = int(nconstraints)
    allocate(calc%cons(int(nconstraints)))
    do i = 1,calc%nconstraints
      nullify(calc%cons(i)%freezeptr)
      calc%cons(i)%frozenatms = .false.
      read(unit,iostat=io,iomsg=message) wall_type,wall_n,wall_subtype,calc%cons(i)%wscal
      if (io /= 0) return
      if (wall_n /= int(expected_atoms,int32)) then
        io = 49
        message = 'capsule calculator NCI wall atom count differs from metadata'
        return
      end if
      calc%cons(i)%type = int(wall_type)
      calc%cons(i)%n = int(wall_n)
      calc%cons(i)%subtype = int(wall_subtype)
      call read_bool(unit,calc%cons(i)%active,io,message)
      if (io == 0) call read_bool(unit,calc%cons(i)%auto_nci_wall,io,message)
      if (io == 0) call read_int_1d(unit,calc%cons(i)%atms,io,message,expected_atoms)
      if (io == 0) call read_real_1d(unit,calc%cons(i)%ref,io,message,3)
      if (io == 0) call read_real_1d(unit,calc%cons(i)%fc,io,message,2)
      if (io /= 0) return
    end do
    if (.not.supported_process_nci_walls(calc,expected_atoms)) then
      io = 49
      message = 'decoded capsule NCI wall violates the supported scientific contract'
    end if
  end subroutine read_process_nci_walls

  subroutine write_bool(unit,value,io,message)
    integer,intent(in) :: unit
    logical,intent(in) :: value
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int8) :: encoded
    encoded = merge(1_int8,0_int8,value)
    write(unit,iostat=io,iomsg=message) encoded
  end subroutine write_bool

  subroutine read_bool(unit,value,io,message)
    integer,intent(in) :: unit
    logical,intent(out) :: value
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int8) :: encoded
    if (io /= 0) return
    value = .false.
    encoded = 0_int8
    read(unit,iostat=io,iomsg=message) encoded
    if (io /= 0) return
    if (encoded /= 0_int8 .and. encoded /= 1_int8) then
      io = 50
      message = 'invalid logical encoding in capsule'
      return
    end if
    value = encoded == 1_int8
  end subroutine read_bool

  subroutine read_required_absent(unit,label,io,message)
    integer,intent(in) :: unit
    character(len=*),intent(in) :: label
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    logical :: encoded_present
    if (io /= 0) return
    encoded_present = .false.
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0) return
    if (encoded_present) then
      io = 53
      message = 'unsupported capsule payload is present: '//trim(label)
    end if
  end subroutine read_required_absent

  subroutine write_alloc_string(unit,value,io,message)
    integer,intent(in) :: unit
    character(len=:),allocatable,intent(in) :: value
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    call write_bool(unit,allocated(value),io,message)
    if (io == 0 .and. allocated(value)) call write_fixed_string(unit,value,io,message)
  end subroutine write_alloc_string

  subroutine read_alloc_string(unit,value,io,message)
    integer,intent(in) :: unit
    character(len=:),allocatable,intent(out) :: value
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    logical :: encoded_present
    encoded_present = .false.
    call read_bool(unit,encoded_present,io,message)
    if (io == 0 .and. encoded_present) call read_fixed_string(unit,value,io,message)
  end subroutine read_alloc_string

  subroutine write_fixed_string(unit,value,io,message)
    integer,intent(in) :: unit
    character(len=*),intent(in) :: value
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    write(unit,iostat=io,iomsg=message) int(len(value),int64)
    if (io == 0 .and. len(value) > 0) write(unit,iostat=io,iomsg=message) value
  end subroutine write_fixed_string

  subroutine read_fixed_string(unit,value,io,message)
    integer,intent(in) :: unit
    character(len=:),allocatable,intent(out) :: value
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int64) :: n
    n = 0_int64
    read(unit,iostat=io,iomsg=message) n
    if (io /= 0) return
    if (n < 0_int64 .or. n > 1048576_int64) then
      io = 51
      message = 'invalid string length in capsule'
      return
    end if
    allocate(character(len=int(n)) :: value)
    if (n > 0) read(unit,iostat=io,iomsg=message) value
  end subroutine read_fixed_string

  subroutine write_character_array(unit,value,io,message)
    integer,intent(in) :: unit
    character(len=:),allocatable,intent(in) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer :: i
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    write(unit,iostat=io,iomsg=message) int(size(value),int64),int(len(value),int64)
    if (io == 0 .and. len(value) > 0) then
      do i = 1,size(value)
        write(unit,iostat=io,iomsg=message) value(i)
        if (io /= 0) exit
      end do
    end if
  end subroutine write_character_array

  subroutine read_character_array(unit,value,io,message,expected_n,max_length)
    integer,intent(in) :: unit
    character(len=:),allocatable,intent(out) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer,intent(in),optional :: expected_n,max_length
    logical :: encoded_present
    integer(int64) :: n,l
    integer :: i
    encoded_present = .false.
    n = 0_int64
    l = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n,l
    if (io /= 0) return
    if (.not.valid_shape2(n,l) .or. l > 1048576_int64) then
      io = 52
      message = 'invalid character array dimensions in capsule'
      return
    end if
    if (present(expected_n)) then
      if (n /= int(expected_n,int64)) then
        io = 52; message = 'character array count violates capsule field contract'; return
      end if
    end if
    if (present(max_length)) then
      if (l > int(max_length,int64)) then
        io = 52; message = 'character array width violates capsule field contract'; return
      end if
    end if
    allocate(character(len=int(l)) :: value(int(n)))
    if (l > 0) then
      do i = 1,int(n)
        read(unit,iostat=io,iomsg=message) value(i)
        if (io /= 0) exit
      end do
    end if
  end subroutine read_character_array

  subroutine write_fragment_string_array(unit,value,io,message)
    integer,intent(in) :: unit
    type(alloc_string),allocatable,intent(in) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer :: i
    if (io /= 0) return
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    if (.not.valid_elements(int(size(value),int64))) then
      io = 96
      message = 'cannot encode an oversized GFN-FF fragment array'
      return
    end if
    write(unit,iostat=io,iomsg=message) int(size(value),int64)
    do i = 1,size(value)
      if (.not.allocated(value(i)%value)) then
        io = 96
        message = 'cannot encode an unallocated GFN-FF fragment string'
        return
      end if
      call write_fixed_string(unit,value(i)%value,io,message)
      if (io /= 0) return
    end do
  end subroutine write_fragment_string_array

  subroutine read_fragment_string_array(unit,value,io,message)
    integer,intent(in) :: unit
    type(alloc_string),allocatable,intent(out) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    logical :: encoded_present
    integer(int64) :: n
    integer :: i
    if (io /= 0) return
    encoded_present = .false.
    n = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n
    if (io /= 0) return
    if (.not.valid_elements(n)) then
      io = 96
      message = 'invalid GFN-FF fragment array length in capsule'
      return
    end if
    allocate(value(int(n)))
    do i = 1,int(n)
      call read_fixed_string(unit,value(i)%value,io,message)
      if (io /= 0) return
    end do
  end subroutine read_fragment_string_array

  subroutine write_int_1d(unit,value,io,message)
    integer,intent(in) :: unit
    integer,allocatable,intent(in) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32),allocatable :: encoded(:)
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    write(unit,iostat=io,iomsg=message) int(size(value),int64)
    if (io == 0 .and. size(value) > 0) then
      encoded = int(value,int32)
      write(unit,iostat=io,iomsg=message) encoded
    end if
  end subroutine write_int_1d

  subroutine read_int_1d(unit,value,io,message,expected_n)
    integer,intent(in) :: unit
    integer,allocatable,intent(out) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer,intent(in),optional :: expected_n
    logical :: encoded_present
    integer(int64) :: n
    integer(int32),allocatable :: encoded(:)
    encoded_present = .false.
    n = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n
    if (io /= 0) return
    if (.not.valid_elements(n)) then
      io = 60
      message = 'invalid integer vector length in capsule'
      return
    end if
    if (present(expected_n)) then
      if (n /= int(expected_n,int64)) then
        io = 60
        message = 'integer vector length violates capsule field contract'
        return
      end if
    end if
    allocate(encoded(int(n)),value(int(n)))
    if (n > 0) read(unit,iostat=io,iomsg=message) encoded
    if (io == 0) value = int(encoded)
  end subroutine read_int_1d

  subroutine write_int_2d(unit,value,io,message)
    integer,intent(in) :: unit
    integer,allocatable,intent(in) :: value(:,:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int32),allocatable :: encoded(:,:)
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    write(unit,iostat=io,iomsg=message) int(size(value,1),int64),int(size(value,2),int64)
    if (io == 0 .and. size(value) > 0) then
      encoded = int(value,int32)
      write(unit,iostat=io,iomsg=message) encoded
    end if
  end subroutine write_int_2d

  subroutine read_int_2d(unit,value,io,message,expected_n1,expected_n2)
    integer,intent(in) :: unit
    integer,allocatable,intent(out) :: value(:,:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer,intent(in),optional :: expected_n1,expected_n2
    logical :: encoded_present
    integer(int64) :: n1,n2
    integer(int32),allocatable :: encoded(:,:)
    encoded_present = .false.
    n1 = 0_int64
    n2 = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n1,n2
    if (io /= 0) return
    if (.not.valid_shape2(n1,n2)) then
      io = 61
      message = 'invalid integer matrix shape in capsule'
      return
    end if
    if (present(expected_n1)) then
      if (n1 /= int(expected_n1,int64)) then
        io = 61; message = 'integer matrix first dimension violates field contract'; return
      end if
    end if
    if (present(expected_n2)) then
      if (n2 /= int(expected_n2,int64)) then
        io = 61; message = 'integer matrix second dimension violates field contract'; return
      end if
    end if
    allocate(encoded(int(n1),int(n2)),value(int(n1),int(n2)))
    if (n1*n2 > 0) read(unit,iostat=io,iomsg=message) encoded
    if (io == 0) value = int(encoded)
  end subroutine read_int_2d

  subroutine write_logical_1d(unit,value,io,message)
    integer,intent(in) :: unit
    logical,allocatable,intent(in) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer(int8),allocatable :: encoded(:)
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    write(unit,iostat=io,iomsg=message) int(size(value),int64)
    if (io == 0 .and. size(value) > 0) then
      allocate(encoded(size(value)),source=0_int8)
      where(value) encoded = 1_int8
      write(unit,iostat=io,iomsg=message) encoded
    end if
  end subroutine write_logical_1d

  subroutine read_logical_1d(unit,value,io,message,expected_n)
    integer,intent(in) :: unit
    logical,allocatable,intent(out) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer,intent(in),optional :: expected_n
    logical :: encoded_present
    integer(int64) :: n
    integer(int8),allocatable :: encoded(:)
    encoded_present = .false.
    n = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n
    if (io /= 0) return
    if (.not.valid_elements(n)) then
      io = 62
      message = 'invalid logical vector length in capsule'
      return
    end if
    if (present(expected_n)) then
      if (n /= int(expected_n,int64)) then
        io = 62; message = 'logical vector length violates capsule field contract'; return
      end if
    end if
    allocate(encoded(int(n)),value(int(n)))
    if (n > 0) read(unit,iostat=io,iomsg=message) encoded
    if (io /= 0) return
    if (any(encoded /= 0_int8 .and. encoded /= 1_int8)) then
      io = 63
      message = 'invalid logical vector encoding in capsule'
      return
    end if
    value = encoded == 1_int8
  end subroutine read_logical_1d

  subroutine write_real_1d(unit,value,io,message)
    integer,intent(in) :: unit
    real(wp),allocatable,intent(in) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    write(unit,iostat=io,iomsg=message) int(size(value),int64)
    if (io == 0 .and. size(value) > 0) write(unit,iostat=io,iomsg=message) value
  end subroutine write_real_1d

  subroutine read_real_1d(unit,value,io,message,expected_n)
    integer,intent(in) :: unit
    real(wp),allocatable,intent(out) :: value(:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer,intent(in),optional :: expected_n
    logical :: encoded_present
    integer(int64) :: n
    encoded_present = .false.
    n = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n
    if (io /= 0) return
    if (.not.valid_elements(n)) then
      io = 64
      message = 'invalid real vector length in capsule'
      return
    end if
    if (present(expected_n)) then
      if (n /= int(expected_n,int64)) then
        io = 64; message = 'real vector length violates capsule field contract'; return
      end if
    end if
    allocate(value(int(n)))
    if (n > 0) read(unit,iostat=io,iomsg=message) value
  end subroutine read_real_1d

  subroutine write_real_2d(unit,value,io,message)
    integer,intent(in) :: unit
    real(wp),allocatable,intent(in) :: value(:,:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    write(unit,iostat=io,iomsg=message) int(size(value,1),int64),int(size(value,2),int64)
    if (io == 0 .and. size(value) > 0) write(unit,iostat=io,iomsg=message) value
  end subroutine write_real_2d

  subroutine read_real_2d(unit,value,io,message,expected_n1,expected_n2)
    integer,intent(in) :: unit
    real(wp),allocatable,intent(out) :: value(:,:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer,intent(in),optional :: expected_n1,expected_n2
    logical :: encoded_present
    integer(int64) :: n1,n2
    encoded_present = .false.
    n1 = 0_int64
    n2 = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n1,n2
    if (io /= 0) return
    if (.not.valid_shape2(n1,n2)) then
      io = 65
      message = 'invalid real matrix shape in capsule'
      return
    end if
    if (present(expected_n1)) then
      if (n1 /= int(expected_n1,int64)) then
        io = 65; message = 'real matrix first dimension violates field contract'; return
      end if
    end if
    if (present(expected_n2)) then
      if (n2 /= int(expected_n2,int64)) then
        io = 65; message = 'real matrix second dimension violates field contract'; return
      end if
    end if
    allocate(value(int(n1),int(n2)))
    if (n1*n2 > 0) read(unit,iostat=io,iomsg=message) value
  end subroutine read_real_2d

  subroutine write_real_3d(unit,value,io,message)
    integer,intent(in) :: unit
    real(wp),allocatable,intent(in) :: value(:,:,:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    call write_bool(unit,allocated(value),io,message)
    if (io /= 0 .or. .not.allocated(value)) return
    write(unit,iostat=io,iomsg=message) int(size(value,1),int64), &
    & int(size(value,2),int64),int(size(value,3),int64)
    if (io == 0 .and. size(value) > 0) write(unit,iostat=io,iomsg=message) value
  end subroutine write_real_3d

  subroutine read_real_3d(unit,value,io,message,expected_n1,expected_n2,expected_n3)
    integer,intent(in) :: unit
    real(wp),allocatable,intent(out) :: value(:,:,:)
    integer,intent(inout) :: io
    character(len=*),intent(inout) :: message
    integer,intent(in),optional :: expected_n1,expected_n2,expected_n3
    logical :: encoded_present
    integer(int64) :: n1,n2,n3
    encoded_present = .false.
    n1 = 0_int64
    n2 = 0_int64
    n3 = 0_int64
    call read_bool(unit,encoded_present,io,message)
    if (io /= 0 .or. .not.encoded_present) return
    read(unit,iostat=io,iomsg=message) n1,n2,n3
    if (io /= 0) return
    if (.not.valid_shape3(n1,n2,n3)) then
      io = 66
      message = 'invalid real tensor shape in capsule'
      return
    end if
    if (present(expected_n1)) then
      if (n1 /= int(expected_n1,int64)) then
        io = 66; message = 'real tensor first dimension violates field contract'; return
      end if
    end if
    if (present(expected_n2)) then
      if (n2 /= int(expected_n2,int64)) then
        io = 66; message = 'real tensor second dimension violates field contract'; return
      end if
    end if
    if (present(expected_n3)) then
      if (n3 /= int(expected_n3,int64)) then
        io = 66; message = 'real tensor third dimension violates field contract'; return
      end if
    end if
    allocate(value(int(n1),int(n2),int(n3)))
    if (n1*n2*n3 > 0) read(unit,iostat=io,iomsg=message) value
  end subroutine read_real_3d

  logical function valid_elements(n)
    integer(int64),intent(in) :: n
    valid_elements = n >= 0_int64 .and. n <= max_capsule_elements
  end function valid_elements

  logical function valid_shape2(n1,n2)
    integer(int64),intent(in) :: n1,n2
    valid_shape2 = valid_elements(n1) .and. valid_elements(n2)
    if (valid_shape2 .and. n1 > 0) valid_shape2 = n2 <= max_capsule_elements/n1
  end function valid_shape2

  logical function valid_shape3(n1,n2,n3)
    integer(int64),intent(in) :: n1,n2,n3
    integer(int64) :: product
    valid_shape3 = valid_shape2(n1,n2) .and. valid_elements(n3)
    if (.not.valid_shape3) return
    product = n1*n2
    if (product > 0) valid_shape3 = n3 <= max_capsule_elements/product
  end function valid_shape3

end module mtd_process_capsule
