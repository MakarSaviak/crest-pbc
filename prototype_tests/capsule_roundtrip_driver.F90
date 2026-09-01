program capsule_roundtrip_driver
  use iso_fortran_env,only:int8,int64
  use crest_parameters,only:wp
  use strucrd,only:coord
  use crest_calculator,only:calcdata,constraint,jobtype
  use constraints,only:calc_constraint
  use dynamics_module,only:mddata,type_mtd,cv_rmsd
  use mtd_process_capsule,only:mtd_capsule_meta,write_mtd_capsule, &
  & read_mtd_capsule,verify_mtd_capsule_roundtrip
  implicit none
  type(mtd_capsule_meta) :: meta,decoded_meta
  type(coord) :: mol,decoded_mol
  type(mddata) :: md,decoded_md
  type(calcdata) :: calc,decoded_calc
  type(constraint) :: wall_copy
  character(len=4096) :: capsule,truncated,endian_bad,trailing,v4_header,rejected
  character(len=1024) :: message
  integer :: i,j,io
  real(wp) :: parent_wall_energy,decoded_wall_energy
  real(wp),allocatable :: parent_wall_gradient(:,:),decoded_wall_gradient(:,:)

  if (command_argument_count() /= 1) error stop 'usage: capsule_roundtrip_driver CAPSULE'
  call get_command_argument(1,capsule)

  mol%nat = 938
  allocate(mol%at(mol%nat),mol%xyz(3,mol%nat))
  mol%at = 6
  do i = 1,mol%nat
    mol%xyz(:,i) = [real(i,wp),real(2*i,wp),real(-i,wp)]*1.0e-4_wp
  end do
  mol%energy = -123.456789012345_wp
  mol%comment = 'capsule roundtrip sentinel'

  md%requested = .true.
  md%shake = .false.
  nullify(md%shk%freezeptr)
  md%md_index = 1
  md%input_structure_id = 1
  md%bias_configuration_id = 1
  md%termination_status = -1
  md%simtype = type_mtd
  md%trajectoryfile = '/tmp/capsule-worker/crest_1.trj'
  md%restartfile = '/tmp/capsule-worker/crest_1.mdrestart'
  md%length_ps = 1.5_wp
  md%length_steps = 3000
  md%tstep = 0.5_wp
  md%dumpstep = 50.0_wp
  md%sdump = 100
  md%printstep = 200
  md%md_hmass = 4.0_wp
  md%tsoll = 298.15_wp
  md%thermostat = .true.
  md%thermotype = 'berendsen'
  md%thermo_damp = 500.0_wp
  md%samerand = .false.
  md%blockl = 100
  md%maxblock = 30
  md%npot = 1
  allocate(md%mtd(1),md%cvtype(1))
  md%cvtype(1) = cv_rmsd
  md%mtd(1)%mtdtype = cv_rmsd
  md%mtd(1)%kpush = 0.03_wp
  md%mtd(1)%alpha = 0.5_wp
  md%mtd(1)%com_bias = .true.
  md%mtd(1)%com_factor = 0.001_wp
  md%mtd(1)%com_width = 0.05_wp
  md%mtd(1)%com_mass_weighted = .true.
  md%mtd(1)%cvdump_fs = 50.0_wp
  allocate(md%mtd(1)%atinclude(mol%nat),source=.false.)
  md%mtd(1)%atinclude(809:938) = .true.

  calc%id = 0
  calc%nfreeze = 808
  allocate(calc%freezelist(mol%nat),source=.false.)
  calc%freezelist(1:808) = .true.
  calc%nconstraints = 1
  allocate(calc%cons(1))
  do j = 1,calc%nconstraints
    calc%cons(j)%active = .true.
    calc%cons(j)%auto_nci_wall = .true.
    calc%cons(j)%type = 5
    calc%cons(j)%n = mol%nat
    calc%cons(j)%subtype = 1
    calc%cons(j)%wscal = 1.0_wp
    allocate(calc%cons(j)%atms(mol%nat),calc%cons(j)%ref(3),calc%cons(j)%fc(2))
    do i = 1,mol%nat
      calc%cons(j)%atms(i) = i
    end do
    calc%cons(j)%fc = [298.15_wp,50.0_wp]
    call calc%cons(j)%addfreeze(calc%freezelist)
  end do
  calc%cons(1)%ref = [0.20_wp,0.21_wp,0.22_wp]
  calc%ncalculations = 1
  allocate(calc%calcs(1))
  calc%calcs(1)%id = jobtype%gfnff
  calc%calcs(1)%active = .true.
  calc%calcs(1)%weight = 1.0_wp
  calc%calcs(1)%chrg = 0
  calc%calcs(1)%uhf = 0
  calc%calcs(1)%apiclean = .false.
  calc%calcs(1)%restart = .true.
  calc%calcs(1)%restartfile = '/tmp/capsule-worker/gfnff_topo'
  allocate(calc%calcs(1)%gff_fragments(2))
  allocate(character(len=16) :: calc%calcs(1)%gff_fragments(1)%value)
  allocate(character(len=16) :: calc%calcs(1)%gff_fragments(2)%value)
  calc%calcs(1)%gff_fragments(1)%value(:) = '1-808'
  calc%calcs(1)%gff_fragments(2)%value(:) = '809-938'

  meta%worker_index = 1
  meta%parent_pid = 12345
  meta%inner_threads = 8
  meta%expected_atoms = 938
  meta%workdir = '/tmp/capsule-worker'
  meta%topology_file = '/tmp/capsule-worker/gfnff_topo'
  meta%result_file = '/tmp/capsule-worker/result.bin'

  call write_mtd_capsule(trim(capsule),meta,mol,md,calc,io,message)
  if (io /= 0) error stop 'capsule write failed: '//trim(message)
  call verify_mtd_capsule_roundtrip(trim(capsule),meta,mol,md,calc,io,message)
  if (io /= 0) error stop 'capsule roundtrip failed: '//trim(message)
  call read_mtd_capsule(trim(capsule),decoded_meta,decoded_mol,decoded_md, &
  & decoded_calc,io,message)
  if (io /= 0) error stop 'capsule decode failed: '//trim(message)
  if (decoded_meta%worker_index /= 1 .or. decoded_meta%inner_threads /= 8) &
  & error stop 'decoded metadata mismatch'
  if (decoded_mol%nat /= mol%nat .or. any(decoded_mol%at /= mol%at) .or. &
  & any(decoded_mol%xyz /= mol%xyz)) error stop 'decoded molecule mismatch'
  if (decoded_md%md_index /= md%md_index .or. decoded_md%samerand .or. &
  & decoded_md%mtd(1)%com_factor /= md%mtd(1)%com_factor .or. &
  & any(decoded_md%mtd(1)%atinclude .neqv. md%mtd(1)%atinclude)) &
  & error stop 'decoded MTD state mismatch'
  if (decoded_calc%nfreeze /= calc%nfreeze .or. &
  & any(decoded_calc%freezelist .neqv. calc%freezelist) .or. &
  & len(decoded_calc%calcs(1)%gff_fragments(1)%value) /= 16 .or. &
  & len(decoded_calc%calcs(1)%gff_fragments(2)%value) /= 16 .or. &
  & decoded_calc%calcs(1)%gff_fragments(1)%value /= &
  & calc%calcs(1)%gff_fragments(1)%value .or. &
  & decoded_calc%calcs(1)%gff_fragments(2)%value /= &
  & calc%calcs(1)%gff_fragments(2)%value) &
  & error stop 'decoded calculator state mismatch'
  if (decoded_calc%nconstraints /= 1) error stop 'decoded NCI wall count mismatch'
  if (.not.allocated(decoded_calc%cons)) error stop 'decoded NCI wall is absent'
  if (size(decoded_calc%cons) /= 1) error stop 'decoded NCI wall size mismatch'
  allocate(parent_wall_gradient(3,mol%nat),decoded_wall_gradient(3,mol%nat))
  do j = 1,calc%nconstraints
    if (decoded_calc%cons(j)%type /= calc%cons(j)%type .or. &
    & decoded_calc%cons(j)%n /= calc%cons(j)%n .or. &
    & decoded_calc%cons(j)%subtype /= calc%cons(j)%subtype .or. &
    & decoded_calc%cons(j)%wscal /= calc%cons(j)%wscal .or. &
    & (decoded_calc%cons(j)%active .neqv. calc%cons(j)%active) .or. &
    & (decoded_calc%cons(j)%auto_nci_wall .neqv. calc%cons(j)%auto_nci_wall) .or. &
    & any(decoded_calc%cons(j)%atms /= calc%cons(j)%atms) .or. &
    & any(decoded_calc%cons(j)%ref /= calc%cons(j)%ref) .or. &
    & any(decoded_calc%cons(j)%fc /= calc%cons(j)%fc)) &
    & error stop 'decoded NCI wall field or order mismatch'
    if (associated(decoded_calc%cons(j)%freezeptr) .or. &
    & decoded_calc%cons(j)%frozenatms) &
    & error stop 'decoded NCI wall retained a process-local freeze association'
    call calc_constraint(mol%nat,mol%xyz,calc%cons(j),parent_wall_energy, &
    & parent_wall_gradient)
    call decoded_calc%cons(j)%addfreeze(decoded_calc%freezelist)
    call calc_constraint(decoded_mol%nat,decoded_mol%xyz,decoded_calc%cons(j), &
    & decoded_wall_energy,decoded_wall_gradient)
    if (decoded_wall_energy /= parent_wall_energy .or. &
    & any(decoded_wall_gradient /= parent_wall_gradient)) &
    & error stop 'decoded NCI wall energy or full gradient is not bitwise exact'
    if (any(parent_wall_gradient(:,1:808) /= 0.0_wp) .or. &
    & .not.any(parent_wall_gradient(:,809:938) /= 0.0_wp)) &
    & error stop 'NCI wall frozen-host skip was not preserved'
  end do
  deallocate(parent_wall_gradient,decoded_wall_gradient)
  truncated = trim(capsule)//'.truncated'
  endian_bad = trim(capsule)//'.bad-endian'
  trailing = trim(capsule)//'.trailing-byte'
  v4_header = trim(capsule)//'.v4-header'
  rejected = trim(capsule)//'.must-not-write'
  call make_capsule_variant(trim(capsule),trim(truncated),.false.)
  call read_mtd_capsule(trim(truncated),decoded_meta,decoded_mol,decoded_md, &
  & decoded_calc,io,message)
  if (io == 0) error stop 'truncated capsule was accepted'
  call make_capsule_variant(trim(capsule),trim(endian_bad),.true.)
  call read_mtd_capsule(trim(endian_bad),decoded_meta,decoded_mol,decoded_md, &
  & decoded_calc,io,message)
  if (io == 0) error stop 'wrong-endian capsule was accepted'
  call make_v4_header_variant(trim(capsule),trim(v4_header))
  call read_mtd_capsule(trim(v4_header),decoded_meta,decoded_mol,decoded_md, &
  & decoded_calc,io,message)
  if (io == 0) error stop 'V4 capsule header was accepted by the V5 reader'
  call make_trailing_variant(trim(capsule),trim(trailing))
  call read_mtd_capsule(trim(trailing),decoded_meta,decoded_mol,decoded_md, &
  & decoded_calc,io,message)
  if (io == 0) error stop 'capsule with one trailing byte was accepted'
  decoded_meta = meta
  decoded_meta%worker_index = 2
  call verify_mtd_capsule_roundtrip(trim(capsule),decoded_meta,mol,md,calc,io,message)
  if (io == 0) error stop 'source-vs-decoded metadata mismatch was accepted'
  allocate(md%blocke(1),source=0.0_wp)
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'preallocated dynamics scratch was accepted'
  deallocate(md%blocke)
  md%shake = .true.
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'active SHAKE state was accepted'
  md%shake = .false.
  md%sdump = 0
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'zero trajectory dump interval was accepted'
  md%sdump = 100
  md%mtd(1)%cvdump = 1
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'nonfresh MTD counter was accepted'
  md%mtd(1)%cvdump = 0
  deallocate(calc%calcs(1)%gff_fragments)
  allocate(calc%calcs(1)%gff_fragments(1))
  allocate(character(len=16) :: calc%calcs(1)%gff_fragments(1)%value)
  calc%calcs(1)%gff_fragments(1)%value(:) = '1-808'
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io /= 0) error stop 'single-fragment calculator state was rejected'
  calc%calcs(1)%gff_fragments(1)%value(:) = ''
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'empty GFN-FF fragment text was accepted'
  deallocate(calc%calcs(1)%gff_fragments)
  allocate(calc%calcs(1)%gff_fragments(2))
  allocate(character(len=16) :: calc%calcs(1)%gff_fragments(1)%value)
  allocate(character(len=16) :: calc%calcs(1)%gff_fragments(2)%value)
  calc%calcs(1)%gff_fragments(1)%value(:) = '1-808'
  calc%calcs(1)%gff_fragments(2)%value(:) = '809-938'
  deallocate(calc%calcs(1)%gff_fragments(2)%value)
  allocate(character(len=7) :: calc%calcs(1)%gff_fragments(2)%value)
  calc%calcs(1)%gff_fragments(2)%value(:) = '809-938'
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io /= 0) error stop 'unequal GFN-FF fragment string widths were rejected'
  deallocate(calc%calcs(1)%gff_fragments(2)%value)
  allocate(character(len=16) :: calc%calcs(1)%gff_fragments(2)%value)
  calc%calcs(1)%gff_fragments(2)%value(:) = '809-938'
  calc%calcs(1)%refgeo = 'relative-refgeo.xyz'
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'unsupported relative refgeo was accepted'
  deallocate(calc%calcs(1)%refgeo)
  allocate(calc%calcs(1)%ff_dat)
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'allocated GFN-FF runtime state was silently dropped'
  deallocate(calc%calcs(1)%ff_dat)
  allocate(calc%eweight(1),source=calc%calcs(1)%weight)
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'allocated calculator weight cache was silently dropped'
  deallocate(calc%eweight)
  do i = -2,11
    if (i == 5) cycle
    calc%cons(1)%type = i
    call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
    if (io == 0) error stop 'a non-log-Fermi constraint type was accepted'
  end do
  calc%cons(1)%type = 5
  calc%cons(1)%atms(mol%nat) = mol%nat-1
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'a noncanonical NCI wall atom list was accepted'
  calc%cons(1)%atms(mol%nat) = mol%nat
  calc%cons(1)%auto_nci_wall = .false.
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'an unmarked NCI wall was accepted'
  calc%cons(1)%auto_nci_wall = .true.
  wall_copy = calc%cons(1)
  deallocate(calc%cons)
  calc%nconstraints = 2
  allocate(calc%cons(2))
  calc%cons(1) = wall_copy
  calc%cons(2) = wall_copy
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'a duplicate two-wall calculator was accepted'
  deallocate(calc%cons)
  calc%nconstraints = 0
  call write_mtd_capsule(trim(rejected),meta,mol,md,calc,io,message)
  if (io == 0) error stop 'a calculator without the required NCI wall was accepted'
  calc%nconstraints = 1
  allocate(calc%cons(1))
  calc%cons(1) = wall_copy
  call run_generic_case('generic-A',17,1)
  call run_generic_case('generic-B',47,2)
  write(*,'(a)') 'CAPSULE_NCI_WALLS_EXACT_GATE_PASS'
  write(*,'(a)') 'CAPSULE_GENERIC_CASES_PASS'
  write(*,'(a)') 'CAPSULE_NEGATIVE_GATES_PASS'
  write(*,'(a)') 'CAPSULE_ROUNDTRIP_PASS'
contains
  subroutine run_generic_case(label,nat,case_id)
    character(len=*),intent(in) :: label
    integer,intent(in) :: nat,case_id
    type(mtd_capsule_meta) :: case_meta,case_decoded_meta
    type(coord) :: case_mol,case_decoded_mol
    type(mddata) :: case_md,case_decoded_md
    type(calcdata) :: case_calc,case_decoded_calc
    character(len=4096) :: case_capsule
    character(len=1024) :: case_message
    integer :: case_i,case_io

    case_mol%nat = nat
    allocate(case_mol%at(nat),case_mol%xyz(3,nat))
    case_mol%at = 6
    do case_i = 1,nat
      case_mol%xyz(:,case_i) = [real(case_i,wp),real(-2*case_i,wp), &
      & real(case_i+3,wp)]*1.0e-3_wp
    end do
    case_mol%energy = -real(nat,wp)
    case_mol%comment = trim(label)//' generic capsule case'

    case_md%requested = .true.
    case_md%shake = .false.
    nullify(case_md%shk%freezeptr)
    case_md%md_index = 1
    case_md%input_structure_id = case_id
    case_md%bias_configuration_id = case_id+10
    case_md%termination_status = -1
    case_md%simtype = type_mtd
    case_md%trajectoryfile = '/tmp/capsule-worker/crest_1.trj'
    case_md%restartfile = '/tmp/capsule-worker/crest_1.mdrestart'
    case_md%length_ps = 0.2_wp
    case_md%length_steps = 400
    case_md%tstep = 0.5_wp
    case_md%dumpstep = 20.0_wp
    case_md%sdump = 20
    case_md%printstep = 40
    case_md%md_hmass = 4.0_wp
    case_md%tsoll = 298.15_wp
    case_md%thermostat = .true.
    case_md%thermotype = 'berendsen'
    case_md%thermo_damp = 100.0_wp
    case_md%samerand = .false.
    case_md%blockl = 20
    case_md%maxblock = 20
    case_md%npot = 1
    allocate(case_md%mtd(1),case_md%cvtype(1))
    case_md%cvtype(1) = cv_rmsd
    case_md%mtd(1)%mtdtype = cv_rmsd
    case_md%mtd(1)%kpush = 0.025_wp
    case_md%mtd(1)%alpha = 0.4_wp
    case_md%mtd(1)%cvdump_fs = 20.0_wp
    allocate(case_md%mtd(1)%atinclude(nat),source=.false.)
    if (case_id == 1) then
      case_md%mtd(1)%atinclude([2,7,11,16]) = .true.
      case_md%mtd(1)%com_bias = .false.
      case_md%mtd(1)%com_factor = 0.0_wp
      case_md%mtd(1)%com_width = 0.0_wp
      case_md%mtd(1)%com_mass_weighted = .false.
    else
      case_md%mtd(1)%atinclude([3,8,12,28,40]) = .true.
      case_md%mtd(1)%com_bias = .true.
      case_md%mtd(1)%com_factor = 0.002_wp
      case_md%mtd(1)%com_width = 0.07_wp
      case_md%mtd(1)%com_mass_weighted = .false.
    end if

    case_calc%id = 0
    if (case_id == 1) then
      case_calc%nfreeze = 0
    else
      allocate(case_calc%freezelist(nat),source=.false.)
      case_calc%freezelist([2,5,9,17,31,46]) = .true.
      case_calc%nfreeze = count(case_calc%freezelist)
    end if
    case_calc%nconstraints = 1
    allocate(case_calc%cons(1))
    case_calc%cons(1)%active = .true.
    case_calc%cons(1)%auto_nci_wall = .true.
    case_calc%cons(1)%type = 5
    case_calc%cons(1)%n = nat
    case_calc%cons(1)%subtype = 1
    case_calc%cons(1)%wscal = 1.0_wp
    allocate(case_calc%cons(1)%atms(nat),case_calc%cons(1)%ref(3),case_calc%cons(1)%fc(2))
    do case_i = 1,nat
      case_calc%cons(1)%atms(case_i) = case_i
    end do
    case_calc%cons(1)%ref = [0.2_wp,0.3_wp,0.4_wp]
    case_calc%cons(1)%fc = [298.15_wp,50.0_wp]
    if (allocated(case_calc%freezelist)) call case_calc%cons(1)%addfreeze(case_calc%freezelist)
    case_calc%ncalculations = 1
    allocate(case_calc%calcs(1))
    case_calc%calcs(1)%id = jobtype%gfnff
    case_calc%calcs(1)%active = .true.
    case_calc%calcs(1)%weight = 1.0_wp
    case_calc%calcs(1)%chrg = 0
    case_calc%calcs(1)%uhf = 0
    case_calc%calcs(1)%apiclean = .false.
    case_calc%calcs(1)%restart = .true.
    case_calc%calcs(1)%restartfile = '/tmp/capsule-worker/gfnff_topo'
    if (case_id == 1) then
      allocate(case_calc%calcs(1)%gff_fragments(2))
      allocate(character(len=3) :: case_calc%calcs(1)%gff_fragments(1)%value)
      allocate(character(len=7) :: case_calc%calcs(1)%gff_fragments(2)%value)
      case_calc%calcs(1)%gff_fragments(1)%value = '1-8'
      case_calc%calcs(1)%gff_fragments(2)%value = '9,11-17'
    else
      allocate(case_calc%calcs(1)%gff_fragments(3))
      allocate(character(len=4) :: case_calc%calcs(1)%gff_fragments(1)%value)
      allocate(character(len=8) :: case_calc%calcs(1)%gff_fragments(2)%value)
      allocate(character(len=21) :: case_calc%calcs(1)%gff_fragments(3)%value)
      case_calc%calcs(1)%gff_fragments(1)%value = '1-12'
      case_calc%calcs(1)%gff_fragments(2)%value = '13,15-22'
      case_calc%calcs(1)%gff_fragments(3)%value = '14,23-37,40-42,44-47'
    end if

    case_meta%worker_index = 1
    case_meta%parent_pid = 12345
    case_meta%inner_threads = 2
    case_meta%expected_atoms = nat
    case_meta%workdir = '/tmp/capsule-worker'
    case_meta%topology_file = '/tmp/capsule-worker/gfnff_topo'
    case_meta%result_file = '/tmp/capsule-worker/result.bin'
    case_capsule = '/tmp/'//trim(label)//'.capsule'
    call write_mtd_capsule(trim(case_capsule),case_meta,case_mol,case_md,case_calc,case_io,case_message)
    if (case_io /= 0) error stop trim(label)//' capsule write failed: '//trim(case_message)
    call verify_mtd_capsule_roundtrip(trim(case_capsule),case_meta,case_mol,case_md,case_calc, &
    & case_io,case_message)
    if (case_io /= 0) error stop trim(label)//' capsule roundtrip failed: '//trim(case_message)
    call read_mtd_capsule(trim(case_capsule),case_decoded_meta,case_decoded_mol,case_decoded_md, &
    & case_decoded_calc,case_io,case_message)
    if (case_io /= 0) error stop trim(label)//' capsule decode failed: '//trim(case_message)
    if (case_decoded_meta%expected_atoms /= nat .or. case_decoded_mol%nat /= nat .or. &
    & any(case_decoded_mol%at /= case_mol%at) .or. any(case_decoded_mol%xyz /= case_mol%xyz)) then
      error stop trim(label)//' molecule roundtrip mismatch'
    end if
    if (case_decoded_calc%nfreeze /= case_calc%nfreeze .or. &
    & (allocated(case_decoded_calc%freezelist) .neqv. allocated(case_calc%freezelist))) then
      error stop trim(label)//' freeze-state roundtrip mismatch'
    end if
    if (allocated(case_calc%freezelist)) then
      if (any(case_decoded_calc%freezelist .neqv. case_calc%freezelist)) &
      & error stop trim(label)//' freeze-mask roundtrip mismatch'
    end if
    if (case_decoded_md%mtd(1)%com_bias .neqv. case_md%mtd(1)%com_bias .or. &
    & case_decoded_md%mtd(1)%com_mass_weighted .neqv. case_md%mtd(1)%com_mass_weighted .or. &
    & case_decoded_md%mtd(1)%com_factor /= case_md%mtd(1)%com_factor .or. &
    & case_decoded_md%mtd(1)%com_width /= case_md%mtd(1)%com_width .or. &
    & any(case_decoded_md%mtd(1)%atinclude .neqv. case_md%mtd(1)%atinclude)) then
      error stop trim(label)//' MTD state roundtrip mismatch'
    end if
    if (size(case_decoded_calc%calcs(1)%gff_fragments) /= &
    & size(case_calc%calcs(1)%gff_fragments)) error stop trim(label)//' fragment count mismatch'
    do case_i = 1,size(case_calc%calcs(1)%gff_fragments)
      if (len(case_decoded_calc%calcs(1)%gff_fragments(case_i)%value) /= &
      & len(case_calc%calcs(1)%gff_fragments(case_i)%value) .or. &
      & case_decoded_calc%calcs(1)%gff_fragments(case_i)%value /= &
      & case_calc%calcs(1)%gff_fragments(case_i)%value) then
        error stop trim(label)//' fragment string roundtrip mismatch'
      end if
    end do
    if (case_decoded_calc%nconstraints /= 1 .or. .not.allocated(case_decoded_calc%cons) .or. &
    & case_decoded_calc%cons(1)%n /= nat .or. &
    & any(case_decoded_calc%cons(1)%atms /= case_calc%cons(1)%atms) .or. &
    & any(case_decoded_calc%cons(1)%ref /= case_calc%cons(1)%ref) .or. &
    & any(case_decoded_calc%cons(1)%fc /= case_calc%cons(1)%fc)) then
      error stop trim(label)//' NCI wall roundtrip mismatch'
    end if
  end subroutine run_generic_case

  subroutine make_capsule_variant(source,destination,damage_endian)
    character(len=*),intent(in) :: source,destination
    logical,intent(in) :: damage_endian
    integer(int64) :: bytes_count
    integer(int8),allocatable :: bytes(:)
    integer :: input_unit,output_unit,local_io
    bytes_count=0_int64
    inquire(file=source,size=bytes_count,iostat=local_io)
    if (local_io /= 0) error stop 'cannot inspect test capsule'
    if (bytes_count < 41_int64) error stop 'test capsule is too short'
    allocate(bytes(int(bytes_count)))
    open(newunit=input_unit,file=source,access='stream',form='unformatted', &
    & status='old',action='read')
    read(input_unit) bytes
    close(input_unit)
    open(newunit=output_unit,file=destination,access='stream',form='unformatted', &
    & status='replace',action='write')
    if (damage_endian) then
      bytes(37) = ieor(bytes(37),1_int8)
      write(output_unit) bytes
    else
      write(output_unit) bytes(1:size(bytes)-1)
    end if
    close(output_unit)
  end subroutine make_capsule_variant

  subroutine make_v4_header_variant(source,destination)
    character(len=*),intent(in) :: source,destination
    integer(int64) :: bytes_count
    integer(int8),allocatable :: bytes(:)
    integer :: input_unit,output_unit,local_io

    bytes_count=0_int64
    inquire(file=source,size=bytes_count,iostat=local_io)
    if (local_io /= 0 .or. bytes_count < 36_int64) &
    & error stop 'cannot inspect capsule for V4 rejection test'
    allocate(bytes(int(bytes_count)))
    open(newunit=input_unit,file=source,access='stream',form='unformatted', &
    & status='old',action='read')
    read(input_unit) bytes
    close(input_unit)
    ! A V4 header has V4 magic followed by a native-endian int32 value of 4.
    bytes(28) = int(iachar('4'),int8)
    bytes(33:36) = 0_int8
    bytes(33) = 4_int8
    open(newunit=output_unit,file=destination,access='stream',form='unformatted', &
    & status='replace',action='write')
    write(output_unit) bytes
    close(output_unit)
  end subroutine make_v4_header_variant

  subroutine make_trailing_variant(source,destination)
    character(len=*),intent(in) :: source,destination
    integer(int64) :: bytes_count
    integer(int8),allocatable :: bytes(:)
    integer :: input_unit,output_unit,local_io
    bytes_count=0_int64
    inquire(file=source,size=bytes_count,iostat=local_io)
    if (local_io /= 0) error stop 'cannot inspect test capsule'
    if (bytes_count < 1_int64) error stop 'test capsule is empty'
    allocate(bytes(int(bytes_count)))
    open(newunit=input_unit,file=source,access='stream',form='unformatted', &
    & status='old',action='read')
    read(input_unit) bytes
    close(input_unit)
    open(newunit=output_unit,file=destination,access='stream',form='unformatted', &
    & status='replace',action='write')
    write(output_unit) bytes
    write(output_unit) 0_int8
    close(output_unit)
  end subroutine make_trailing_variant
end program capsule_roundtrip_driver
