program fragment_string_copy_oracle
  use iso_fortran_env,only:output_unit
  implicit none

  type :: alloc_string
    character(len=:),allocatable :: value
  end type alloc_string

  type :: calculation_settings_min
    integer :: id=0
    type(alloc_string),allocatable :: gff_fragments(:)
  end type calculation_settings_min

  type :: calcdata_min
    integer :: ncalculations=0
    type(calculation_settings_min),allocatable :: calcs(:)
  end type calcdata_min

  type(calculation_settings_min) :: level_a,level_b
  type(calcdata_min) :: calc_a,calc_b,calc_c
  type(calcdata_min),allocatable :: workers(:)
  integer :: i

  level_a%id=9
  allocate(level_a%gff_fragments(3))
  level_a%gff_fragments(1)%value='1-808'
  level_a%gff_fragments(2)%value='809-938'
  level_a%gff_fragments(3)%value='1000,1002-1017,2001'

  level_b=level_a
  call assert_level('direct calculation_settings assignment',level_b)

  calc_a%ncalculations=1
  allocate(calc_a%calcs(1))
  calc_a%calcs(1)=level_a
  call assert_calc('calcs element assignment',calc_a)

  calc_b=calc_a
  call assert_calc('nested calcdata assignment',calc_b)
  do i=1,20
    calc_c=calc_b
    call assert_calc('repeated nested calcdata assignment A',calc_c)
    calc_b=calc_c
    call assert_calc('repeated nested calcdata assignment B',calc_b)
  end do

  allocate(workers(8),source=calc_a)
  do i=1,size(workers)
    call assert_calc('allocate(source=calcdata) element',workers(i))
  end do

  ! Prove deep-copy independence rather than pointer aliasing.
  level_a%gff_fragments(2)%value='changed'
  call assert_level('deep-copy independence',level_b)
  call assert_calc('nested deep-copy independence',calc_b)
  do i=1,size(workers)
    call assert_calc('worker deep-copy independence',workers(i))
  end do

  ! Main-program allocatables have static lifetime in GFortran; release them
  ! explicitly so LeakSanitizer can distinguish the copy test from exit-time
  ! retained allocations.
  deallocate(workers)
  if (allocated(calc_a%calcs)) deallocate(calc_a%calcs)
  if (allocated(calc_b%calcs)) deallocate(calc_b%calcs)
  if (allocated(calc_c%calcs)) deallocate(calc_c%calcs)
  if (allocated(level_a%gff_fragments)) deallocate(level_a%gff_fragments)
  if (allocated(level_b%gff_fragments)) deallocate(level_b%gff_fragments)

  write(output_unit,'(a)') 'FRAGMENT_STRING_COPY_ORACLE PASS'

contains

  subroutine assert_level(label,level)
    character(len=*),intent(in) :: label
    type(calculation_settings_min),intent(in) :: level
    if (.not.allocated(level%gff_fragments)) error stop trim(label)//': unallocated array'
    if (size(level%gff_fragments) /= 3) error stop trim(label)//': wrong array size'
    if (.not.allocated(level%gff_fragments(1)%value)) error stop trim(label)//': value 1 absent'
    if (.not.allocated(level%gff_fragments(2)%value)) error stop trim(label)//': value 2 absent'
    if (.not.allocated(level%gff_fragments(3)%value)) error stop trim(label)//': value 3 absent'
    if (level%gff_fragments(1)%value /= '1-808') error stop trim(label)//': value 1 changed'
    if (level%gff_fragments(2)%value /= '809-938') error stop trim(label)//': value 2 changed'
    if (level%gff_fragments(3)%value /= '1000,1002-1017,2001') &
    & error stop trim(label)//': value 3 changed'
  end subroutine assert_level

  subroutine assert_calc(label,calc)
    character(len=*),intent(in) :: label
    type(calcdata_min),intent(in) :: calc
    if (calc%ncalculations /= 1) error stop trim(label)//': wrong declared count'
    if (.not.allocated(calc%calcs)) error stop trim(label)//': calcs absent'
    if (size(calc%calcs) /= 1) error stop trim(label)//': wrong calcs size'
    call assert_level(label,calc%calcs(1))
  end subroutine assert_calc

end program fragment_string_copy_oracle
