program fragment_string_real_copy_oracle
  use iso_fortran_env,only:int8,output_unit
  use calc_type,only:alloc_string,calculation_settings,calcdata
  implicit none

  type(calculation_settings) :: level_a,level_b,level_c
  type(calcdata) :: calc_a,calc_b,calc_c,calc_add
  type(calcdata),allocatable :: workers(:)
  integer(int8),allocatable :: churn(:)
  integer :: i,j

  call set_level(level_a,101,[character(len=20) :: '1-808', &
  &                           '809-938','1000,1002-1017,2001'])
  call set_level(level_b,102,[character(len=12) :: '10-20', &
  &                           '30,32-45','901-938'])
  call set_level(level_c,103,[character(len=7) :: '1-808','809-938'])

  ! Real calculation_settings intrinsic assignment.
  calc_a%ncalculations = 0
  level_c = level_a
  call assert_level('real calculation_settings assignment',level_c,level_a)

  ! Real calcdata%add first-add path: self%calcs(1)=cal.
  call calc_a%add(level_a)
  call assert_calc_one('real first add',calc_a,level_a)

  ! Real additional-add path: callist(1:j)=self%calcs and callist(i)=cal.
  call calc_a%add(level_b)
  call assert_calc_two('real second add',calc_a,level_a,level_b)

  ! Real nested calcdata intrinsic assignment.
  calc_b = calc_a
  call assert_calc_two('real nested calcdata assignment',calc_b,level_a,level_b)

  ! Repeatedly overwrite already-allocated nested values while perturbing the
  ! heap, then exercise both add paths again in a fresh calcdata object.
  do i = 1,500
    allocate(churn(1024+mod(7919*i,131071)))
    churn = int(mod(i,127),int8)

    calc_c = calc_b
    call assert_calc_two('heap-churn nested copy A',calc_c,level_a,level_b)
    calc_b = calc_c
    call assert_calc_two('heap-churn nested copy B',calc_b,level_a,level_b)

    call calc_add%reset()
    call calc_add%add(level_a)
    call assert_calc_one('heap-churn first add',calc_add,level_a)
    call calc_add%add(level_b)
    call assert_calc_two('heap-churn additional add',calc_add,level_a,level_b)

    deallocate(churn)
  end do

  ! Array allocation with SOURCE invokes another real nested intrinsic copy.
  allocate(workers(16),source=calc_a)
  do i = 1,size(workers)
    call assert_calc_two('allocate(source=real calcdata)',workers(i),level_a,level_b)
  end do

  ! Prove the successful copies own independent character allocations.
  level_a%gff_fragments(2)%value = 'changed'
  calc_a%calcs(1)%gff_fragments(3)%value = 'also changed'
  call assert_calc_two('nested deep-copy independence',calc_b,level_c,level_b)
  do j = 1,size(workers)
    call assert_calc_two('worker deep-copy independence',workers(j),level_c,level_b)
  end do

  ! Main-program allocatables have static lifetime in GFortran; release them
  ! explicitly so LeakSanitizer audits the copy operations rather than
  ! reporting harmless exit-time retention.
  deallocate(workers)
  call calc_a%reset()
  call calc_b%reset()
  call calc_c%reset()
  call calc_add%reset()
  if (allocated(level_a%gff_fragments)) deallocate(level_a%gff_fragments)
  if (allocated(level_b%gff_fragments)) deallocate(level_b%gff_fragments)
  if (allocated(level_c%gff_fragments)) deallocate(level_c%gff_fragments)

  write(output_unit,'(a)') 'FRAGMENT_STRING_REAL_COPY_ORACLE PASS'

contains

  subroutine set_level(level,id,values)
    type(calculation_settings),intent(out) :: level
    integer,intent(in) :: id
    character(len=*),intent(in) :: values(:)
    integer :: k
    level%id = id
    allocate(level%gff_fragments(size(values)))
    do k = 1,size(values)
      level%gff_fragments(k)%value = values(k)
    end do
  end subroutine set_level

  subroutine assert_level(label,actual,expected)
    character(len=*),intent(in) :: label
    type(calculation_settings),intent(in) :: actual,expected
    integer :: k
    if (actual%id /= expected%id) error stop trim(label)//': id changed'
    if (allocated(actual%gff_fragments) .neqv. allocated(expected%gff_fragments)) &
    & error stop trim(label)//': array allocation changed'
    if (.not.allocated(actual%gff_fragments)) return
    if (size(actual%gff_fragments) /= size(expected%gff_fragments)) &
    & error stop trim(label)//': array size changed'
    do k = 1,size(actual%gff_fragments)
      if (allocated(actual%gff_fragments(k)%value) .neqv. &
      & allocated(expected%gff_fragments(k)%value)) &
      & error stop trim(label)//': element allocation changed'
      if (.not.allocated(actual%gff_fragments(k)%value)) cycle
      if (len(actual%gff_fragments(k)%value) /= len(expected%gff_fragments(k)%value)) &
      & error stop trim(label)//': element length changed'
      if (actual%gff_fragments(k)%value /= expected%gff_fragments(k)%value) &
      & error stop trim(label)//': element payload changed'
    end do
  end subroutine assert_level

  subroutine assert_calc_one(label,calc,expected_a)
    character(len=*),intent(in) :: label
    type(calcdata),intent(in) :: calc
    type(calculation_settings),intent(in) :: expected_a
    if (calc%ncalculations /= 1) error stop trim(label)//': wrong declared count'
    if (.not.allocated(calc%calcs)) error stop trim(label)//': calcs absent'
    if (size(calc%calcs) /= 1) error stop trim(label)//': wrong calcs size'
    call assert_level(label,calc%calcs(1),expected_a)
  end subroutine assert_calc_one

  subroutine assert_calc_two(label,calc,expected_a,expected_b)
    character(len=*),intent(in) :: label
    type(calcdata),intent(in) :: calc
    type(calculation_settings),intent(in) :: expected_a,expected_b
    if (calc%ncalculations /= 2) error stop trim(label)//': wrong declared count'
    if (.not.allocated(calc%calcs)) error stop trim(label)//': calcs absent'
    if (size(calc%calcs) /= 2) error stop trim(label)//': wrong calcs size'
    call assert_level(label,calc%calcs(1),expected_a)
    call assert_level(label,calc%calcs(2),expected_b)
  end subroutine assert_calc_two

end program fragment_string_real_copy_oracle
