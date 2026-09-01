program nci_wall_lifecycle_driver
  use crest_parameters,only:wp
  use crest_calculator,only:calcdata,constraint
  implicit none

  type(calcdata) :: calc,copied
  type(constraint) :: user_constraint,auto_wall
  type(constraint),allocatable :: detached(:),none(:)
  integer :: i

  calc%nfreeze=2
  allocate(calc%freezelist(4),source=.false.)
  calc%freezelist(1:2)=.true.

  user_constraint%type=1
  user_constraint%n=2
  allocate(user_constraint%atms(2),source=[3,4])
  allocate(user_constraint%ref(1),source=2.0_wp)
  allocate(user_constraint%fc(1),source=0.5_wp)
  call user_constraint%addfreeze(calc%freezelist)
  call calc%add(user_constraint)

  auto_wall%type=5
  auto_wall%n=4
  auto_wall%auto_nci_wall=.true.
  allocate(auto_wall%atms(4),auto_wall%ref(3),auto_wall%fc(2))
  do i=1,4
    auto_wall%atms(i)=i
  end do
  auto_wall%ref=[10.0_wp,11.0_wp,12.0_wp]
  auto_wall%fc=[298.15_wp,50.0_wp]
  call auto_wall%addfreeze(calc%freezelist)
  call calc%add(auto_wall)

  if (calc%nconstraints/=2) error stop 'setup constraint count mismatch'
  copied=calc
  call copied%detach_auto_nci_walls(detached)
  if (size(detached)/=1) error stop 'automatic wall was not detached'
  if (copied%nconstraints/=1) error stop 'user constraint count changed on detach'
  if (copied%cons(1)%auto_nci_wall) error stop 'automatic wall remained attached'
  if (copied%cons(1)%type/=user_constraint%type) error stop 'user constraint order changed'
  if (associated(detached(1)%freezeptr)) error stop 'detached wall retained freeze pointer'
  if (detached(1)%frozenatms) error stop 'detached wall retained frozen-atoms state'
  if (.not.associated(copied%cons(1)%freezeptr)) then
    error stop 'retained user constraint was not rebound'
  end if
  calc%freezelist(1)=.false.
  if (.not.copied%cons(1)%freezeptr(1)) then
    error stop 'retained user constraint points at source freeze storage'
  end if
  calc%freezelist(1)=.true.
  copied%freezelist(1)=.false.
  if (copied%cons(1)%freezeptr(1)) then
    error stop 'retained user constraint does not point at copied freeze storage'
  end if
  copied%freezelist(1)=.true.

  call copied%attach_auto_nci_walls(detached)
  if (copied%nconstraints/=2) error stop 'automatic wall was not restored'
  if (.not.copied%cons(2)%auto_nci_wall) error stop 'restored wall lost provenance'
  if (.not.associated(copied%cons(2)%freezeptr)) error stop 'restored wall was not rebound'
  if (.not.copied%cons(2)%frozenatms) error stop 'restored wall lost frozen-atoms state'
  if (any(copied%cons(2)%ref/=auto_wall%ref)) error stop 'restored wall radii changed'
  if (any(copied%cons(2)%atms/=auto_wall%atms)) error stop 'restored wall atoms changed'
  copied%freezelist(1)=.false.
  if (copied%cons(2)%freezeptr(1)) error stop 'restored wall points at stale freeze storage'
  copied%freezelist(1)=.true.
  if (.not.copied%cons(2)%freezeptr(1)) error stop 'restored wall missed mask restoration'

  call copied%detach_auto_nci_walls(detached)
  call copied%detach_auto_nci_walls(none)
  if (size(none)/=0) error stop 'zero-wall detach returned a wall'
  if (copied%nconstraints/=1) error stop 'zero-wall detach changed user constraints'
  call copied%attach_auto_nci_walls(detached)

  call copied%detach_auto_nci_walls(detached)
  deallocate(copied%cons)
  copied%nconstraints=0
  allocate(copied%cons(0))
  call copied%attach_auto_nci_walls(detached)
  if (copied%nconstraints/=1) error stop 'zero-sized array restore failed'
  if (.not.allocated(copied%cons)) error stop 'zero-sized array restore left no array'
  if (size(copied%cons)/=1) error stop 'zero-sized array restore has wrong size'

  write(*,'(a)') 'NCI_WALL_LIFECYCLE_PASS'
end program nci_wall_lifecycle_driver
