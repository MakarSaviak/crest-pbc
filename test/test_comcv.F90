module test_comcv
  use testdrive, only : new_unittest, unittest_type, error_type, test_failed
  use crest_parameters, only : wp
  use strucrd, only : coord
  use metadynamics_module, only : mtdpot, cv_rmsd, calc_com_mtd
  use ls_rmsd, only : rmsd
  use atmasses, only : ams
  use omp_lib, only : omp_set_dynamic, omp_set_num_threads
  implicit none
  private

  public :: collect_comcv

  real(wp), parameter :: etol = 5.0e-13_wp
  real(wp), parameter :: gtol = 2.0e-9_wp

contains

  subroutine collect_comcv(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("rigid translation: RMSD zero, COM active", test_rigid_translation), &
      new_unittest("COM disabled is exact no-op", test_disabled), &
      new_unittest("mass-weighted energy and gradient", test_mass_weighted), &
      new_unittest("finite-difference COM gradient", test_finite_difference), &
      new_unittest("uniform weights and multiple hills", test_uniform_multiple), &
      new_unittest("one/two-thread agreement", test_thread_agreement) &
    ]
  end subroutine collect_comcv

  subroutine init_case(mol,pot,nref,mass_weighted)
    type(coord), intent(out) :: mol
    type(mtdpot), intent(out) :: pot
    integer, intent(in) :: nref
    logical, intent(in) :: mass_weighted
    real(wp), parameter :: ref(3,5) = reshape([ &
      0.0_wp, 0.0_wp, 0.0_wp, &
      1.0_wp, 0.0_wp, 0.0_wp, &
      0.0_wp, 1.2_wp, 0.1_wp, &
      0.2_wp, 0.1_wp, 1.4_wp, &
      1.1_wp, 1.0_wp, 0.7_wp ], [3,5])
    real(wp), parameter :: shift(3) = [0.40_wp,-0.20_wp,0.30_wp]
    integer :: i

    mol%nat = 5
    allocate(mol%at(5), mol%xyz(3,5))
    mol%at = [6,1,6,7,8]
    mol%xyz = ref
    do i = 2,5
      mol%xyz(:,i) = mol%xyz(:,i)+shift
    end do

    pot%mtdtype = cv_rmsd
    pot%ncur = nref
    pot%nmax = nref
    pot%maxsave = nref
    pot%com_bias = .true.
    pot%com_factor = 0.0010_wp
    pot%com_width = 0.05_wp
    pot%com_mass_weighted = mass_weighted
    pot%damp = 1.0_wp
    allocate(pot%atinclude(5), source=.false.)
    pot%atinclude(2:5) = .true.
    allocate(pot%cvxyz(3,5,nref), source=0.0_wp)
    do i = 1,nref
      pot%cvxyz(:,:,i) = ref
    end do
  end subroutine init_case

  subroutine expected_one(mol,pot,iref,e,g)
    type(coord), intent(in) :: mol
    type(mtdpot), intent(in) :: pot
    integer, intent(in) :: iref
    real(wp), intent(out) :: e
    real(wp), intent(out) :: g(3,mol%nat)
    real(wp) :: rnow(3),rref(3),dr(3),w,wsum,scale
    integer :: j

    rnow = 0.0_wp
    rref = 0.0_wp
    wsum = 0.0_wp
    g = 0.0_wp
    do j = 1,mol%nat
      if (.not.pot%atinclude(j)) cycle
      if (pot%com_mass_weighted) then
        w = ams(mol%at(j))
      else
        w = 1.0_wp
      end if
      wsum = wsum+w
      rnow = rnow+w*mol%xyz(:,j)
      rref = rref+w*pot%cvxyz(:,j,iref)
    end do
    rnow = rnow/wsum
    rref = rref/wsum
    dr = rnow-rref
    scale = 1.0_wp
    if (iref == pot%ncur) scale = pot%damp
    e = pot%com_factor*scale*exp(-pot%com_width*dot_product(dr,dr))
    do j = 1,mol%nat
      if (.not.pot%atinclude(j)) cycle
      if (pot%com_mass_weighted) then
        w = ams(mol%at(j))/wsum
      else
        w = 1.0_wp/wsum
      end if
      g(:,j) = w*(-2.0_wp*pot%com_width*e*dr)
    end do
  end subroutine expected_one


  subroutine test_rigid_translation(error)
    type(error_type), allocatable, intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp) :: ecom,gcom(3,5)
    real(wp) :: u(3,3),xcenter(3),ycenter(3),rmsdval,grmsd(3,4)

    call init_case(mol,pot,1,.true.)
    call rmsd(4,mol%xyz(:,2:5),pot%cvxyz(:,2:5,1),1,u,xcenter,ycenter, &
      rmsdval,.false.,grmsd)
    call calc_com_mtd(mol,pot,ecom,gcom)

    ! The quaternion implementation returns a residual of order sqrt(epsilon)
    ! for an exact rigid translation, so use a tolerance just above that scale.
    if (abs(rmsdval) > 5.0e-8_wp) then
      call test_failed(error,"rigid translation was not removed by aligned RMSD")
      return
    end if
    if (ecom <= 0.0_wp .or. maxval(abs(gcom(:,2:5))) <= 0.0_wp) then
      call test_failed(error,"rigid translation did not activate COM bias")
      return
    end if
    if (any(gcom(:,1) /= 0.0_wp)) then
      call test_failed(error,"COM gradient leaked onto unselected atom")
    end if
  end subroutine test_rigid_translation

  subroutine test_disabled(error)
    type(error_type), allocatable, intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp) :: e,g(3,5)

    call init_case(mol,pot,1,.true.)
    pot%com_bias = .false.
    call calc_com_mtd(mol,pot,e,g)
    if (e /= 0.0_wp .or. any(g /= 0.0_wp)) then
      call test_failed(error,"disabled COM path changed energy or gradient")
    end if
  end subroutine test_disabled

  subroutine test_mass_weighted(error)
    type(error_type), allocatable, intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp) :: e,eref,g(3,5),gref(3,5),gsum(3),expected_sum(3)

    call init_case(mol,pot,1,.true.)
    call calc_com_mtd(mol,pot,e,g)
    call expected_one(mol,pot,1,eref,gref)

    if (abs(e-eref) > etol) then
      call test_failed(error,"mass-weighted COM energy mismatch")
      return
    end if
    if (maxval(abs(g-gref)) > etol) then
      call test_failed(error,"mass-weighted COM gradient mismatch")
      return
    end if
    if (any(g(:,1) /= 0.0_wp)) then
      call test_failed(error,"unselected atom received COM gradient")
      return
    end if
    gsum = sum(g,dim=2)
    expected_sum = sum(gref,dim=2)
    if (maxval(abs(gsum-expected_sum)) > etol) then
      call test_failed(error,"selected gradient sum mismatch")
    end if
  end subroutine test_mass_weighted

  subroutine test_finite_difference(error)
    type(error_type), allocatable, intent(out) :: error
    type(coord) :: mol,mp,mm
    type(mtdpot) :: pot
    real(wp) :: e,g(3,5),ep,em,gwork(3,5),fd,h
    integer :: atom,axis

    call init_case(mol,pot,1,.true.)
    call calc_com_mtd(mol,pot,e,g)
    h = 1.0e-6_wp
    do atom = 2,5
      do axis = 1,3
        mp = mol
        mm = mol
        mp%xyz(axis,atom) = mp%xyz(axis,atom)+h
        mm%xyz(axis,atom) = mm%xyz(axis,atom)-h
        call calc_com_mtd(mp,pot,ep,gwork)
        call calc_com_mtd(mm,pot,em,gwork)
        fd = (ep-em)/(2.0_wp*h)
        if (abs(fd-g(axis,atom)) > gtol) then
          call test_failed(error,"central finite-difference gradient mismatch")
          return
        end if
      end do
    end do
  end subroutine test_finite_difference

  subroutine test_uniform_multiple(error)
    type(error_type), allocatable, intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp) :: e,e1,e2,g(3,5),g1(3,5),g2(3,5)
    integer :: j

    call init_case(mol,pot,2,.false.)
    pot%damp = 0.35_wp
    do j = 2,5
      pot%cvxyz(:,j,2) = pot%cvxyz(:,j,2)+[0.10_wp,0.25_wp,-0.15_wp]
    end do
    call calc_com_mtd(mol,pot,e,g)
    call expected_one(mol,pot,1,e1,g1)
    call expected_one(mol,pot,2,e2,g2)
    if (abs(e-(e1+e2)) > etol) then
      call test_failed(error,"multiple COM hills did not sum correctly")
      return
    end if
    if (maxval(abs(g-(g1+g2))) > etol) then
      call test_failed(error,"multiple COM gradients did not sum correctly")
      return
    end if
    if (any(g(:,1) /= 0.0_wp)) then
      call test_failed(error,"uniform-weight COM touched unselected atom")
    end if
  end subroutine test_uniform_multiple

  subroutine test_thread_agreement(error)
    type(error_type), allocatable, intent(out) :: error
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp) :: e1,e2,g1(3,5),g2(3,5)
    integer :: i,j

    call init_case(mol,pot,4,.true.)
    pot%damp = 0.7_wp
    do i = 2,4
      do j = 2,5
        pot%cvxyz(:,j,i) = pot%cvxyz(:,j,i) + &
          real(i-1,wp)*[0.05_wp,-0.03_wp,0.02_wp]
      end do
    end do
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(1)
    call calc_com_mtd(mol,pot,e1,g1)
    call omp_set_num_threads(2)
    call calc_com_mtd(mol,pot,e2,g2)
    if (abs(e1-e2) > 1.0e-12_wp .or. maxval(abs(g1-g2)) > 1.0e-12_wp) then
      call test_failed(error,"one-thread and two-thread COM results disagree")
    end if
  end subroutine test_thread_agreement

end module test_comcv
