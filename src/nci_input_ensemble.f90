!================================================================================!
! Multi-input NCI placements for the first iMTD-GC iteration.
!================================================================================!
module nci_input_ensemble
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use crest_parameters,only:wp,stdout,bohr
  use crest_data,only:systemdata
  use crest_calculator,only:calcdata
  use optimize_module,only:optimize_geometry
  use strucrd,only:coord,rdensembleparam,rdensemble,checkcoordtype,wrensemble
  use dynamics_module,only:mddata
  use ls_rmsd,only:rmsd
  implicit none
  private

  public :: load_nci_input_ensemble
  public :: write_nci_input_starts
  public :: write_nci_mtd_jobs

  real(wp),parameter :: transform_tolerance = 1.0e-7_wp
  real(wp),parameter :: fixed_atom_tolerance = 1.0e-8_wp

contains

  subroutine load_nci_input_ensemble(env,prepared_ref,mols,ninputs)
    type(systemdata),intent(inout) :: env
    type(coord),intent(in) :: prepared_ref
    type(coord),allocatable,intent(out) :: mols(:)
    integer,intent(out) :: ninputs

    type(coord) :: canonical_ref
    integer :: filetype,nat,natmax,iframe,iatom
    integer :: ntopo,nfixed,nelec
    integer,allocatable :: nats(:),ats(:,:),topo0(:),topoi(:)
    character(len=512),allocatable :: comments(:)
    real(wp),allocatable :: xyz_ang(:,:,:),raw_bohr(:,:,:),xyz_bohr(:,:,:)
    real(wp),allocatable :: rmsd_gradient(:,:)
    real(wp) :: rotation(3,3),raw_center(3),target_center(3),det
    real(wp) :: fit_rmsd,aligned_rmsd,max_fixed_delta,max_prepared_delta
    logical :: is_xyz

    ninputs = 1
    if (.not.env%NCI) return
    if (.not.allocated(env%inputcoords)) return

    call checkcoordtype(trim(env%inputcoords),filetype)
    is_xyz = filetype == 2
    if (.not.is_xyz) return

    call rdensembleparam(trim(env%inputcoords),natmax,ninputs)
    if (ninputs <= 1) then
      ninputs = 1
      return
    end if

    write(stdout,'(/,1x,a,i0,a)') 'Detected ',ninputs, &
      ' positional XYZ frames for the first NCI MTD stage.'

    allocate(nats(ninputs),ats(natmax,ninputs),comments(ninputs))
    allocate(xyz_ang(3,natmax,ninputs),source=0.0_wp)
    nats=0; ats=0
    call rdensemble(trim(env%inputcoords),natmax,ninputs,nats,ats,xyz_ang,comments)

    nat=nats(1)
    if (nat /= prepared_ref%nat .or. nat /= env%ref%nat) then
      error stop 'Multi-input NCI: frame 1 atom count differs from parsed input.'
    end if
    do iframe=1,ninputs
      if (nats(iframe) /= nat) error stop 'Multi-input NCI: atom counts differ.'
      if (any(ats(1:nat,iframe) /= ats(1:nat,1))) then
        error stop 'Multi-input NCI: elements or atom ordering differ.'
      end if
    end do
    if (any(ats(1:nat,1) /= prepared_ref%at)) then
      error stop 'Multi-input NCI: elements/order differ from parsed reference.'
    end if
    if (.not.all(ieee_is_finite(xyz_ang(:,1:nat,1:ninputs)))) then
      error stop 'Multi-input NCI: supplied frames contain non-finite coordinates.'
    end if

    call validate_global_selections(env,nat)
    nelec=sum(ats(1:nat,1))-env%chrg
    if (nelec < 0 .or. mod(nelec-env%uhf,2) /= 0) then
      error stop 'Multi-input NCI: incompatible charge and unpaired-electron count.'
    end if

    max_fixed_delta=0.0_wp
    nfixed=0
    if (allocated(env%calc%freezelist)) then
      nfixed=count(env%calc%freezelist)
      do iframe=2,ninputs
        do iatom=1,nat
          if (env%calc%freezelist(iatom)) then
            max_fixed_delta=max(max_fixed_delta, &
              maxval(abs(xyz_ang(:,iatom,iframe)-xyz_ang(:,iatom,1))))
          end if
        end do
      end do
      if (max_fixed_delta > fixed_atom_tolerance) then
        error stop 'Multi-input NCI: frozen-atom coordinates differ between frames.'
      end if
    end if

    ntopo=nat*(nat+1)/2
    allocate(topo0(ntopo),topoi(ntopo))
    call quicktopo(nat,ats(1:nat,1),xyz_ang(:,1:nat,1)/bohr,ntopo,topo0)
    do iframe=2,ninputs
      call quicktopo(nat,ats(1:nat,iframe),xyz_ang(:,1:nat,iframe)/bohr,ntopo,topoi)
      if (any(topoi /= topo0)) then
        error stop 'Multi-input NCI: supplied frames have incompatible topology.'
      end if
    end do

    ! The ordinary CREST input path has already written its canonical,
    ! transformed but unoptimized first frame to `coord`.  The NCI trial
    ! optimization may subsequently have changed prepared_ref internally, so
    ! raw frame 1 must be fitted to `coord`, not to prepared_ref.  One proper
    ! rigid transform is then applied unchanged to all supplied placements.
    call canonical_ref%open('coord')
    if (canonical_ref%nat /= nat) then
      error stop 'Multi-input NCI: canonical reference atom count differs.'
    end if
    if (any(canonical_ref%at /= ats(1:nat,1))) then
      error stop 'Multi-input NCI: canonical reference elements/order differ.'
    end if
    if (.not.all(ieee_is_finite(canonical_ref%xyz))) then
      error stop 'Multi-input NCI: canonical reference contains non-finite coordinates.'
    end if

    allocate(raw_bohr(3,nat,ninputs),xyz_bohr(3,nat,ninputs))
    allocate(rmsd_gradient(3,nat),source=0.0_wp)
    raw_bohr=xyz_ang(:,1:nat,:)/bohr
    call rmsd(nat,raw_bohr(:,:,1),canonical_ref%xyz,1,rotation, &
      raw_center,target_center,fit_rmsd,.false.,rmsd_gradient)

    if (.not.ieee_is_finite(fit_rmsd) .or. fit_rmsd < 0.0_wp) then
      error stop 'Multi-input NCI: canonical-reference rigid fit failed.'
    end if
    if (.not.all(ieee_is_finite(rotation)) .or. &
        .not.all(ieee_is_finite(raw_center)) .or. &
        .not.all(ieee_is_finite(target_center))) then
      error stop 'Multi-input NCI: canonical transform is non-finite.'
    end if

    det=determinant3(rotation)
    if (.not.ieee_is_finite(det) .or. abs(det-1.0_wp) > transform_tolerance) then
      error stop 'Multi-input NCI: canonical transform is not a proper rotation.'
    end if

    do iframe=1,ninputs
      do iatom=1,nat
        xyz_bohr(:,iatom,iframe)=matmul(rotation, &
          raw_bohr(:,iatom,iframe)-raw_center)+target_center
      end do
    end do

    aligned_rmsd=sqrt(sum((xyz_bohr(:,:,1)-canonical_ref%xyz)**2)/real(nat,wp))
    max_prepared_delta=maxval(abs(xyz_bohr(:,:,1)-canonical_ref%xyz))
    write(stdout,'(1x,a,es12.4,a)') 'Rigid-fit solver RMSD          : ', &
      fit_rmsd,' Bohr'
    write(stdout,'(1x,a,es12.4,a)') 'Canonical-reference RMSD       : ', &
      aligned_rmsd,' Bohr'
    write(stdout,'(1x,a,es12.4,a)') 'Canonical-reference max delta  : ', &
      max_prepared_delta,' Bohr'
    write(stdout,'(1x,a,f16.12)') 'Common transform determinant   : ',det
    if (aligned_rmsd > transform_tolerance .or. &
        max_prepared_delta > transform_tolerance) then
      error stop 'Multi-input NCI: transformed frame 1 differs from canonical reference.'
    end if
    call canonical_ref%deallocate()

    allocate(mols(ninputs))
    mols(1)=prepared_ref
    do iframe=2,ninputs
      mols(iframe)%nat=nat
      mols(iframe)%at=ats(1:nat,iframe)
      mols(iframe)%xyz=xyz_bohr(:,:,iframe)
      mols(iframe)%chrg=env%chrg
      mols(iframe)%uhf=env%uhf
    end do

    if (env%preopt) call preoptimize_additional_inputs(env,mols,topo0)

    write(stdout,'(1x,a,i0)') 'Validated NCI input structures : ',ninputs
    write(stdout,'(1x,a,i0)') 'Frozen atoms checked           : ',nfixed
    deallocate(rmsd_gradient,xyz_bohr,raw_bohr,topoi,topo0,xyz_ang,comments,ats,nats)
  end subroutine load_nci_input_ensemble

  subroutine validate_global_selections(env,nat)
    type(systemdata),intent(in) :: env
    integer,intent(in) :: nat
    if (.not.allocated(env%includeRMSD)) then
      error stop 'Multi-input NCI: RMSD/COM selection is not initialized.'
    end if
    if (size(env%includeRMSD) /= nat) then
      error stop 'Multi-input NCI: RMSD/COM selection has wrong size.'
    end if
    if (count(env%includeRMSD /= 0) == 0) then
      error stop 'Multi-input NCI: RMSD/COM selection is empty.'
    end if
    if (allocated(env%calc%freezelist)) then
      if (size(env%calc%freezelist) /= nat) then
        error stop 'Multi-input NCI: fixed-atom selection has wrong size.'
      end if
    end if
  end subroutine validate_global_selections

  subroutine preoptimize_additional_inputs(env,mols,input_topology)
    type(systemdata),intent(inout) :: env
    type(coord),intent(inout) :: mols(:)
    integer,intent(in) :: input_topology(:)
    type(coord) :: molopt
    type(calcdata) :: tmpcalc
    real(wp),allocatable :: grad(:,:)
    integer,allocatable :: topo(:)
    real(wp) :: energy
    integer :: iframe,io,ntopo,outer_threads,inner_threads
    logical :: pr,wr

    ntopo=size(input_topology)
    allocate(grad(3,mols(1)%nat),topo(ntopo))
    pr=.false.; wr=.false.
    call new_ompautoset(env,'serial',0,outer_threads,inner_threads)
    write(stdout,'(/,1x,a,i0,a)') 'Preoptimizing ',size(mols)-1, &
      ' additional NCI inputs with one GFN-FF core each.'
    do iframe=2,size(mols)
      tmpcalc=env%calc
      tmpcalc%optlev=-1
      molopt=mols(iframe)
      grad=0.0_wp
      call optimize_geometry(mols(iframe),molopt,tmpcalc,energy,grad,pr,wr,io)
      if (io /= 0) error stop 'Multi-input NCI: initial optimization failed.'
      call quicktopo(molopt%nat,molopt%at,molopt%xyz,ntopo,topo)
      if (any(topo /= input_topology)) then
        error stop 'Multi-input NCI: topology changed during initial optimization.'
      end if
      mols(iframe)=molopt
      call molopt%deallocate()
    end do
    call new_ompautoset(env,'max',0,outer_threads,inner_threads)
    deallocate(topo,grad)
  end subroutine preoptimize_additional_inputs

  pure real(wp) function determinant3(matrix)
    real(wp),intent(in) :: matrix(3,3)
    determinant3=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
      -matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
      +matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
  end function determinant3

  subroutine write_nci_input_starts(mols)
    type(coord),intent(in) :: mols(:)
    call wrensemble('crest_nci_input_starts.xyz',size(mols),mols)
  end subroutine write_nci_input_starts

  subroutine write_nci_mtd_jobs(mddats)
    type(mddata),intent(in) :: mddats(:)
    integer :: i,ich
    character(len=:),allocatable :: trajectory
    open(newunit=ich,file='crest_nci_mtd_jobs.tsv',status='replace',action='write')
    write(ich,'(a)') 'job_index'//achar(9)//'input_placement'//achar(9)// &
      'bias_configuration'//achar(9)//'trajectory'//achar(9)//'termination_status'
    do i=1,size(mddats)
      trajectory=''
      if (allocated(mddats(i)%trajectoryfile)) trajectory=mddats(i)%trajectoryfile
      write(ich,'(i0,a,i0,a,i0,a,a,a,i0)') mddats(i)%md_index,achar(9), &
        mddats(i)%input_structure_id,achar(9),mddats(i)%bias_configuration_id, &
        achar(9),trim(trajectory),achar(9),mddats(i)%termination_status
    end do
    close(ich)
  end subroutine write_nci_mtd_jobs
end module nci_input_ensemble
