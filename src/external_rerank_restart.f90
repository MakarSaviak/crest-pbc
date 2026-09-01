!================================================================================!
! Explicit staged restart after external ranking of iteration-1 source conformers.
!================================================================================!
module external_rerank_restart
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use, intrinsic :: iso_fortran_env,only:iostat_end
  use crest_parameters,only:wp,stdout,bohr
  use crest_data,only:systemdata
  use strucrd,only:coord,rdensembleparam,rdensemble,get_atlist
  implicit none
  private

  character(len=*),parameter :: checkpoint_file='crest.external_rerank'
  integer,parameter :: nci_wall_protocol_one_wall=1
  integer,parameter :: max_external_rerank_biases=65536

  type,public :: external_rerank_state
    integer :: nci_wall_protocol=0
    integer :: completed_mtd=0
    integer :: target_mtd=0
    integer :: nmetadyn=0
    real(wp) :: elowest=0.0_wp
    real(wp) :: eprivious=0.0_wp
    character(len=512) :: settings_file=''
    character(len=512) :: constraints_file=''
    character(len=512) :: archive_file=''
    real(wp),allocatable :: metadfac(:),metadexp(:)
  end type external_rerank_state

  public :: external_rerank_checkpoint_exists
  public :: write_external_rerank_checkpoint
  public :: read_external_rerank_checkpoint
  public :: restore_external_rerank_state
  public :: validate_external_rerank_seed
  public :: remove_external_rerank_checkpoint

contains

  logical function external_rerank_checkpoint_exists() result(exists)
    inquire(file=checkpoint_file,exist=exists)
  end function external_rerank_checkpoint_exists

  subroutine write_external_rerank_checkpoint(env,archive_file,completed_mtd)
    type(systemdata),intent(in) :: env
    character(len=*),intent(in) :: archive_file
    integer,intent(in) :: completed_mtd
    integer :: ich,i,io,close_io
    logical :: ex

    if (.not.env%NCI) then
      error stop '**ERROR** external reranking is supported only for NCI searches'
    end if
    if (completed_mtd/=1) then
      error stop '**ERROR** external-rerank checkpoint must follow exactly MTD iteration 1'
    end if
    if (env%Maxrestart<=completed_mtd) then
      error stop '**ERROR** external reranking requires an MTD iteration after the checkpoint'
    end if
    if (.not.allocated(env%input_settings_file)) then
      error stop '**ERROR** external reranking requires an original TOML settings file'
    end if
    if (len_trim(env%input_settings_file)==0) then
      error stop '**ERROR** external reranking requires a nonempty TOML settings path'
    end if
    inquire(file=trim(archive_file),exist=ex)
    if (.not.ex) error stop '**ERROR** iteration archive missing before external checkpoint'
    if (env%nmetadyn<1 .or. env%nmetadyn>max_external_rerank_biases) then
      error stop '**ERROR** invalid MTD bias count at external checkpoint'
    end if
    if (.not.allocated(env%metadfac)) then
      error stop '**ERROR** MTD bias factors are absent at external checkpoint'
    end if
    if (.not.allocated(env%metadexp)) then
      error stop '**ERROR** MTD bias exponents are absent at external checkpoint'
    end if
    if (size(env%metadfac)<env%nmetadyn) then
      error stop '**ERROR** MTD bias-factor array is too short at external checkpoint'
    end if
    if (size(env%metadexp)<env%nmetadyn) then
      error stop '**ERROR** MTD bias-exponent array is too short at external checkpoint'
    end if
    if (.not.ieee_is_finite(env%elowest) .or. .not.ieee_is_finite(env%eprivious)) then
      error stop '**ERROR** nonfinite energy state at external checkpoint'
    end if
    if (.not.all(ieee_is_finite(env%metadfac(1:env%nmetadyn)))) then
      error stop '**ERROR** nonfinite MTD bias factor at external checkpoint'
    end if
    if (.not.all(ieee_is_finite(env%metadexp(1:env%nmetadyn)))) then
      error stop '**ERROR** nonfinite MTD bias exponent at external checkpoint'
    end if
    call require_exact_live_nci_wall(env)

    open(newunit=ich,file=checkpoint_file,status='replace',action='write',iostat=io)
    if (io/=0) error stop '**ERROR** cannot create external-rerank checkpoint'
    write(ich,'(a)',iostat=io) 'CREST_EXTERNAL_RERANK_V3'
    if (io==0) write(ich,'(i0)',iostat=io) nci_wall_protocol_one_wall
    if (io==0) write(ich,'(i0)',iostat=io) completed_mtd
    if (io==0) write(ich,'(i0)',iostat=io) env%Maxrestart
    if (io==0) write(ich,'(i0)',iostat=io) env%nmetadyn
    if (io==0) write(ich,'(es26.17)',iostat=io) env%elowest
    if (io==0) write(ich,'(es26.17)',iostat=io) env%eprivious
    if (io==0) write(ich,'(a)',iostat=io) trim(env%input_settings_file)
    if (io==0) write(ich,'(a)',iostat=io) trim(env%constraints)
    if (io==0) write(ich,'(a)',iostat=io) trim(archive_file)
    do i=1,env%nmetadyn
      if (io==0) write(ich,'(2es26.17)',iostat=io) env%metadfac(i),env%metadexp(i)
    end do
    close(ich,iostat=close_io)
    if (io/=0 .or. close_io/=0) error stop '**ERROR** failed to write external-rerank checkpoint'
  end subroutine write_external_rerank_checkpoint

  subroutine read_external_rerank_checkpoint(state)
    type(external_rerank_state),intent(out) :: state
    integer :: ich,io,i,allocstat
    character(len=1024) :: record
    logical :: ex

    inquire(file=checkpoint_file,exist=ex)
    if (.not.ex) error stop '**ERROR** external-rerank checkpoint not found'
    open(newunit=ich,file=checkpoint_file,status='old',action='read',iostat=io)
    if (io/=0) error stop '**ERROR** cannot open external-rerank checkpoint'
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** unreadable external-rerank checkpoint format'
    end if
    if (trim(record)/='CREST_EXTERNAL_RERANK_V3') then
      close(ich)
      error stop '**ERROR** invalid external-rerank checkpoint format'
    end if
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank checkpoint protocol'
    end if
    read(record,*,iostat=io) state%nci_wall_protocol
    if (io/=0) then
      close(ich)
      error stop '**ERROR** invalid external-rerank checkpoint protocol'
    end if
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank completed-iteration state'
    end if
    read(record,*,iostat=io) state%completed_mtd
    if (io/=0) then
      close(ich)
      error stop '**ERROR** invalid external-rerank completed-iteration state'
    end if
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank target-iteration state'
    end if
    read(record,*,iostat=io) state%target_mtd
    if (io/=0) then
      close(ich)
      error stop '**ERROR** invalid external-rerank target-iteration state'
    end if
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank bias count'
    end if
    read(record,*,iostat=io) state%nmetadyn
    if (io/=0) then
      close(ich)
      error stop '**ERROR** invalid external-rerank bias count'
    end if
    if (state%nmetadyn<1 .or. state%nmetadyn>max_external_rerank_biases) then
      close(ich)
      error stop '**ERROR** invalid external-rerank bias count'
    end if
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank lowest-energy state'
    end if
    read(record,*,iostat=io) state%elowest
    if (io/=0) then
      close(ich)
      error stop '**ERROR** invalid external-rerank lowest-energy state'
    end if
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank previous-energy state'
    end if
    read(record,*,iostat=io) state%eprivious
    if (io/=0) then
      close(ich)
      error stop '**ERROR** invalid external-rerank previous-energy state'
    end if
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank settings path'
    end if
    if (len_trim(record)>len(state%settings_file)) then
      close(ich)
      error stop '**ERROR** external-rerank settings path is too long'
    end if
    state%settings_file=trim(record)
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank constraints path'
    end if
    if (len_trim(record)>len(state%constraints_file)) then
      close(ich)
      error stop '**ERROR** external-rerank constraints path is too long'
    end if
    state%constraints_file=trim(record)
    read(ich,'(a)',iostat=io) record
    if (io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank archive path'
    end if
    if (len_trim(record)>len(state%archive_file)) then
      close(ich)
      error stop '**ERROR** external-rerank archive path is too long'
    end if
    state%archive_file=trim(record)
    if (state%nci_wall_protocol/=nci_wall_protocol_one_wall) then
      close(ich)
      error stop '**ERROR** external-rerank checkpoint uses an incompatible NCI wall protocol'
    end if
    if (state%completed_mtd/=1 .or. state%target_mtd<=state%completed_mtd) then
      close(ich)
      error stop '**ERROR** invalid external-rerank iteration state'
    end if
    if (.not.ieee_is_finite(state%elowest) .or. .not.ieee_is_finite(state%eprivious)) then
      close(ich)
      error stop '**ERROR** nonfinite external-rerank energy state'
    end if
    if (len_trim(state%settings_file)==0 .or. len_trim(state%archive_file)==0) then
      close(ich)
      error stop '**ERROR** external-rerank checkpoint lacks required paths'
    end if
    allocate(state%metadfac(state%nmetadyn),stat=allocstat)
    if (allocstat/=0) then
      close(ich)
      error stop '**ERROR** cannot allocate external-rerank bias-factor state'
    end if
    allocate(state%metadexp(state%nmetadyn),stat=allocstat)
    if (allocstat/=0) then
      deallocate(state%metadfac)
      close(ich)
      error stop '**ERROR** cannot allocate external-rerank bias-exponent state'
    end if
    do i=1,state%nmetadyn
      read(ich,'(a)',iostat=io) record
      if (io/=0) then
        close(ich)
        error stop '**ERROR** truncated external-rerank bias state'
      end if
      read(record,*,iostat=io) state%metadfac(i),state%metadexp(i)
      if (io/=0) then
        close(ich)
        error stop '**ERROR** invalid external-rerank bias state'
      end if
    end do
    read(ich,'(a)',iostat=io) record
    close(ich)
    if (io/=iostat_end) then
      error stop '**ERROR** external-rerank checkpoint has trailing or unreadable records'
    end if
    if (.not.all(ieee_is_finite(state%metadfac))) then
      error stop '**ERROR** nonfinite external-rerank bias factor'
    end if
    if (.not.all(ieee_is_finite(state%metadexp))) then
      error stop '**ERROR** nonfinite external-rerank bias exponent'
    end if
  end subroutine read_external_rerank_checkpoint

  subroutine restore_external_rerank_state(env,state)
    type(systemdata),intent(inout) :: env
    type(external_rerank_state),intent(in) :: state
    if (.not.env%NCI) then
      error stop '**ERROR** external-rerank restart is supported only for NCI searches'
    end if
    if (state%nci_wall_protocol/=nci_wall_protocol_one_wall .or. state%completed_mtd/=1) then
      error stop '**ERROR** invalid external-rerank restart protocol state'
    end if
    if (state%target_mtd<=state%completed_mtd) then
      error stop '**ERROR** invalid external-rerank restart target iteration'
    end if
    if (state%nmetadyn<1 .or. state%nmetadyn>max_external_rerank_biases) then
      error stop '**ERROR** invalid external-rerank restart bias count'
    end if
    if (.not.allocated(state%metadfac)) then
      error stop '**ERROR** external-rerank restart bias factors are absent'
    end if
    if (.not.allocated(state%metadexp)) then
      error stop '**ERROR** external-rerank restart bias exponents are absent'
    end if
    if (size(state%metadfac)/=state%nmetadyn) then
      error stop '**ERROR** external-rerank restart bias-factor size mismatch'
    end if
    if (size(state%metadexp)/=state%nmetadyn) then
      error stop '**ERROR** external-rerank restart bias-exponent size mismatch'
    end if
    if (.not.ieee_is_finite(state%elowest) .or. .not.ieee_is_finite(state%eprivious)) then
      error stop '**ERROR** nonfinite external-rerank restart energy state'
    end if
    if (.not.all(ieee_is_finite(state%metadfac)) .or. &
      & .not.all(ieee_is_finite(state%metadexp))) then
      error stop '**ERROR** nonfinite external-rerank restart bias state'
    end if
    call require_exact_live_nci_wall(env)

    ! Only the compact MTD state is restored.  The NCI wall has already been
    ! reparsed and ensemble-finalized in the normal setup path and is never
    ! serialized in this checkpoint.
    call env%deallocate()
    call env%allocate(state%nmetadyn)
    env%metadfac=state%metadfac
    env%metadexp=state%metadexp
    env%metadynset=.true.
    env%Maxrestart=state%target_mtd
    env%nmetadyn=state%nmetadyn
    env%elowest=state%elowest
    env%eprivious=state%eprivious
    if(len_trim(state%constraints_file)>0) env%constraints=trim(state%constraints_file)
  end subroutine restore_external_rerank_state

  subroutine require_exact_live_nci_wall(env)
    type(systemdata),intent(in) :: env

    if (.not.env%NCI) then
      error stop '**ERROR** external reranking requires an NCI calculation'
    end if
    if (env%calc%nconstraints/=1) then
      error stop '**ERROR** external reranking requires exactly one automatic NCI wall'
    end if
    if (.not.allocated(env%calc%cons)) then
      error stop '**ERROR** external-rerank automatic NCI wall is absent'
    end if
    if (size(env%calc%cons)/=1) then
      error stop '**ERROR** external-rerank NCI wall array/count mismatch'
    end if
    if (.not.env%calc%cons(1)%is_exact_auto_nci_wall(env%ref%nat)) then
      error stop '**ERROR** external-rerank automatic NCI wall schema is invalid'
    end if
  end subroutine require_exact_live_nci_wall

  subroutine validate_external_rerank_seed(env,state,matched_frame)
    type(systemdata),intent(in) :: env
    type(external_rerank_state),intent(in) :: state
    integer,intent(out),optional :: matched_frame
    type(coord) :: seed
    integer :: natmax,nall,i,match
    integer,allocatable :: nats(:),ats(:,:)
    real(wp),allocatable :: xyz_ang(:,:,:)
    character(len=512),allocatable :: comments(:)
    character(len=:),allocatable :: fragments(:)
    logical :: have_fragments,ok,ex

    inquire(file=trim(state%archive_file),exist=ex)
    if(.not.ex) error stop '**ERROR** iteration-1 archive not found for external restart'
    call env%ref%to(seed)
    if(allocated(seed%lat)) error stop '**ERROR** external-rerank seed must be finite/non-periodic'

    call rdensembleparam(trim(state%archive_file),natmax,nall)
    if(nall<1) error stop '**ERROR** iteration-1 archive is empty'
    allocate(nats(nall),ats(natmax,nall),xyz_ang(3,natmax,nall),comments(nall))
    nats=0; ats=0; xyz_ang=0.0_wp
    call rdensemble(trim(state%archive_file),natmax,nall,nats,ats,xyz_ang,comments)
    call get_consistent_gff_fragments(env,fragments,have_fragments)

    match=0
    do i=1,nall
      if(nats(i)/=seed%nat) cycle
      if(any(ats(1:seed%nat,i)/=seed%at)) cycle
      if(have_fragments) then
        ok=fragment_topology_equal(seed%nat,seed%at,seed%xyz, &
          xyz_ang(:,1:seed%nat,i)/bohr,fragments)
      else
        ok=full_topology_equal(seed%nat,seed%at,seed%xyz, &
          xyz_ang(:,1:seed%nat,i)/bohr)
      end if
      if(ok) then
        match=i
        exit
      end if
    end do
    if(match==0) error stop '**ERROR** external-rerank seed topology does not match iteration-1 archive'

    write(stdout,'(/,1x,a)') 'External-rerank seed validation:'
    write(stdout,'(3x,a,i0)') 'compatible iteration-1 frame: ',match
    if(have_fragments) then
      write(stdout,'(3x,a,i0)') 'GFN-FF fragment groups used: ',size(fragments)
      write(stdout,'(3x,a)') 'inter-fragment contacts excluded from topology matching'
    else
      write(stdout,'(3x,a)') 'full-system topology compared (no fragments supplied)'
    end if
    write(stdout,'(3x,a)') 'coordinates and energies were not required to match'
    if(present(matched_frame)) matched_frame=match

    call seed%deallocate()
    deallocate(comments,xyz_ang,ats,nats)
    if(allocated(fragments)) deallocate(fragments)
  end subroutine validate_external_rerank_seed

  logical function full_topology_equal(nat,at,xyz_a,xyz_b) result(equal)
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz_a(3,nat),xyz_b(3,nat)
    integer :: ntopo
    integer,allocatable :: topo_a(:),topo_b(:)
    ntopo=nat*(nat+1)/2
    allocate(topo_a(ntopo),topo_b(ntopo))
    call quicktopo(nat,at,xyz_a,ntopo,topo_a)
    call quicktopo(nat,at,xyz_b,ntopo,topo_b)
    equal=all(topo_a==topo_b)
    deallocate(topo_b,topo_a)
  end function full_topology_equal

  logical function fragment_topology_equal(nat,at,xyz_a,xyz_b,fragments) result(equal)
    integer,intent(in) :: nat,at(nat)
    real(wp),intent(in) :: xyz_a(3,nat),xyz_b(3,nat)
    character(len=*),intent(in) :: fragments(:)
    integer,allocatable :: frag_ids(:),subat(:),topo_a(:),topo_b(:)
    real(wp),allocatable :: sub_a(:,:),sub_b(:,:)
    logical,allocatable :: selected(:)
    integer :: f,nf,nsub,i,k,ntopo

    allocate(frag_ids(nat),source=0)
    do f=1,size(fragments)
      call get_atlist(nat,selected,trim(fragments(f)),at)
      where(selected) frag_ids=f
      deallocate(selected)
    end do
    nf=maxval(frag_ids)
    equal=.false.
    do f=0,nf
      nsub=count(frag_ids==f)
      if(nsub<=1) cycle
      allocate(subat(nsub),sub_a(3,nsub),sub_b(3,nsub))
      k=0
      do i=1,nat
        if(frag_ids(i)/=f) cycle
        k=k+1
        subat(k)=at(i)
        sub_a(:,k)=xyz_a(:,i)
        sub_b(:,k)=xyz_b(:,i)
      end do
      ntopo=nsub*(nsub+1)/2
      allocate(topo_a(ntopo),topo_b(ntopo))
      call quicktopo(nsub,subat,sub_a,ntopo,topo_a)
      call quicktopo(nsub,subat,sub_b,ntopo,topo_b)
      if(any(topo_a/=topo_b)) then
        deallocate(topo_b,topo_a,sub_b,sub_a,subat,frag_ids)
        return
      end if
      deallocate(topo_b,topo_a,sub_b,sub_a,subat)
    end do
    equal=.true.
    deallocate(frag_ids)
  end function fragment_topology_equal

  subroutine get_consistent_gff_fragments(env,fragments,found)
    type(systemdata),intent(in) :: env
    character(len=:),allocatable,intent(out) :: fragments(:)
    logical,intent(out) :: found
    integer :: i,j,first,width
    found=.false.; first=0
    do i=1,env%calc%ncalculations
      if(.not.allocated(env%calc%calcs(i)%gff_fragments)) cycle
      if(first==0) then
        width=0
        do j=1,size(env%calc%calcs(i)%gff_fragments)
          if(.not.allocated(env%calc%calcs(i)%gff_fragments(j)%value)) then
            error stop '**ERROR** unallocated GFN-FF fragment string'
          end if
          width=max(width,len(env%calc%calcs(i)%gff_fragments(j)%value))
        end do
        if(width<1) error stop '**ERROR** empty GFN-FF fragment string array'
        allocate(character(len=width) :: fragments(size(env%calc%calcs(i)%gff_fragments)))
        do j=1,size(fragments)
          fragments(j)=env%calc%calcs(i)%gff_fragments(j)%value
        end do
        first=i; found=.true.
      else
        if(size(fragments)/=size(env%calc%calcs(i)%gff_fragments)) then
          error stop '**ERROR** inconsistent GFN-FF fragments across calculation levels'
        end if
        do j=1,size(fragments)
          if(.not.allocated(env%calc%calcs(i)%gff_fragments(j)%value)) then
            error stop '**ERROR** unallocated GFN-FF fragment string'
          end if
          if(trim(fragments(j))/=trim(env%calc%calcs(i)%gff_fragments(j)%value)) then
            error stop '**ERROR** inconsistent GFN-FF fragments across calculation levels'
          end if
        end do
      end if
    end do
  end subroutine get_consistent_gff_fragments

  subroutine remove_external_rerank_checkpoint()
    integer :: ich,io
    logical :: ex
    inquire(file=checkpoint_file,exist=ex)
    if(.not.ex) return
    open(newunit=ich,file=checkpoint_file,status='old',iostat=io)
    if(io==0) close(ich,status='delete')
  end subroutine remove_external_rerank_checkpoint
end module external_rerank_restart
