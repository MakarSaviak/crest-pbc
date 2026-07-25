!================================================================================!
! Explicit staged restart after external ranking of iteration-1 source conformers.
!================================================================================!
module external_rerank_restart
  use crest_parameters,only:wp,stdout,bohr
  use crest_data,only:systemdata
  use strucrd,only:coord,rdensembleparam,rdensemble,get_atlist
  implicit none
  private

  character(len=*),parameter :: checkpoint_file='crest.external_rerank'

  type,public :: external_rerank_state
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
    integer :: ich,i
    logical :: ex

    if (.not.allocated(env%input_settings_file)) then
      error stop '**ERROR** external reranking requires an original TOML settings file'
    end if
    inquire(file=trim(archive_file),exist=ex)
    if (.not.ex) error stop '**ERROR** iteration archive missing before external checkpoint'
    if (env%nmetadyn<1 .or. .not.allocated(env%metadfac) .or. .not.allocated(env%metadexp)) then
      error stop '**ERROR** no MTD bias state available for external checkpoint'
    end if
    if (size(env%metadfac)<env%nmetadyn .or. size(env%metadexp)<env%nmetadyn) then
      error stop '**ERROR** inconsistent MTD bias arrays at external checkpoint'
    end if

    open(newunit=ich,file=checkpoint_file,status='replace',action='write')
    write(ich,'(a)') 'CREST_EXTERNAL_RERANK_V2'
    write(ich,'(i0)') completed_mtd
    write(ich,'(i0)') env%Maxrestart
    write(ich,'(i0)') env%nmetadyn
    write(ich,'(es26.17)') env%elowest
    write(ich,'(es26.17)') env%eprivious
    write(ich,'(a)') trim(env%input_settings_file)
    write(ich,'(a)') trim(env%constraints)
    write(ich,'(a)') trim(archive_file)
    do i=1,env%nmetadyn
      write(ich,'(2es26.17)') env%metadfac(i),env%metadexp(i)
    end do
    close(ich)
  end subroutine write_external_rerank_checkpoint

  subroutine read_external_rerank_checkpoint(state)
    type(external_rerank_state),intent(out) :: state
    integer :: ich,io,i
    character(len=64) :: magic
    logical :: ex

    inquire(file=checkpoint_file,exist=ex)
    if (.not.ex) error stop '**ERROR** external-rerank checkpoint not found'
    open(newunit=ich,file=checkpoint_file,status='old',action='read')
    read(ich,'(a)',iostat=io) magic
    if (io/=0 .or. trim(magic)/='CREST_EXTERNAL_RERANK_V2') then
      error stop '**ERROR** invalid external-rerank checkpoint format'
    end if
    read(ich,*,iostat=io) state%completed_mtd
    if(io==0) read(ich,*,iostat=io) state%target_mtd
    if(io==0) read(ich,*,iostat=io) state%nmetadyn
    if(io==0) read(ich,*,iostat=io) state%elowest
    if(io==0) read(ich,*,iostat=io) state%eprivious
    if(io==0) read(ich,'(a)',iostat=io) state%settings_file
    if(io==0) read(ich,'(a)',iostat=io) state%constraints_file
    if(io==0) read(ich,'(a)',iostat=io) state%archive_file
    if(io/=0) then
      close(ich)
      error stop '**ERROR** truncated external-rerank checkpoint header'
    end if
    if(state%completed_mtd<1 .or. state%target_mtd<=state%completed_mtd) then
      close(ich)
      error stop '**ERROR** invalid external-rerank iteration state'
    end if
    if(state%nmetadyn<1) then
      close(ich)
      error stop '**ERROR** invalid external-rerank bias count'
    end if
    allocate(state%metadfac(state%nmetadyn),state%metadexp(state%nmetadyn))
    do i=1,state%nmetadyn
      read(ich,*,iostat=io) state%metadfac(i),state%metadexp(i)
      if(io/=0) exit
    end do
    close(ich)
    if(io/=0) error stop '**ERROR** truncated external-rerank bias state'
    if(len_trim(state%settings_file)==0 .or. len_trim(state%archive_file)==0) then
      error stop '**ERROR** external-rerank checkpoint lacks required paths'
    end if
  end subroutine read_external_rerank_checkpoint

  subroutine restore_external_rerank_state(env,state)
    type(systemdata),intent(inout) :: env
    type(external_rerank_state),intent(in) :: state
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
    integer :: i,j,first
    found=.false.; first=0
    do i=1,env%calc%ncalculations
      if(.not.allocated(env%calc%calcs(i)%gff_fragments)) cycle
      if(first==0) then
        fragments=env%calc%calcs(i)%gff_fragments
        first=i; found=.true.
      else
        if(size(fragments)/=size(env%calc%calcs(i)%gff_fragments)) then
          error stop '**ERROR** inconsistent GFN-FF fragments across calculation levels'
        end if
        do j=1,size(fragments)
          if(trim(fragments(j))/=trim(env%calc%calcs(i)%gff_fragments(j))) then
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
