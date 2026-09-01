!================================================================================!
! In-memory carrier and fail-closed reader for canonical process-MTD trajectories.
!================================================================================!
module crest_poststage_ensemble
  use iso_fortran_env, only : wp => real64, int64, iostat_end
  use iso_c_binding, only : c_char, c_double, c_int, c_int64_t, c_loc, &
  & c_null_char, c_ptr
  use ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use crest_data, only : status_normal, status_ioerr, status_args, status_input
  use strucrd, only : e2i, grepenergy, i2e
  implicit none
  private

  integer, parameter, public :: poststage_comment_length = 128

  type, public :: poststage_ensemble
    ! Coordinates are in Angstrom until the optimization caller explicitly
    ! performs the same conversion as the legacy formatted-reader path.
    integer :: nat = 0
    integer :: nall = 0
    integer, allocatable :: at(:)
    real(wp), allocatable :: xyz(:,:,:)
    real(wp), allocatable :: eread(:)
    character(len=poststage_comment_length), allocatable :: comments(:)
  contains
    procedure :: clear => clear_poststage_ensemble
    procedure :: valid => valid_poststage_ensemble
  end type poststage_ensemble

  public :: move_poststage_ensemble
  public :: round_fixed_decimal_10
  public :: write_canonical_ensemble_fast
  public :: read_canonical_trajectory_slice
  public :: read_canonical_trajectory_slice_fast
  public :: read_canonical_trajectory_slice_legacy

  interface
    function c_read_canonical_xyz_fast(path,expected_at,nat,expected_frames,xyz, &
    & energies,comments,comment_len,file_bytes,read_seconds,decode_seconds, &
    & error_frame,error_atom,message,message_len) &
    & bind(C,name='crest_read_canonical_xyz_fast') result(rc)
      import :: c_ptr, c_int, c_int64_t, c_double
      type(c_ptr), value :: path
      type(c_ptr), value :: expected_at
      integer(c_int), value :: nat
      integer(c_int), value :: expected_frames
      type(c_ptr), value :: xyz
      type(c_ptr), value :: energies
      type(c_ptr), value :: comments
      integer(c_int), value :: comment_len
      integer(c_int64_t), intent(out) :: file_bytes
      real(c_double), intent(out) :: read_seconds
      real(c_double), intent(out) :: decode_seconds
      integer(c_int), intent(out) :: error_frame
      integer(c_int), intent(out) :: error_atom
      type(c_ptr), value :: message
      integer(c_int), value :: message_len
      integer(c_int) :: rc
    end function c_read_canonical_xyz_fast

    function c_write_canonical_xyz_fast(path,atomic_numbers,nat,nframes,xyz, &
    & comments,comment_len,file_bytes,write_seconds,error_frame,error_atom, &
    & message,message_len) &
    & bind(C,name='crest_write_canonical_xyz_fast') result(rc)
      import :: c_ptr, c_int, c_int64_t, c_double
      type(c_ptr), value :: path
      type(c_ptr), value :: atomic_numbers
      integer(c_int), value :: nat
      integer(c_int), value :: nframes
      type(c_ptr), value :: xyz
      type(c_ptr), value :: comments
      integer(c_int), value :: comment_len
      integer(c_int64_t), intent(out) :: file_bytes
      real(c_double), intent(out) :: write_seconds
      integer(c_int), intent(out) :: error_frame
      integer(c_int), intent(out) :: error_atom
      type(c_ptr), value :: message
      integer(c_int), value :: message_len
      integer(c_int) :: rc
    end function c_write_canonical_xyz_fast
  end interface

contains

  subroutine round_fixed_decimal_10(value,result,status,used_fallback)
    ! Reproduce the numeric value obtained from an F20.10 formatted
    ! write/read without paying the formatted-I/O cost for ordinary values.
    ! Values close enough to a decimal half-way boundary fall back to the
    ! runtime formatter/parser, preserving its tie-breaking semantics.
    real(wp), intent(in) :: value
    real(wp), intent(out) :: result
    integer, intent(out) :: status
    logical, intent(out) :: used_fallback

    real(wp), parameter :: scale = 1.0e10_wp
    real(wp) :: scaled,nearest,distance_to_half,tolerance
    character(len=32) :: buffer
    integer :: io

    status = 0
    used_fallback = .false.
    scaled = value*scale
    nearest = anint(scaled)
    distance_to_half = 0.5_wp-abs(scaled-nearest)
    tolerance = max(8.0_wp*spacing(scaled),1.0e-10_wp)

    if (.not.ieee_is_finite(value) .or. abs(value) >= 1.0e8_wp .or. &
    & distance_to_half <= tolerance) then
      used_fallback = .true.
      write(buffer,'(f20.10)',iostat=io) value
      if (io /= 0) then
        status = io
        result = value
        return
      end if
      read(buffer,*,iostat=io) result
      if (io /= 0) then
        status = io
        result = value
      end if
    else
      result = nearest/scale
    end if
  end subroutine round_fixed_decimal_10

  subroutine write_canonical_ensemble_fast(path,at,xyz,comments,file_bytes, &
  & write_seconds,status,message)
    character(len=*), intent(in) :: path
    integer, intent(in) :: at(:)
    real(wp), contiguous, target, intent(in) :: xyz(:,:,:)
    character(len=poststage_comment_length), contiguous, intent(in) :: comments(:)
    integer(int64), intent(out) :: file_bytes
    real(wp), intent(out) :: write_seconds
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    character(kind=c_char), allocatable, target :: c_path(:),comment_bytes(:),c_message(:)
    integer(c_int), allocatable, target :: c_at(:)
    integer(c_int) :: rc,error_frame,error_atom
    integer(c_int64_t) :: c_file_bytes
    real(c_double) :: c_write_seconds
    integer :: i,frame,first,last,npath,nframes

    status = status_normal
    message = ''
    file_bytes = 0_int64
    write_seconds = 0.0_wp
    npath = len_trim(path)
    nframes = size(comments)
    if (npath < 1 .or. size(at) < 1 .or. nframes < 1) then
      status = status_args
      message = 'invalid fast canonical writer arguments'
      return
    end if
    if (size(xyz,1) /= 3 .or. size(xyz,2) /= size(at) .or. &
    & size(xyz,3) /= nframes) then
      status = status_args
      message = 'inconsistent fast canonical writer dimensions'
      return
    end if

    allocate(c_path(npath+1),c_at(size(at)), &
    & comment_bytes(poststage_comment_length*nframes),c_message(1024))
    do i = 1,npath
      c_path(i) = path(i:i)
    end do
    c_path(npath+1) = c_null_char
    c_at = int(at,c_int)
    do frame = 1,nframes
      first = (frame-1)*poststage_comment_length+1
      last = first+poststage_comment_length-1
      comment_bytes(first:last) = transfer(comments(frame),comment_bytes(first:last))
    end do
    c_message = c_null_char

    rc = c_write_canonical_xyz_fast(c_loc(c_path(1)),c_loc(c_at(1)), &
    & int(size(at),c_int),int(nframes,c_int),c_loc(xyz(1,1,1)), &
    & c_loc(comment_bytes(1)),int(poststage_comment_length,c_int), &
    & c_file_bytes,c_write_seconds,error_frame,error_atom,c_loc(c_message(1)), &
    & int(size(c_message),c_int))

    file_bytes = int(c_file_bytes,int64)
    write_seconds = real(c_write_seconds,wp)
    select case(int(rc))
    case(0)
      status = status_normal
    case(1)
      status = status_args
    case(2,4)
      status = status_ioerr
    case default
      status = status_input
    end select
    if (status /= status_normal) then
      do i = 1,min(len(message),size(c_message))
        if (c_message(i) == c_null_char) exit
        message(i:i) = c_message(i)
      end do
      if (len_trim(message) == 0) then
        write(message,'(a,i0,a,i0)') 'fast canonical writer failed at frame ', &
        & int(error_frame),', atom ',int(error_atom)
      end if
    end if
  end subroutine write_canonical_ensemble_fast

  subroutine clear_poststage_ensemble(self)
    class(poststage_ensemble), intent(inout) :: self

    if (allocated(self%at)) deallocate(self%at)
    if (allocated(self%xyz)) deallocate(self%xyz)
    if (allocated(self%eread)) deallocate(self%eread)
    if (allocated(self%comments)) deallocate(self%comments)
    self%nat = 0
    self%nall = 0
  end subroutine clear_poststage_ensemble

  logical function valid_poststage_ensemble(self) result(valid)
    class(poststage_ensemble), intent(in) :: self

    valid = .false.
    if (self%nat < 1 .or. self%nall < 1) return
    if (.not.allocated(self%at)) return
    if (.not.allocated(self%xyz)) return
    if (.not.allocated(self%eread)) return
    if (.not.allocated(self%comments)) return
    if (size(self%at) /= self%nat) return
    if (size(self%xyz,1) /= 3) return
    if (size(self%xyz,2) /= self%nat) return
    if (size(self%xyz,3) /= self%nall) return
    if (size(self%eread) /= self%nall) return
    if (size(self%comments) /= self%nall) return
    valid = .true.
  end function valid_poststage_ensemble

  subroutine move_poststage_ensemble(from,to)
    type(poststage_ensemble), intent(inout) :: from
    type(poststage_ensemble), intent(inout) :: to

    ! FROM and TO must be distinct objects.  Moving each allocatable component
    ! avoids a multi-gigabyte intrinsic-assignment copy of the coordinate data.
    call to%clear()
    to%nat = from%nat
    to%nall = from%nall
    if (allocated(from%at)) call move_alloc(from%at,to%at)
    if (allocated(from%xyz)) call move_alloc(from%xyz,to%xyz)
    if (allocated(from%eread)) call move_alloc(from%eread,to%eread)
    if (allocated(from%comments)) call move_alloc(from%comments,to%comments)
    from%nat = 0
    from%nall = 0
  end subroutine move_poststage_ensemble


  subroutine read_canonical_trajectory_slice(path,expected_at,expected_frames, &
  &                                           xyz,eread,comments,file_bytes, &
  &                                           status,message,read_seconds,decode_seconds)
    character(len=*), intent(in) :: path
    integer, intent(in) :: expected_at(:)
    integer, intent(in) :: expected_frames
    real(wp), contiguous, target, intent(out) :: xyz(:,:,:)
    real(wp), contiguous, target, intent(out) :: eread(:)
    character(len=poststage_comment_length), contiguous, intent(out) :: comments(:)
    integer(int64), intent(out) :: file_bytes
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(wp), intent(out), optional :: read_seconds,decode_seconds

    call read_canonical_trajectory_slice_fast(path,expected_at,expected_frames, &
    & xyz,eread,comments,file_bytes,status,message,read_seconds,decode_seconds)
  end subroutine read_canonical_trajectory_slice

  subroutine read_canonical_trajectory_slice_fast(path,expected_at,expected_frames, &
  &                                                xyz,eread,comments,file_bytes, &
  &                                                status,message,read_seconds, &
  &                                                decode_seconds)
    ! Parse a complete process-MTD trajectory without invoking the formatted
    ! Fortran runtime for every record.  One private byte buffer is allocated
    ! per active parser worker and coordinates are written directly into the
    ! caller-owned final array slice.
    character(len=*), intent(in) :: path
    integer, intent(in) :: expected_at(:)
    integer, intent(in) :: expected_frames
    real(wp), contiguous, target, intent(out) :: xyz(:,:,:)
    real(wp), contiguous, target, intent(out) :: eread(:)
    character(len=poststage_comment_length), contiguous, intent(out) :: comments(:)
    integer(int64), intent(out) :: file_bytes
    integer, intent(out) :: status
    character(len=*), intent(out) :: message
    real(wp), intent(out), optional :: read_seconds,decode_seconds

    character(kind=c_char), allocatable, target :: c_path(:),comment_bytes(:)
    character(kind=c_char), allocatable, target :: c_message(:)
    integer(c_int), allocatable, target :: expected_at_c(:)
    integer(c_int) :: rc,error_frame,error_atom
    integer(c_int64_t) :: file_bytes_c
    real(c_double) :: read_seconds_c,decode_seconds_c
    integer :: i,frame,first,last,path_length,message_length

    status = status_normal
    message = ''
    file_bytes = 0_int64
    if (present(read_seconds)) read_seconds = 0.0_wp
    if (present(decode_seconds)) decode_seconds = 0.0_wp

    if (len_trim(path) < 1) then
      status = status_args
      message = 'canonical trajectory path is empty'
      return
    end if
    if (size(expected_at) < 1) then
      status = status_args
      message = 'canonical trajectory expected atom list is empty'
      return
    end if
    if (expected_frames < 1) then
      status = status_args
      message = 'canonical trajectory expected frame count is not positive'
      return
    end if
    if (size(xyz,1) /= 3 .or. size(xyz,2) /= size(expected_at) .or. &
    &   size(xyz,3) /= expected_frames .or. size(eread) /= expected_frames .or. &
    &   size(comments) /= expected_frames) then
      status = status_args
      message = 'canonical trajectory destination slice has inconsistent dimensions'
      return
    end if

    path_length = len_trim(path)
    message_length = max(256,len(message)+1)
    allocate(c_path(path_length+1),expected_at_c(size(expected_at)), &
    & comment_bytes(poststage_comment_length*expected_frames), &
    & c_message(message_length),stat=i)
    if (i /= 0) then
      status = status_ioerr
      message = 'fast canonical trajectory parser workspace allocation failed'
      return
    end if
    do i = 1,path_length
      c_path(i) = path(i:i)
    end do
    c_path(path_length+1) = c_null_char
    expected_at_c = int(expected_at,c_int)
    c_message = c_null_char

    rc = c_read_canonical_xyz_fast(c_loc(c_path(1)),c_loc(expected_at_c(1)), &
    & int(size(expected_at),c_int),int(expected_frames,c_int), &
    & c_loc(xyz(1,1,1)),c_loc(eread(1)),c_loc(comment_bytes(1)), &
    & int(poststage_comment_length,c_int),file_bytes_c,read_seconds_c, &
    & decode_seconds_c,error_frame,error_atom,c_loc(c_message(1)), &
    & int(message_length,c_int))

    file_bytes = int(file_bytes_c,int64)
    if (present(read_seconds)) read_seconds = real(read_seconds_c,wp)
    if (present(decode_seconds)) decode_seconds = real(decode_seconds_c,wp)

    select case (int(rc))
    case (0)
      status = status_normal
    case (1)
      status = status_args
    case (2,4)
      status = status_ioerr
    case default
      status = status_input
    end select
    if (status == status_normal) then
      do frame = 1,expected_frames
        first = (frame-1)*poststage_comment_length+1
        last = first+poststage_comment_length-1
        comments(frame) = transfer(comment_bytes(first:last),comments(frame))
      end do
    else
      message = ''
      do i = 1,min(len(message),size(c_message))
        if (c_message(i) == c_null_char) exit
        message(i:i) = c_message(i)
      end do
      if (len_trim(message) == 0) then
        write(message,'(a,i0,a,i0)') 'fast trajectory parser failed at frame ', &
        & int(error_frame),', atom ',int(error_atom)
      end if
    end if
  end subroutine read_canonical_trajectory_slice_fast

  subroutine read_canonical_trajectory_slice_legacy(path,expected_at,expected_frames, &
  &                                                  xyz,eread,comments,file_bytes, &
  &                                                  status,message)
    ! Read one already-produced canonical process-MTD trajectory into a caller-
    ! owned slice.  The routine performs no allocation and has no ERROR STOP;
    ! callers can therefore invoke it concurrently on disjoint array sections.
    character(len=*), intent(in) :: path
    integer, intent(in) :: expected_at(:)
    integer, intent(in) :: expected_frames
    real(wp), contiguous, intent(out) :: xyz(:,:,:)
    real(wp), contiguous, intent(out) :: eread(:)
    character(len=poststage_comment_length), contiguous, intent(out) :: comments(:)
    integer(int64), intent(out) :: file_bytes
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    integer :: unit,raw_status,close_status,comment_status,lexical_status
    integer :: frame,atom,natoms,atomic_number
    real(wp) :: coordinates(3),unassigned_coordinate,comment_energy
    logical :: exists,opened
    character(len=6) :: symbol
    character(len=32) :: energy_label,equals,trailing_token
    character(len=512) :: line
    character(len=1024) :: local_message,close_message,detail

    status = status_normal
    message = ''
    file_bytes = 0_int64
    unit = -1
    opened = .false.
    unassigned_coordinate = ieee_value(0.0_wp,ieee_quiet_nan)

    if (len_trim(path) < 1) then
      status = status_args
      message = 'canonical trajectory path is empty'
      return
    end if
    if (size(expected_at) < 1) then
      status = status_args
      message = 'canonical trajectory expected atom list is empty'
      return
    end if
    if (any(expected_at < 1) .or. any(expected_at > 118)) then
      status = status_args
      message = 'canonical trajectory expected atom list is outside 1:118'
      return
    end if
    if (expected_frames < 1) then
      status = status_args
      message = 'canonical trajectory expected frame count is not positive'
      return
    end if
    if (size(xyz,1) /= 3 .or. size(xyz,2) /= size(expected_at) .or. &
    &   size(xyz,3) /= expected_frames .or. size(eread) /= expected_frames .or. &
    &   size(comments) /= expected_frames) then
      status = status_args
      message = 'canonical trajectory destination slice has inconsistent dimensions'
      return
    end if

    local_message = ''
    inquire(file=trim(path),exist=exists,size=file_bytes,iostat=raw_status, &
    &       iomsg=local_message)
    if (raw_status /= 0) then
      status = status_ioerr
      message = 'cannot inspect canonical trajectory '//trim(path)//': '// &
      &         trim(local_message)
      return
    end if
    if (.not.exists) then
      status = status_ioerr
      message = 'canonical trajectory is missing: '//trim(path)
      return
    end if
    if (file_bytes < 1_int64) then
      status = status_input
      message = 'canonical trajectory is empty: '//trim(path)
      return
    end if

    local_message = ''
    open(newunit=unit,file=trim(path),status='old',action='read',form='formatted', &
    &    iostat=raw_status,iomsg=local_message)
    if (raw_status /= 0) then
      status = status_ioerr
      message = 'cannot open canonical trajectory '//trim(path)//': '// &
      &         trim(local_message)
      return
    end if
    opened = .true.

    do frame = 1,expected_frames
      local_message = ''
      read(unit,*,iostat=raw_status,iomsg=local_message) natoms
      if (raw_status == iostat_end) then
        status = status_input
        write(detail,'(a,i0,a,i0)') 'canonical trajectory ended before frame ', &
        & frame,' of ',expected_frames
        message = trim(detail)
        goto 900
      else if (raw_status /= 0) then
        status = status_input
        write(detail,'(a,i0,2a)') 'malformed atom count at canonical frame ', &
        & frame,': ',trim(local_message)
        message = trim(detail)
        goto 900
      end if
      if (natoms /= size(expected_at)) then
        status = status_input
        write(detail,'(a,i0,a,i0,a,i0)') 'canonical frame ',frame, &
        & ' has atom count ',natoms,'; expected ',size(expected_at)
        message = trim(detail)
        goto 900
      end if

      local_message = ''
      read(unit,'(a)',iostat=raw_status,iomsg=local_message) line
      if (raw_status /= 0) then
        status = status_input
        write(detail,'(a,i0,2a)') 'cannot read comment at canonical frame ', &
        & frame,': ',trim(local_message)
        message = trim(detail)
        goto 900
      end if
      energy_label = ''
      equals = ''
      comment_energy = 0.0_wp
      read(line,*,iostat=comment_status) energy_label,equals,comment_energy
      if (comment_status /= 0 .or. trim(energy_label) /= 'Epot' .or. &
      &   trim(equals) /= '=' .or. .not.ieee_is_finite(comment_energy)) then
        status = status_input
        write(detail,'(a,i0)') 'noncanonical Epot comment at frame ',frame
        message = trim(detail)
        goto 900
      end if
      comments(frame) = trim(line)
      eread(frame) = grepenergy(line)
      if (.not.ieee_is_finite(eread(frame))) then
        status = status_input
        write(detail,'(a,i0)') 'nonfinite parsed energy at frame ',frame
        message = trim(detail)
        goto 900
      end if

      do atom = 1,natoms
        local_message = ''
        read(unit,'(a)',iostat=raw_status,iomsg=local_message) line
        if (raw_status /= 0) then
          status = status_input
          write(detail,'(a,i0,a,i0,2a)') 'cannot read canonical frame ',frame, &
          & ', atom ',atom,': ',trim(local_message)
          message = trim(detail)
          goto 900
        end if
        symbol = ''
        coordinates = unassigned_coordinate
        trailing_token = ''
        read(line,*,iostat=lexical_status) symbol,coordinates,trailing_token
        if (lexical_status /= iostat_end .or. &
        & trim(symbol) /= trim(i2e(expected_at(atom),'nc')) .or. &
        & .not.all(ieee_is_finite(coordinates))) then
          status = status_input
          write(detail,'(a,i0,a,i0)') 'noncanonical coordinate at frame ', &
          & frame,', atom ',atom
          message = trim(detail)
          goto 900
        end if
        atomic_number = e2i(symbol)
        if (atomic_number /= expected_at(atom) .or. &
        &   trim(symbol) /= trim(i2e(expected_at(atom),'nc'))) then
          status = status_input
          write(detail,'(a,i0,a,i0)') 'atom identity/order mismatch at frame ', &
          & frame,', atom ',atom
          message = trim(detail)
          goto 900
        end if
        xyz(:,atom,frame) = coordinates
      end do
    end do

    ! A successful fixed-count read must end here.  A complete additional frame,
    ! malformed trailing record, or partial record is a contract violation.
    local_message = ''
    read(unit,*,iostat=raw_status,iomsg=local_message) natoms
    if (raw_status == 0) then
      status = status_input
      message = 'canonical trajectory contains more frames than expected'
      goto 900
    else if (raw_status /= iostat_end) then
      status = status_input
      message = 'canonical trajectory has malformed trailing content: '// &
      &         trim(local_message)
      goto 900
    end if

900 continue
    if (opened) then
      close_message = ''
      close(unit,iostat=close_status,iomsg=close_message)
      opened = .false.
      if (close_status /= 0 .and. status == status_normal) then
        status = status_ioerr
        message = 'cannot close canonical trajectory '//trim(path)//': '// &
        &         trim(close_message)
      end if
    end if
  end subroutine read_canonical_trajectory_slice_legacy

end module crest_poststage_ensemble
