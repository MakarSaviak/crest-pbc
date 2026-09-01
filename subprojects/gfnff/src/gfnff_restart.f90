!================================================================================!
! This file is part of gfnff.
!
! Copyright (C) 2023 Philipp Pracht
!
! gfnff is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! gfnff is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with gfnff. If not, see <https://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------!
!> The original (unmodified) source code can be found under the GNU LGPL 3.0 license
!> Copyright (C) 2019-2020 Sebastian Ehlert, Sebastian Spicher, Stefan Grimme
!> at https://github.com/grimme-lab/xtb
!================================================================================!
module gfnff_restart
  use iso_fortran_env,only:wp => real64,stdout => output_unit,iostat_end
  use gfnff_data_types
  use gfnff_param
  implicit none
  private
  public :: write_restart_gff,read_restart_gff

  integer,parameter :: i8 = selected_int_kind(18)
  integer(i8),parameter :: restart_extension_magic = int(z'47464653544F5253',i8)
  integer(i8),parameter :: restart_hyb_extension_magic = int(z'4746464859423031',i8)
  integer,parameter :: special_torsion_rows = 6

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine write_restart_gff(fname,nat,version,topo,iostat_out)
    implicit none
    type(TGFFTopology),intent(in) :: topo
    character(len=*),intent(in) :: fname
    integer,intent(in)  :: nat
    integer,intent(in)  :: version
    integer,intent(out),optional :: iostat_out
    integer :: ich ! file handle
    integer :: err
    integer :: close_err
    integer(i8) :: sTorsl_rows,sTorsl_columns,hyb_length
    logical :: have_sTorsl,have_hyb

    err = 0
    open (file=fname,newunit=ich,status='replace',action='write',form='unformatted',iostat=err)
    if (err /= 0) then
      if (present(iostat_out)) iostat_out = err
      return
    end if

!>--- Dimensions
    write (ich,iostat=err) int(version,i8),int(nat,i8)
    if (err == 0) write (ich,iostat=err) topo%nbond,topo%nangl,topo%ntors,   &
             &  topo%nathbH,topo%nathbAB,topo%natxbAB,topo%nbatm, &
             &  topo%nfrag,topo%nsystem,topo%maxsystem
    if (err == 0) write (ich,iostat=err) topo%nbond_blist,topo%nbond_vbond,topo%nangl_alloc, &
             &  topo%ntors_alloc,topo%bond_hb_nr,topo%b_max
!>--- Arrays Integers
    if (err == 0) write (ich,iostat=err) topo%nb,topo%bpair,topo%blist,topo%alist, &
             &  topo%tlist,topo%b3list,topo%fraglist,topo%hbatHl,topo%hbatABl, &
             &  topo%xbatABl,topo%ispinsyst,topo%nspinsyst,topo%bond_hb_AH, &
             &  topo%bond_hb_B,topo%bond_hb_Bn,topo%nr_hb
!>--- Arrays Reals
    if (err == 0) write (ich,iostat=err) topo%vbond,topo%vangl,topo%vtors,topo%chieeq, &
             &  topo%gameeq,topo%alpeeq,topo%alphanb,topo%qa,  &
             &  topo%xyze0,topo%zetac6, &
             &  topo%qfrag,topo%hbbas,topo%hbaci
!>--- Trailing extensions: static special C#C torsions and hybridization were
!>    not part of the legacy restart records.  Keep the legacy parameterization
!>    version intact, and retain the original extension before appending hyb.
    if (err == 0) write (ich,iostat=err) restart_extension_magic
    have_sTorsl = allocated(topo%sTorsl)
    if (have_sTorsl) then
      sTorsl_rows = int(size(topo%sTorsl,1),i8)
      sTorsl_columns = int(size(topo%sTorsl,2),i8)
    else
      sTorsl_rows = 0_i8
      sTorsl_columns = 0_i8
    end if
    if (err == 0) write (ich,iostat=err) have_sTorsl,sTorsl_rows,sTorsl_columns
    if (err == 0 .and. have_sTorsl) write (ich,iostat=err) topo%sTorsl
    if (err == 0) write (ich,iostat=err) restart_hyb_extension_magic
    have_hyb = allocated(topo%hyb)
    if (have_hyb) then
      hyb_length = int(size(topo%hyb),i8)
    else
      hyb_length = 0_i8
    end if
    if (err == 0) write (ich,iostat=err) have_hyb,hyb_length
    if (err == 0 .and. have_hyb) write (ich,iostat=err) topo%hyb
    close (ich,iostat=close_err)
    if (err == 0 .and. close_err /= 0) err = close_err
    if (present(iostat_out)) iostat_out = err
  end subroutine write_restart_gff
!========================================================================================!

  subroutine read_restart_gff(fname,n,version,success,verbose,topo)
    implicit none
    character(len=*),parameter :: source = 'restart_read_restart_gff'
    type(TGFFTopology),intent(inout) :: topo
    character(len=*),intent(in) :: fname
    integer,intent(in)  :: n
    integer,intent(in)  :: version
    logical,intent(out) :: success
    logical,intent(in)  :: verbose

    integer(i8) :: iver8,nat8
    integer(i8) :: extension_magic,sTorsl_rows,sTorsl_columns
    integer(i8) :: hyb_extension_magic,hyb_length

    integer :: ich ! file handle
    integer :: err
    integer :: alloc_err
    logical :: have_sTorsl,have_hyb

    success = .false.
    open (file=fname,newunit=ich,status='old',action='read', &
         & form='unformatted',iostat=err)
    if (err == 0) then
!>--- read the first byte, which identify the calculation specs
      read (ich,iostat=err) iver8,nat8
      if (err .eq. 0) then
        if (iver8 .ne. int(version,i8).and.verbose) &
           &  write (stdout,'("Version number missmatch in restart file.",a)') source
        if (nat8 .ne. n.and.verbose) then
          write (stdout,'("Atom number missmatch in restart file.",a)') source
          success = .false.
          close (ich)
          return
        else if (iver8 .eq. int(version)) then
          success = .true.
          read (ich) topo%nbond,topo%nangl,topo%ntors, &
                  &  topo%nathbH,topo%nathbAB,topo%natxbAB,topo%nbatm, &
                  &  topo%nfrag,topo%nsystem,topo%maxsystem
          read (ich) topo%nbond_blist,topo%nbond_vbond,topo%nangl_alloc, &
                  &  topo%ntors_alloc,topo%bond_hb_nr,topo%b_max
!>--- allocate some memory now
          call gfnff_param_alloc(topo,n)
          if (.not.allocated(topo%ispinsyst)) allocate (topo%ispinsyst(n,topo%maxsystem),source=0)
          if (.not.allocated(topo%nspinsyst)) allocate (topo%nspinsyst(topo%maxsystem),source=0)
          read (ich) topo%nb,topo%bpair,topo%blist,topo%alist, &
             & topo%tlist,topo%b3list,topo%fraglist,topo%hbatHl,topo%hbatABl, &
             & topo%xbatABl,topo%ispinsyst,topo%nspinsyst,topo%bond_hb_AH, &
             & topo%bond_hb_B,topo%bond_hb_Bn,topo%nr_hb
          read (ich) topo%vbond,topo%vangl,topo%vtors,topo%chieeq, &
             & topo%gameeq,topo%alpeeq,topo%alphanb,topo%qa, &
             & topo%xyze0,topo%zetac6,&
             & topo%qfrag,topo%hbbas,topo%hbaci
          if (allocated(topo%sTorsl)) deallocate(topo%sTorsl)
          if (allocated(topo%hyb)) deallocate(topo%hyb)
          ! A legacy file ends exactly here.  A new file carries the trailing
          ! extension, which must either be complete and valid or fail closed.
          read (ich,iostat=err) extension_magic
          if (err == iostat_end) then
            err = 0
          else if (err /= 0) then
            if (verbose) write (stdout,'("Malformed restart extension.",a)') source
            success = .false.
          else if (extension_magic /= restart_extension_magic) then
            if (verbose) write (stdout,'("Unknown restart extension.",a)') source
            success = .false.
          else
            read (ich,iostat=err) have_sTorsl,sTorsl_rows,sTorsl_columns
            if (err /= 0) then
              if (verbose) write (stdout,'("Malformed restart extension.",a)') source
              success = .false.
            else if (have_sTorsl) then
              if (sTorsl_rows /= int(special_torsion_rows,i8) .or. &
              & sTorsl_columns < 1_i8 .or. &
              & sTorsl_columns > int(huge(0),i8)) then
                if (verbose) write (stdout,'("Invalid special torsion restart dimensions.",a)') source
                success = .false.
              else
                allocate(topo%sTorsl(special_torsion_rows,int(sTorsl_columns)),stat=alloc_err)
                if (alloc_err /= 0) then
                  if (verbose) write (stdout,'("Cannot allocate special torsion restart state.",a)') source
                  success = .false.
                else
                  read (ich,iostat=err) topo%sTorsl
                  if (err /= 0) then
                    deallocate(topo%sTorsl)
                    if (verbose) write (stdout,'("Malformed special torsion restart state.",a)') source
                    success = .false.
                  end if
                end if
              end if
            else if (sTorsl_rows /= 0_i8 .or. sTorsl_columns /= 0_i8) then
              if (verbose) write (stdout,'("Invalid empty special torsion restart state.",a)') source
              success = .false.
            end if
            if (success) then
              ! Existing first-generation extensions end here.  A second,
              ! separately marked trailing record preserves hybridization.
              read (ich,iostat=err) hyb_extension_magic
              if (err == iostat_end) then
                err = 0
              else if (err /= 0) then
                if (verbose) write (stdout,'("Malformed hybridization restart extension.",a)') source
                success = .false.
              else if (hyb_extension_magic /= restart_hyb_extension_magic) then
                if (verbose) write (stdout,'("Unknown hybridization restart extension.",a)') source
                success = .false.
              else
                read (ich,iostat=err) have_hyb,hyb_length
                if (err /= 0) then
                  if (verbose) write (stdout,'("Malformed hybridization restart extension.",a)') source
                  success = .false.
                else if (have_hyb) then
                  if (hyb_length /= int(n,i8)) then
                    if (verbose) write (stdout,'("Invalid hybridization restart dimensions.",a)') source
                    success = .false.
                  else
                    allocate(topo%hyb(n),stat=alloc_err)
                    if (alloc_err /= 0) then
                      if (verbose) write (stdout,'("Cannot allocate hybridization restart state.",a)') source
                      success = .false.
                    else
                      read (ich,iostat=err) topo%hyb
                      if (err /= 0) then
                        deallocate(topo%hyb)
                        if (verbose) write (stdout,'("Malformed hybridization restart state.",a)') source
                        success = .false.
                      end if
                    end if
                  end if
                else if (hyb_length /= 0_i8) then
                  if (verbose) write (stdout,'("Invalid empty hybridization restart state.",a)') source
                  success = .false.
                end if
              end if
            end if
          end if
        else
          if (verbose) &
            write (stdout,'("Dimension missmatch in restart file.",a)') source
          success = .false.
        end if
      else
        if (verbose) &
          write (stdout,'("Dimension missmatch in restart file.",a)') source
        success = .false.
      end if
      close (ich)
    end if

  end subroutine read_restart_gff

!========================================================================================!
end module gfnff_restart
