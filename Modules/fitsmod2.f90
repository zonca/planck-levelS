!-----------------------------------------------------------------------------
!
!  This file is part of the Planck simulation package
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version.

!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.

!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!-----------------------------------------------------------------------------

!  The Planck simulation package is being developed at the Max-Planck-Institut
!  fuer Astrophysik and financially supported by the Deutsches Zentrum fuer
!  Luft- und Raumfahrt (DLR).

!  Copyright (C) 2010-2013 Max-Planck-Society
!  \author Martin Reinecke

!! Low level FITS file access routines
module fitsmod2

use planck_config
use planck_types
use ls_misc_utils
use linebufmod
implicit none
private

! FITS data types
!! use if you want to create a binary column or image of 1 byte integer data
integer(i4b), parameter, public :: FITS_INT1=7
!! use if you want to create a binary column or image of 2 byte integer data
integer(i4b), parameter, public :: FITS_INT2=16
!! use if you want to create a binary column or image of 4 byte integer data
integer(i4b), parameter, public :: FITS_INT4=32
!! use if you want to create a binary column or image of 8 byte integer data
integer(i4b), parameter, public :: FITS_INT8=64
!! use if you want to create a binary column or image of 4 byte real data
integer(i4b), parameter, public :: FITS_REAL4=-32
!! use if you want to create a binary column or image of 8 byte real data
integer(i4b), parameter, public :: FITS_REAL8=-64
!! use if you want to create a binary column of character data
integer(i4b), parameter, public :: FITS_CHAR=8

! FITS file access types
!! symbolic constant specifying readonly file access
integer(i4b), parameter, public :: FITS_READONLY=0
!! symbolic constant specifying read/write file access
integer(i4b), parameter, public :: FITS_READWRITE=1

! FITS HDU types
!! symbolic constant specifying an image HDU
integer(i4b), parameter, public :: FITS_IMAGE=0
!! symbolic constant specifying an ASCII table HDU
integer(i4b), parameter, public :: FITS_ASCTAB=1
!! symbolic constant specifying a binary table HDU
integer(i4b), parameter, public :: FITS_BINTAB=2

! constant for yet undefined values
integer(i4b), parameter, public :: FITS_INVALID=-123456

!! data type containing information about a single column in a FITS table
type fitscolumn
  !! column name
  character(len=80) :: name=''
  !! physical unit of the column
  character(len=80) :: unit=''
  !! repetition count for numerical data types,
  !! string length for character columns
  integer(i4b) :: repcount=1
  !! data type of the column
  integer(i4b) :: type=FITS_INVALID
end type

public fitshandle, fitscolumn, fits_open, fits_close, fits_insert_bintab, &
  fits_create, fits_delete, fits_goto_hdu, fits_delete_chdu, fits_add_key, &
  fits_get_key, fits_write_column, fits_read_column, fits_add_comment, &
  fits_write_image, fits_read_image, fits_insert_image, fits_insert_asctab, &
  fits_read_header, fits_write_header, fits_key_present, fits_fast_nelms, &
  fits_update_key, fits_delete_key, fits_num_hdus, fitstype2type, type2fitstype

interface fits_add_key
  module procedure fits_add_key_int1, fits_add_key_int2, fits_add_key_int4, &
                   fits_add_key_int8, fits_add_key_real4, fits_add_key_real8, &
                   fits_add_key_string, fits_add_key_bool
end interface

interface fits_update_key
  module procedure fits_upd_key_int1, fits_upd_key_int2, fits_upd_key_int4, &
                   fits_upd_key_int8, fits_upd_key_real4, fits_upd_key_real8, &
                   fits_upd_key_string, fits_upd_key_bool
end interface

interface fits_get_key
  module procedure fits_get_key_int1, fits_get_key_int2, fits_get_key_int4, &
                   fits_get_key_int8, fits_get_key_real4, fits_get_key_real8, &
                   fits_get_key_string, fits_get_key_bool
end interface

interface fits_write_column
  module procedure fits_write_column_real4, fits_write_column_real8, &
                   fits_write_column_int1, fits_write_column_int2, &
                   fits_write_column_int4, fits_write_column_int8, &
                   fits_write_column_string
end interface

interface fits_read_column
  module procedure fits_read_column_real4, fits_read_column_real8, &
                   fits_read_column_int1, fits_read_column_int2, &
                   fits_read_column_int4, fits_read_column_int8, &
                   fits_read_column_string
end interface

interface fits_write_image
  module procedure fits_write_image_real4, fits_write_image_real8, &
                   fits_write_image_int4
end interface

interface fits_read_image
  module procedure fits_read_image_real4, fits_read_image_real8, &
                   fits_read_image_int4
end interface

!! data type containing information about the currently connected HDU
type fitshandle
  !! logical unit number, internal use only
  integer(i4b) :: lun=-1
  !! READ-ONLY: type of current HDU
  integer(i4b) :: hdutype=FITS_INVALID

  !! READ-ONLY: image dimensions (size of the array = number of dimensions)
  !! (only defined if current HDU is an image)
  integer(i4b), dimension(:), pointer :: axes=>NULL()
  !! READ-ONLY: data type of the image
  !! (only defined if current HDU is an image)
  integer(i4b) :: datatype=FITS_INVALID

  !! READ-ONLY: information about the columns
  !! (size of the array = number of columns)
  !! (only defined if current HDU is a table)
  type(fitscolumn), dimension(:), pointer :: columns=>NULL()
  !! READ-ONLY: current number of rows in the table
  !! (only defined if current HDU is a table)
  integer(i8b) :: nrows=0
end type

! internal parameter describing the maximum length of a comment
integer, parameter :: comlen=2000

contains

function type2char (type)
  integer(i4b), intent(in) :: type
  character :: type2char

  select case (type)
    case (FITS_REAL4)
      type2char='E'
    case (FITS_REAL8)
      type2char='D'
    case (FITS_INT1)
      type2char='B'
    case (FITS_INT2)
      type2char='I'
    case (FITS_INT4)
      type2char='J'
    case (FITS_INT8)
      type2char='K'
    case (FITS_CHAR)
      type2char='A'
    case default
      call exit_with_status (1, 'FITS: wrong data type')
  end select
end function

subroutine type2asciiform (type,form)
  integer(i4b), intent(in) :: type
  character(len=*),intent(out) :: form

  select case (type)
    case (FITS_REAL4)
      form='E14.7'
    case (FITS_REAL8)
      form='D23.15'
    case (FITS_INT1)
      form='I4'
    case (FITS_INT2)
      form='I6'
    case (FITS_INT4)
      form='I11'
    case (FITS_INT8)
      form='I22'
    case default
      call exit_with_status (1, 'FITS: wrong data type')
  end select
end subroutine

function char2type (char)
  character(len=*), intent(in) :: char
  integer(i4b) :: char2type

  select case (char)
    case ('E','e')
      char2type=FITS_REAL4
    case ('D','d')
      char2type=FITS_REAL8
    case ('B','b')
      char2type=FITS_INT1
    case ('I','i')
      char2type=FITS_INT2
    case ('J','j')
      char2type=FITS_INT4
    case ('K','k')
      char2type=FITS_INT8
    case ('A','a')
      char2type=FITS_CHAR
    case default
      char2type=-1234
      call exit_with_status (1, 'FITS: wrong data type')
  end select
end function

function datacode2type (datacode)
  integer, intent(in) :: datacode
  integer(i4b) :: datacode2type

  select case (datacode)
    case (42)
      datacode2type=FITS_REAL4
    case (82)
      datacode2type=FITS_REAL8
    case (11)
      datacode2type=FITS_INT1
    case (21)
      datacode2type=FITS_INT2
    case (41)
      datacode2type=FITS_INT4
    case (81)
      datacode2type=FITS_INT8
    case (16)
      datacode2type=FITS_CHAR
    case default
      datacode2type=-1234
      call exit_with_status (1, 'FITS: wrong datacode')
  end select
end function

function fixkey (input)
  character(len=*), intent(in) :: input
  character(len=filenamelen) fixkey

  fixkey=input
  call ftupch(fixkey)
  if (fixkey/=input) fixkey = "HIERARCH " // input
end function

subroutine clear_data (handle)
  type(fitshandle), intent(inout) :: handle

  if (handle%hdutype==FITS_IMAGE) then
    if (associated(handle%axes)) deallocate(handle%axes)
  endif
  handle%axes=>NULL()
  if ((handle%hdutype==FITS_BINTAB).or.(handle%hdutype==FITS_ASCTAB)) then
    if (associated(handle%columns)) deallocate(handle%columns)
  endif
  handle%columns=>NULL()
  handle%hdutype=FITS_INVALID
  handle%datatype=FITS_INVALID
  handle%nrows=0
end subroutine

subroutine clearfitshandle (handle)
  type(fitshandle), intent(inout) :: handle

  handle%lun=-1
  call clear_data(handle)
end subroutine

subroutine checkFitsErrors(stat)
  integer, intent(inout) :: stat
  character(len=30) errtext
  character(len=80) errmessage

  errtext=''
  errmessage=''
  if (stat<=0) return
  call ftgerr(stat,errtext)
  print *,'FITSIO Error Status =',stat,': ',errtext

  call ftgmsg(errmessage)
  do while (errmessage /= ' ')
    print *,errmessage
    call ftgmsg(errmessage)
  end do

  call exit_with_status (1,'FITS error')
end subroutine checkFitsErrors

subroutine replaceDate (unit)
  integer, intent(in) :: unit
  integer stat
  character(len=80) dt,fd

  stat=0
  call ftpdat(unit,stat)
  call checkFitsErrors(stat)
  call date_and_time(dt)
  fd=dt(1:4)//'-'//dt(5:6)//'-'//dt(7:8)
  call ftukys(unit,'DATE',fd,'current date',stat)
  call checkFitsErrors(stat)
end subroutine replaceDate

subroutine init_image (handle)
  type(fitshandle), intent(inout) :: handle

  integer(i4b) :: stat, dims

  stat=0
  call ftgidt (handle%lun, handle%datatype, stat)
  call ftgidm (handle%lun, dims, stat)
  call checkFitsErrors(stat)
  allocate (handle%axes(dims))
  call ftgisz (handle%lun, dims, handle%axes, stat)
  call checkFitsErrors(stat)
end subroutine

subroutine init_bintab (handle)
  type(fitshandle), intent(inout) :: handle

  integer(i4b) :: stat, ncols, i
  real(dp) :: tscal, tzero
  integer(i4b) :: tnull
  character(len=80) :: typstr, tdisp

  typstr=''
  tdisp=''
  stat=0
  call ftgncl (handle%lun,ncols,stat)
  call ftgnrwll (handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
  allocate (handle%columns(ncols))
  do i=1,ncols
    call ftgbcl (handle%lun,i,handle%columns(i)%name,handle%columns(i)%unit, &
      typstr,handle%columns(i)%repcount,tscal,tzero,tnull,tdisp,stat)
    call checkFitsErrors(stat)
    handle%columns(i)%type=char2type(trim(typstr))
  end do
end subroutine

subroutine init_asctab (handle)
  type(fitshandle), intent(inout) :: handle

  integer(i4b) :: stat, ncols, i, width, decimals
  integer(i4b), allocatable, dimension(:) :: tbcol
  real(dp) :: tscal, tzero
  character(len=80) :: tform, snull, tdisp
  integer :: datacode

  tform=''
  snull=''
  tdisp=''
  stat=0
  call ftgncl (handle%lun,ncols,stat)
  call ftgnrwll (handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
  allocate (handle%columns(ncols))
  allocate (tbcol(ncols))
  do i=1,ncols
    handle%columns(i)%repcount=1
    call ftgacl (handle%lun,i,handle%columns(i)%name,tbcol, &
      handle%columns(i)%unit,tform,tscal,tzero,snull,tdisp,stat)
    call ftasfm (tform,datacode,width,decimals,stat)
    call checkFitsErrors(stat)
    handle%columns(i)%type=datacode2type(datacode)
  end do
  deallocate(tbcol)
end subroutine

subroutine initData (handle)
  type(fitshandle), intent(inout) :: handle

  select case (handle%hdutype)
    case (FITS_IMAGE)
      call init_image (handle)
    case (FITS_ASCTAB)
      call init_asctab (handle)
    case (FITS_BINTAB)
      call init_bintab (handle)
    case default
      call exit_with_status (1,'FITS: unknown HDU type')
  end select
end subroutine

subroutine fits_goto_hdu (handle, num_hdu)
!FIXME: add a special value (e.g. -1) to go to the last HDU in the file
  type(fitshandle), intent(inout) :: handle
  integer(i4b), intent(in) :: num_hdu

  integer(i4b) :: stat

  stat=0
  call clear_data(handle)
  call ftmahd(handle%lun,num_hdu,handle%hdutype,stat)
  call checkFitsErrors(stat)
  call initData(handle)
end subroutine

subroutine fits_delete_chdu (handle)
  type(fitshandle), intent(inout) :: handle

  integer(i4b) :: stat

  stat=0
  call clear_data(handle)
  call ftdhdu (handle%lun,handle%hdutype,stat)
  call checkFitsErrors (stat)
  call initData (handle)
end subroutine

subroutine fits_insert_bintab (handle,columns)
  type(fitshandle), intent(inout) :: handle
  type(fitscolumn), dimension(:), intent(in) :: columns

  integer(i4b) :: stat, i
  character(len=80), dimension(size(columns)) :: ttype,tform,tunit

  call clear_data(handle)
  allocate (handle%columns(size(columns)))
  handle%columns=columns

  do i=1,size(columns)
    ttype(i)=trim(columns(i)%name)
    write (tform(i),*) columns(i)%repcount
    tform(i)=trim(adjustl(tform(i)))//type2char(columns(i)%type)
    tunit(i)=trim(columns(i)%unit)
  end do

  stat=0
  call ftibin (handle%lun,0_i4b,size(columns),ttype,tform,tunit,'xtension', &
               0,stat)
  call checkFitsErrors (stat)
  handle%hdutype=FITS_BINTAB
end subroutine

subroutine fits_insert_asctab (handle,columns)
  type(fitshandle), intent(inout) :: handle
  type(fitscolumn), dimension(:), intent(in) :: columns

  integer(i4b) :: stat, i, rowlen
  character(len=80), dimension(size(columns)) :: ttype,tform,tunit
  integer(i4b), dimension(size(columns)) :: tbcol

  call clear_data(handle)
  allocate (handle%columns(size(columns)))
  handle%columns=columns

  do i=1,size(columns)
    ttype(i)=trim(columns(i)%name)
    call assert(columns(i)%repcount==1,'repcount must be 1 in ASCII tables')
    write (tform(i),*) columns(i)%repcount
    call type2asciiform(columns(i)%type,tform(i))
    tunit(i)=trim(columns(i)%unit)
  end do

  stat=0
  call ftgabc(size(columns),tform,1,rowlen,tbcol,stat)
  call ftitab (handle%lun,rowlen,0_i4b,size(columns),ttype,tbcol,tform,tunit, &
               'xtension', stat)
  call checkFitsErrors (stat)
  handle%hdutype=FITS_BINTAB
end subroutine

subroutine fits_insert_image (handle, datatype, axes)
  type(fitshandle), intent(inout) :: handle
  integer(i4b), intent(in) :: datatype
  integer(i4b), intent(in), dimension(:) :: axes

  integer(i4b) :: stat

  call clear_data(handle)
  allocate (handle%axes(size(axes)))
  handle%axes=axes
  handle%datatype=datatype
  stat=0
  call ftiimg (handle%lun,datatype,size(axes),axes,stat)
  call checkFitsErrors (stat)
  handle%hdutype=FITS_IMAGE
end subroutine

subroutine fits_close (handle)
  type(fitshandle), intent(inout) :: handle

  integer(i4b) :: stat

  stat=0
  call ftclos(handle%lun,stat)
  call ftfiou(handle%lun,stat)
  handle%lun=-1
  call checkFitsErrors(stat)
  call clearfitshandle (handle)
end subroutine

function is_extended_syntax (file_name)
  ! Determine if FILE_NAME uses CFITSIO's extended syntax
  ! to specify the HDU (FILENAME[n] or FILENAME+n).
  character(len=*), intent(in) :: file_name
  logical :: is_extended_syntax
  integer :: idx

  idx = len_trim(file_name)
  is_extended_syntax = &
    ((file_name(idx:idx) .eq. ']') .or. (index(file_name, '+') .gt. 0))
end function

subroutine fits_open (handle, file, rwmode, num_hdu)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: file
  integer(i4b), intent(in), optional :: rwmode, num_hdu

  integer(i4b) :: rwm, nhdu, stat, blocksize

  call clearfitshandle(handle)

  rwm=FITS_READONLY
  if (present(rwmode)) rwm=rwmode

  stat=0
  call ftgiou(handle%lun,stat)

  if (.not. is_extended_syntax(file)) then
    nhdu=2
    if (present(num_hdu)) nhdu=num_hdu
    call ftopen(handle%lun,trim(file),rwm,blocksize,stat)
    call ftmahd(handle%lun,nhdu,handle%hdutype,stat)
  else
    call assert(.not. present(num_hdu), &
      "cannot use 'num_hdu' and extended syntax at the same time")
    ! We are using the extended syntax, so the HDU is encoded in the file name
    call ftnopn(handle%lun,trim(file),rwm,stat)
    ! Determine the HDU type
    call ftghdt(handle%lun,handle%hdutype,stat)
  end if

  call checkFitsErrors(stat)
  call initData(handle)
end subroutine

subroutine fits_delete (file)
  character(len=*), intent(in) :: file

  integer(i4b) :: unit,stat,blocksize

  stat=0
  call ftgiou(unit,stat)
  call ftopen(unit,trim(file), 1,blocksize,stat)
  call ftdelt(unit,stat)
  call ftfiou(unit,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_create (handle, file)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: file

  integer(i4b) :: stat

  stat=0
  call ftgiou(handle%lun,stat)
  call ftinit(handle%lun,trim(file),FITS_READWRITE,stat)
  call ftphps(handle%lun,8,0,(/0,0/),stat)
  call checkFitsErrors(stat)
  call replaceDate(handle%lun)
end subroutine

subroutine split_comment (c1,c2,pos,comment)
  character(len=*), intent(out) :: c1,c2
  integer, intent(in) :: pos
  character(len=*), intent(in), optional :: comment

  c1=''
  c2=''

  if (.not. present(comment)) return

  if (len_trim(comment)<=pos) then
    c1=comment
  else
    c1=comment(1:pos)
    c2=comment(pos+1:)
  endif
end subroutine

subroutine fits_add_key_int1 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i1b), intent(in) :: value

  call fits_add_key (handle,key,int(value,i4b),comment)
end subroutine

subroutine fits_upd_key_int1 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i1b), intent(in) :: value

  call fits_update_key (handle,key,int(value,i4b),comment)
end subroutine

subroutine fits_add_key_int2 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i2b), intent(in) :: value

  call fits_add_key (handle,key,int(value,i4b),comment)
end subroutine

subroutine fits_upd_key_int2 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i2b), intent(in) :: value

  call fits_update_key (handle,key,int(value,i4b),comment)
end subroutine

subroutine fits_add_key_int4 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i4b), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftpkyj (handle%lun,trim(fixkey(key)),value,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_upd_key_int4 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i4b), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftukyj (handle%lun,trim(fixkey(key)),value,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_add_key_int8 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i8b), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftpkyk (handle%lun,trim(fixkey(key)),value,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_upd_key_int8 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  integer(i8b), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftukyk (handle%lun,trim(fixkey(key)),value,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_add_key_real4 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  real(sp), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftpkye (handle%lun,trim(fixkey(key)),value,-15,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_upd_key_real4 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  real(sp), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftukye (handle%lun,trim(fixkey(key)),value,-15,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_add_key_real8 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  real(dp), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftpkyd (handle%lun,trim(fixkey(key)),value,-20,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_upd_key_real8 (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  real(dp), intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftukyd (handle%lun,trim(fixkey(key)),value,-20,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_add_key_string (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  character(len=*), intent(in) :: value

  integer(i4b) :: stat, pos
  character(len=comlen) :: c1,c2

  stat=0
  pos = min(47,65-len_trim(value))
  call split_comment (c1,c2,pos,comment)
  call ftpkys (handle%lun,trim(fixkey(key)),trim(value),trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_upd_key_string (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  character(len=*), intent(in) :: value

  integer(i4b) :: stat, pos
  character(len=comlen) :: c1,c2

  stat=0
  pos = min(47,65-len_trim(value))
  call split_comment (c1,c2,pos,comment)
  call ftukys (handle%lun,trim(fixkey(key)),trim(value),trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_add_key_bool (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  logical, intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftpkyl (handle%lun,trim(fixkey(key)),value,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_upd_key_bool (handle,key,value,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(in), optional :: comment
  logical, intent(in) :: value

  integer(i4b) :: stat
  character(len=comlen) :: c1,c2

  stat=0
  call split_comment (c1,c2,47,comment)
  call ftukyl (handle%lun,trim(fixkey(key)),value,trim(c1),stat)
  if (c2/='') call ftpcom (handle%lun,trim(c2),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_delete_key (handle,key)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key

  integer(i4b) :: stat

  stat=0
  call ftdkey (handle%lun,trim(key),stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_add_comment (handle,comment)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: comment

  integer(i4b) :: stat

  stat=0
  call ftpcom (handle%lun,trim(comment),stat)
  call checkFitsErrors(stat)
end subroutine

function fits_key_present (handle, key)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical fits_key_present

  integer(i4b) :: stat
  character(len=80) card

  fits_key_present=.true.
  stat=0
  card=''
  call ftgcrd(handle%lun,trim(key),card,stat)
  if (stat==202) then
    fits_key_present=.false.
    stat=0
    call ftcmsg
  endif
  call checkFitsErrors(stat)
end function

subroutine fits_get_key_int1 (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i1b), intent(out) :: value

  integer(i4b) :: tval

  call fits_get_key (handle,key,tval)
  value = int(tval,i1b)
end subroutine

subroutine fits_get_key_int2 (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i2b), intent(out) :: value

  integer(i4b) :: tval

  call fits_get_key (handle,key,tval)
  value = int(tval,i2b)
end subroutine

subroutine fits_get_key_int4 (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i4b), intent(out) :: value

  integer(i4b) :: stat
  character(len=80) comment

  comment=''
  stat=0
  call ftgkyj (handle%lun,trim(key),value,comment,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_get_key_int8 (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  integer(i8b), intent(out) :: value

  integer(i4b) :: stat
  character(len=80) comment

  comment=''
  stat=0
  call ftgkyk (handle%lun,trim(key),value,comment,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_get_key_real4 (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(sp), intent(out) :: value

  integer(i4b) :: stat
  character(len=80) comment

  comment=''
  stat=0
  call ftgkye (handle%lun,trim(key),value,comment,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_get_key_real8 (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  real(dp), intent(out) :: value

  integer(i4b) :: stat
  character(len=80) comment

  comment=''
  stat=0
  call ftgkyd (handle%lun,trim(key),value,comment,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_get_key_string (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  character(len=*), intent(out) :: value

  integer(i4b) :: stat
  character(len=80) comment

  comment=''
  stat=0
  value=''
  call ftgkys (handle%lun,trim(key),value,comment,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_get_key_bool (handle,key,value)
  type(fitshandle), intent(inout) :: handle
  character(len=*), intent(in) :: key
  logical, intent(out) :: value

  integer(i4b) :: stat
  character(len=80) comment

  comment=''
  stat=0
  call ftgkyl (handle%lun,trim(key),value,comment,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine check_prep_write (handle, col, frow, felem, offset)
  type(fitshandle), intent(in) :: handle
  integer, intent(in) :: col
  integer(i8b), intent(out) :: frow, felem
  integer(i8b), intent(in), optional :: offset

  integer(i8b) ofs

  ofs=0
  if (present(offset)) ofs=offset

  call assert((handle%hdutype==FITS_BINTAB).or. &
              (handle%hdutype==FITS_ASCTAB), 'wrong HDU type')

  frow  = ofs/handle%columns(col)%repcount+1
  felem = mod(ofs,int(handle%columns(col)%repcount,i8b))+1
end subroutine

subroutine check_prep_read (handle, col, datasize, frow, felem, offset)
  type(fitshandle), intent(in) :: handle
  integer, intent(in) :: col, datasize
  integer(i8b), intent(out) :: frow, felem
  integer(i8b), intent(in), optional :: offset

  integer(i8b) ofs

  ofs=0
  if (present(offset)) ofs=offset

  call assert((handle%hdutype==FITS_BINTAB).or. &
              (handle%hdutype==FITS_ASCTAB), 'wrong HDU type')

  frow  = ofs/handle%columns(col)%repcount+1
  felem = mod(ofs,int(handle%columns(col)%repcount,i8b))+1

  call assert (datasize<=handle%columns(col)%repcount*handle%nrows-ofs, &
    'read_column: array too large', 1)
end subroutine

subroutine fits_write_column_real4 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(in) :: data(:)
  integer(i8b), intent(in), optional :: offset

  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_write (handle, colnum, frow, felem, offset)
  stat=0
  call ftpclell(handle%lun,colnum,frow,felem,size(data),data(:),stat)
  call ftgnrwll(handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_column_real8 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(in) :: data(:)
  integer(i8b), intent(in), optional :: offset

  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_write (handle, colnum, frow, felem, offset)
  stat=0
  call ftpcldll(handle%lun,colnum,frow,felem,size(data),data(:),stat)
  call ftgnrwll(handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_column_int1 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(in) :: data(:)
  integer(i8b), intent(in), optional :: offset

  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_write (handle, colnum, frow, felem, offset)
  stat=0
  call ftpclbll(handle%lun,colnum,frow,felem,size(data),data(:),stat)
  call ftgnrwll(handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_column_int2 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(in) :: data(:)
  integer(i8b), intent(in), optional :: offset

  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_write (handle, colnum, frow, felem, offset)
  stat=0
  call ftpclill(handle%lun,colnum,frow,felem,size(data),data(:),stat)
  call ftgnrwll(handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_column_int4 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(in) :: data(:)
  integer(i8b), intent(in), optional :: offset

  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_write (handle, colnum, frow, felem, offset)
  stat=0
  call ftpcljll(handle%lun,colnum,frow,felem,size(data),data(:),stat)
  call ftgnrwll(handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_column_int8 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(in) :: data(:)
  integer(i8b), intent(in), optional :: offset

  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_write (handle, colnum, frow, felem, offset)
  stat=0
  call ftpclkll(handle%lun,colnum,frow,felem,size(data),data(:),stat)
  call ftgnrwll(handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_column_string (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(in) :: data(:)
  integer(i8b), intent(in), optional :: offset

  integer stat, stringlen
  integer(i8b) ofs, frow, felem

  if (size(data)==0) return
  ofs=0
  if (present(offset)) ofs=offset
  stringlen=handle%columns(colnum)%repcount
  call check_prep_write (handle, colnum, frow, felem, ofs*stringlen)
  call assert(felem==1,'offset error in fits_write_column_string()')

  stat=0
  call ftpclsll(handle%lun,colnum,frow,felem,size(data),data(:),stat)
  call ftgnrwll(handle%lun,handle%nrows,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_column_real4 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(sp), intent(out) :: data(:)
  integer(i8b), intent(in), optional :: offset

  logical anynul
  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_read (handle, colnum, size(data), frow, felem, offset)

  stat=0
  call ftgcvell(handle%lun,colnum,frow,felem,size(data),0._sp,data(:),anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_column_real8 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  real(dp), intent(out) :: data(:)
  integer(i8b), intent(in), optional :: offset

  logical anynul
  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_read (handle, colnum, size(data), frow, felem, offset)

  stat=0
  call ftgcvdll(handle%lun,colnum,frow,felem,size(data),0._dp,data(:),anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_column_int1 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i1b), intent(out) :: data(:)
  integer(i8b), intent(in), optional :: offset

  logical anynul
  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_read (handle, colnum, size(data), frow, felem, offset)

  stat=0
  call ftgcvbll(handle%lun,colnum,frow,felem,size(data),0_i2b,data(:),anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_column_int2 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i2b), intent(out) :: data(:)
  integer(i8b), intent(in), optional :: offset

  logical anynul
  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_read (handle, colnum, size(data), frow, felem, offset)

  stat=0
  call ftgcvill(handle%lun,colnum,frow,felem,size(data),0_i2b,data(:),anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_column_int4 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i4b), intent(out) :: data(:)
  integer(i8b), intent(in), optional :: offset

  logical anynul
  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_read (handle, colnum, size(data), frow, felem, offset)

  stat=0
  call ftgcvjll(handle%lun,colnum,frow,felem,size(data),0_i4b,data(:),anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_column_int8 (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  integer(i8b), intent(out) :: data(:)
  integer(i8b), intent(in), optional :: offset

  logical anynul
  integer stat
  integer(i8b) :: frow, felem

  if (size(data)==0) return
  call check_prep_read (handle, colnum, size(data), frow, felem, offset)

  stat=0
  call ftgcvkll(handle%lun,colnum,frow,felem,size(data),0_i8b,data(:),anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_column_string (handle, colnum, data, offset)
  type(fitshandle), intent(inout) :: handle
  integer, intent(in) :: colnum
  character(len=*), intent(out) :: data(:)
  integer(i8b), intent(in), optional :: offset

  logical anynul
  integer stat, stringlen
  integer(i8b) ofs, frow, felem

  if (size(data)==0) return
  ofs=0
  if (present(offset)) ofs=offset
  stringlen=handle%columns(colnum)%repcount
  call check_prep_read (handle, colnum, size(data)*stringlen, frow, felem, &
                        ofs*stringlen)
  call assert(felem==1,'offset error in fits_read_column_string()')

  stat=0
  data=''
  call ftgcvsll(handle%lun,colnum,frow,felem,size(data),'',data(:),anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_image_real4 (handle, data)
  type(fitshandle), intent(inout) :: handle
  real(sp), intent(in) :: data(:,:)

  integer stat

  call assert (handle%hdutype==FITS_IMAGE,'wrong HDU type')
  call assert (size(handle%axes)==2,'wrong number of dimensions')
  call assert (size(data,1)==handle%axes(1),'wrong size of dimension 1')
  call assert (size(data,2)==handle%axes(2),'wrong size of dimension 2')

  stat=0
  call ftp2de(handle%lun,0,handle%axes(1),handle%axes(1),handle%axes(2), &
              data,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_image_real8 (handle, data)
  type(fitshandle), intent(inout) :: handle
  real(dp), intent(in) :: data(:,:)

  integer stat

  call assert (handle%hdutype==FITS_IMAGE,'wrong HDU type')
  call assert (size(handle%axes)==2,'wrong number of dimensions')
  call assert (size(data,1)==handle%axes(1),'wrong size of dimension 1')
  call assert (size(data,2)==handle%axes(2),'wrong size of dimension 2')

  stat=0
  call ftp2dd(handle%lun,0,handle%axes(1),handle%axes(1),handle%axes(2), &
              data,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_write_image_int4 (handle, data)
  type(fitshandle), intent(inout) :: handle
  integer(i4b), intent(in) :: data(:,:)

  integer stat

  call assert (handle%hdutype==FITS_IMAGE,'wrong HDU type')
  call assert (size(handle%axes)==2,'wrong number of dimensions')
  call assert (size(data,1)==handle%axes(1),'wrong size of dimension 1')
  call assert (size(data,2)==handle%axes(2),'wrong size of dimension 2')

  stat=0
  call ftp2dj(handle%lun,0,handle%axes(1),handle%axes(1), handle%axes(2), &
              data,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_image_real4 (handle, data)
  type(fitshandle), intent(inout) :: handle
  real(sp), intent(out) :: data(:,:)

  integer stat
  logical anynul

  call assert (handle%hdutype==FITS_IMAGE,'wrong HDU type')
  call assert (size(handle%axes)==2,'wrong number of dimensions')
  call assert (size(data,1)==handle%axes(1),'wrong size of dimension 1')
  call assert (size(data,2)==handle%axes(2),'wrong size of dimension 2')

  stat=0
  call ftg2de(handle%lun,0,0._sp,handle%axes(1),handle%axes(1), &
              handle%axes(2),data,anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_image_real8 (handle, data)
  type(fitshandle), intent(inout) :: handle
  real(dp), intent(out) :: data(:,:)

  integer stat
  logical anynul

  call assert (handle%hdutype==FITS_IMAGE,'wrong HDU type')
  call assert (size(handle%axes)==2,'wrong number of dimensions')
  call assert (size(data,1)==handle%axes(1),'wrong size of dimension 1')
  call assert (size(data,2)==handle%axes(2),'wrong size of dimension 2')

  stat=0
  call ftg2dd(handle%lun,0,0._dp,handle%axes(1),handle%axes(1), &
              handle%axes(2),data,anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_image_int4 (handle, data)
  type(fitshandle), intent(inout) :: handle
  integer(i4b), intent(out) :: data(:,:)

  integer stat
  logical anynul

  call assert (handle%hdutype==FITS_IMAGE,'wrong HDU type')
  call assert (size(handle%axes)==2,'wrong number of dimensions')
  call assert (size(data,1)==handle%axes(1),'wrong size of dimension 1')
  call assert (size(data,2)==handle%axes(2),'wrong size of dimension 2')

  stat=0
  call ftg2dj(handle%lun,0,0_i4b,handle%axes(1),handle%axes(1), &
              handle%axes(2),data,anynul,stat)
  call checkFitsErrors(stat)
end subroutine

subroutine fits_read_header (handle, buf)
  type(fitshandle), intent(inout) :: handle
  type(linebuf), intent(inout) :: buf

  integer(i4b) :: i, stat
  character(len=80) cd
  character(len=8), dimension(1), parameter ::  inclist = (/"*"/)
  character(len=8), dimension(23) :: exclist

  exclist = &
    (/"SIMPLE  ","BITPIX  ","NAXIS   ","NAXIS#  ","PCOUNT  ","GCOUNT  ", &
      "EXTEND  ","ORIGIN  ","DATE*   ","TFIELDS ","TTYPE#  ","TFORM#  ", &
      "TUNIT#  ","EXTNAME ","CTYPE#  ","CRVAL#  ","CRPIX#  ","CDELT#  ", &
      "XTENSION","INSTRUME","TELESCOP","PDMTYPE ","TBCOL#  "/)

  cd=''
  stat=0
  call ftgrec(handle%lun,0,cd,stat)
  call checkFitsErrors(stat)

  i=0
  do
    i=i+1
    call ftgnxk(handle%lun,inclist,size(inclist),exclist,size(exclist),cd,stat)
    if (stat>0) exit
    call linebuf_add(buf,cd)
  end do

  call ftcmsg
end subroutine

subroutine fits_write_header (handle, buf)
  type(fitshandle), intent(inout) :: handle
  type(linebuf), intent(in) :: buf

  integer i, stat

  stat=0
  do i=1,buf%size
    select case (trim(buf%data(i)(1:8)))
      case ('END')
!       do nothing
      case default
        call ftprec(handle%lun,buf%data(i),stat)
        call checkFitsErrors(stat)
    end select
  end do
end subroutine

function fits_fast_nelms (handle, col)
  type(fitshandle), intent(in) :: handle
  integer, intent(in) :: col
  integer(i8b) :: fits_fast_nelms

  integer stat,nelms

  stat=0
  call assert (handle%hdutype==FITS_BINTAB, &
               'fits_fast_nelms: wrong HDUTYPE')
  call assert (col<=size(handle%columns),'fits_fast_nelms: bad column number')
  call ftgrsz(handle%lun,nelms,stat)
  call checkFitsErrors(stat)
  call assert (nelms>=1,"bad recommended FITS chunk size")
  fits_fast_nelms=nelms*handle%columns(col)%repcount
end function

function fits_num_hdus (handle)
  type(fitshandle), intent(in) :: handle
  integer fits_num_hdus

  integer stat

  stat=0
  call ftthdu(handle%lun,fits_num_hdus,stat)
  call checkFitsErrors(stat)
end function

function fitstype2type (ftype)
  integer(i4b), intent(in) :: ftype
  integer(i4b) :: fitstype2type

  select case (ftype)
    case (FITS_INT1)
      fitstype2type = PLANCK_INT8
    case (FITS_INT2)
      fitstype2type = PLANCK_INT16
    case (FITS_INT4)
      fitstype2type = PLANCK_INT32
    case (FITS_INT8)
      fitstype2type = PLANCK_INT64
    case (FITS_REAL4)
      fitstype2type = PLANCK_FLOAT32
    case (FITS_REAL8)
      fitstype2type = PLANCK_FLOAT64
    case (FITS_CHAR)
      fitstype2type = PLANCK_STRING
    case default
      fitstype2type = PLANCK_INVALID_TYPE
      call exit_with_status(1,"fitstype2type: unknown FITS type")
  end select
end function

function type2fitstype (type)
  integer(i4b), intent(in) :: type
  integer(i4b) :: type2fitstype

  select case (type)
    case (PLANCK_INT8)
      type2fitstype = FITS_INT1
    case (PLANCK_INT16)
      type2fitstype = FITS_INT2
    case (PLANCK_INT32)
      type2fitstype = FITS_INT4
    case (PLANCK_INT64)
      type2fitstype = FITS_INT8
    case (PLANCK_FLOAT32)
      type2fitstype = FITS_REAL4
    case (PLANCK_FLOAT64)
      type2fitstype = FITS_REAL8
    case (PLANCK_STRING)
      type2fitstype = FITS_CHAR
    case default
      type2fitstype = -1234
      call exit_with_status(1,"type2fitstype: unknown type")
  end select
end function

end module fitsmod2
