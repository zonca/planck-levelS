!-----------------------------------------------------------------------------
!
!  Copyright (C) 2002-2013 Michele Maris
!
!  This file is part of the zodiacal light component of the Planck simulation
!  package.
!
!  This code is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This code is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this code; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!-----------------------------------------------------------------------------
Module Physical_Parameters_CGS

!
! Physical_Parameters_CGS 1.0 By M. Maris - 8 Mar 2001 - 10 Sep 2002 - 8 Jul 2004 - 12 Jul 2004 -
!
! Physical parameters in CGS
!
! :Update: M. Maris : 10 Sep 2002 :
!    Physical_Parameters_CGS_Version introduced,  Physical_Parameters_CGS_Display updated.
!
! :Update: M. Maris : 8 Jul 2004 -
!   Physical parameters are in double precision
!
! :Update: M. Maris : 12 Jul 2004 -

IMPLICIT NONE

Character(len=*), Parameter :: Physical_Parameters_CGS_Version &
 & = 'Physical_Parameters_CGS 1.0 By M. Maris - 8 Mar 2001 - 12 July 2004 -'

Double Precision, Parameter :: CGS_K = 1.380658d-16    ! erg/K
Double Precision, Parameter :: CGS_C = 2.99792458d10   ! cm/sec
Double Precision, Parameter :: CGS_h = 6.6260755d-27   ! erg / sec i.e. erg Hz
Double Precision, Parameter :: CGS_flux_Jansky = 1d23  ! (cm2 sec cm sterad)/erg * Jansky

Double Precision, Parameter :: pi = 3.14159265358979323846234d0 ! PI GREEK

Double Precision, Parameter :: ScaleToMJ = 1.d17 ! 1.e23 Jy/erg cm2 sec * 1e-6 MJy/Jy

CONTAINS

subroutine Physical_Parameters_CGS_Display(UNIT,PFX)
!
! Subroutine to print on the screen the list of Physical Parameters
! in cgs
!
   Integer, Optional, Intent(IN) :: Unit
   Character(len=*), Optional, Intent(IN) :: PFX
   Character(len=10) :: LPFX

   if (present(Unit)) then
      if (present(PFX)) then
         LPFX = PFX
      else
         LPFX=''
      endif
      write(Unit,*) trim(LPFX),'CGS_K           = ',CGS_K,' erg/K'
      write(Unit,*) trim(LPFX),'CGS_C           = ',CGS_C,' cm/sec'
      write(Unit,*) trim(LPFX),'CGS_H           = ',CGS_H,' erg/sec'
      write(Unit,*) trim(LPFX),'CGS_FLUX_JANSKY = ',CGS_FLUX_JANSKY,' (cm2 sec cm sterad)/erg * Jansky'
      write(Unit,*) trim(LPFX),'ScaleToMJ       = ',ScaleToMJ,' 1.e23 Jy/erg cm2 sec * 1e-6 MJy/Jy'
      write(Unit,*) trim(LPFX),'PI              = ',pi,' '
   else
      print*,'CGS_K           = ',CGS_K,' erg/K'
      print*,'CGS_C           = ',CGS_C,' cm/sec'
      print*,'CGS_H           = ',CGS_H,' erg/sec'
      print*,'CGS_FLUX_JANSKY = ',CGS_FLUX_JANSKY,' (cm2 sec cm sterad)/erg * Jansky'
      print*,'ScaleToMJ       = ',ScaleToMJ,' 1.e23 Jy/erg cm2 sec * 1e-6 MJy/Jy'
      print*,'PI              = ',pi,' '
   endif

end subroutine Physical_Parameters_CGS_Display


End Module Physical_Parameters_CGS
