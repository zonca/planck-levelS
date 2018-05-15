!-----------------------------------------------------------------------------
!
!  Copyright (C) 2001-2013 Michele Maris
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
Module ZOD_Debug
!
! ZOD_Debug 0.0 By M. Maris - 8 Mar 2001 -
!
! This module handles trapped exceptions
!
IMPLICIT NONE

Character(len=*), Parameter :: ZOD_Debug_Version = 'ZOD_Debug 0.0 By M. Maris - 8 Mar 2001 -'

Contains

Subroutine Die(Caller,Motivation,Stat)
!
! Close the program producing an error message
!

    Character(len=*),           Intent(IN) :: Caller
    Character(len=*),           Intent(IN) :: Motivation
    Integer         , Optional, Intent(IN) :: Stat

    print*
    print*,'### Error in : ',Caller
    print*,'### Reason   : ',Motivation
    if (present(Stat)) print*,'### Status   : ',Stat
    print*
    stop
End Subroutine Die

End Module ZOD_Debug
