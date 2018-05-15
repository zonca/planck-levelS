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
MODULE ZOD_BlackBody
!
! ZOD_BlackBody 0.0 By M. Maris - 8 Mar 2001 -
!
! Class to handle the Black Body Emissivity for ZODIACAL LIGHT
!

Use ZOD_Debug
Use Physical_Parameters_CGS
implicit none

Character(len=*), Parameter :: ZOD_BlackBody_Version='ZOD_BlackBody 0.0 By M. Maris - 8 Mar 2001 -'

CONTAINS

Double Precision function bbl_cgs(lambda,T)
!
! BB(lambda,T) thermal radiance function cgs
!
! lambda in cm
!
! bbl_cgs = erg/(cm2 sec cm sterad)
!

    Double Precision, Intent(IN) :: Lambda
    Double Precision, Intent(IN) :: T

! Local Variables
    Double Precision bbl, FT, ET

    if (T.le.0) call die('ZOD_BlackBody::bbl_cgs','Negative or null T')
    if (Lambda.le.0) call die('ZOD_BlackBody::bbl_cgs','Negative or null Lambda')

    FT = lambda*lambda*lambda*lambda*lambda
    FT = 2.*cgs_h*cgs_c*cgs_c/FT;
    ET = exp( cgs_h*cgs_c/(lambda*cgs_k*T) );
    ET = ET - 1.

    bbl = FT / ET

    bbl_cgs = bbl
end function bbl_cgs

Double Precision function bbn_cgs(nu,T)
!
! BB(nu,T) thermal radiance function cgs
!
! nu in Hz
!
! bbn_cgs = erg/(cm2 sec cm sterad)
!
    Double Precision, Intent(IN) :: Nu
    Double Precision, Intent(IN) :: T

! Local Variables
    Double Precision bbn, FT, ET

    if (T.le.0) call die('ZOD_BlackBody::bbn_cgs','Negative or null T')
    if (Nu.le.0) call die('ZOD_BlackBody::bbn_cgs','Negative or null Nu')

    FT = Nu*1.d0
    FT = FT*FT*FT

    ET = cgs_c*cgs_c
    FT = 2.d0*cgs_h*FT/ET
    ET = cgs_h*nu/(cgs_k*T)
    ET = exp(ET)
    ET = ET - 1.d0

    bbn = FT / ET

    bbn_cgs = bbn
end function bbn_cgs

Double Precision function bbl_rj_cgs(lambda,T)
!
! BB thermal radiance function cgs in Rayleight-Jeans approx
!
! lambda in cm
! T in K
! rj_cgs = erg/(cm2 sec cm sterad)
!
    Double Precision, Intent(IN) :: Lambda
    Double Precision, Intent(IN) :: T

! Local Variables
    Double Precision rj, FT

    if (T.le.0) call die('ZOD_BlackBody::bbl_rj_cgs','Negative or null T')
    if (Lambda.le.0) call die('ZOD_BlackBody::bbl_rj_cgs','Negative or null Lambda')

    FT = lambda*lambda*lambda*lambda
    rj = 2.*cgs_k*T*cgs_c/FT

    bbl_rj_cgs = rj
end function bbl_rj_cgs

Double Precision function bbn_rj_cgs(nu,T)
!
! BB thermal radiance function cgs in Rayleight-Jeans approx
!
! nu in Hz
! T in K
!
! rj_cgs = erg/(cm2 sec Hz sterad)
!
    Double Precision, Intent(IN) :: Nu
    Double Precision, Intent(IN) :: T

! Local Variables
    Double Precision rj, FT

    if (T.le.0) call die('ZOD_BlackBody::bbn_rj_cgs','Negative or null T')
    if (nu.le.0) call die('ZOD_BlackBody::bbn_rj_cgs','Negative or null Lambda')

    FT = (nu/cgs_c)
    FT = FT*FT
    rj = 2.* cgs_k*T*FT
    bbn_rj_cgs = rj
end function bbn_rj_cgs

END MODULE ZOD_BlackBody
