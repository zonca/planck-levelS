!-----------------------------------------------------------------------------
!
!  Copyright (C) 1999-2013 Daniel Mortlock
!
!  This file is part of the "simmission" component of the Planck simulation
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

module solarsystem_solsys
  use planck_config
  use general_const
  use general_error
  use general_time
  use solarsystem_star
  use solarsystem_orbit
  use solarsystem_planet
  use solarsystem_l2orbit
  implicit none
  private

  public :: sssolsys, ss_solsys_init, ss_planet_init, &
    ss_planet_name, ss_lambda_star

  ! Structure containing full information on the Solar system.
  type sssolsys
    integer :: npl
    type(ssplanet), pointer :: planets(:) => null()
    type(ssstar) :: star
  end type sssolsys

  ! Length conversions between astronomical units (AU), kilometres (KM)
  ! and metres (M).
  real(dp), parameter :: GNAU_M = 1.49597892e11_dp

  ! Solar paramters, mainly taken from the general package, but given
  ! a distinct name here as the GN prefix is indicative of their use
  ! as astronomical units as opposed to the properties of a particular
  ! object.
  real(dp), parameter :: &
    SSM_SUN_KG = 1.9891e30_dp, &
    SSR_SUN_M = 6.9599e8_dp, &
    SSRADFORCE_SUN_N = 4.5d-6 * (fourpi * GNAU_M**2)

  real(dp), parameter :: &
    GNJ2000_CENT = 20.0_dp

  real(dp), parameter :: &
    GNARCSEC_DEG = 0.0002777777777777777777777777777777777777778_dp

  ! Classic Keplerian orbital elements for the nine major planets (including
  ! Pluto), correct for the epoch J2000 (JED 2451545.0), as taken from
  ! Table 5.8.1 of `The Explanatory Supplement To The Astronomical Almanac'.
  ! The units are those given in the above reference, but are immediately
  ! converted to the package's default units (m for length; kg for mass;
  ! radians for angles). In the case of the Earth, parameters are given for
  ! the Earth-Moon system, and the almanac's modifications to the 1976 IAU
  ! values have been adopted.  Each parameter is explained in more detail
  ! below.
  integer, parameter :: SSNPL = 9
  real(dp) :: SST_0_CENT = GNJ2000_CENT

  ! Common English name for the planet.
  character(len=7), parameter :: SSNAME_PL(SSNPL) &
    = (/ &
    'Mercury', &
    'Venus  ', &
    'Earth  ', &
    'Mars   ', &
    'Jupiter', &
    'Saturn ', &
    'Uranus ', &
    'Neptune', &
    'Pluto  ' &
    /)

  ! Ratio of the Sun's mass to that of the planet.
  real(dp), parameter :: SSM_SUN_M_PL(SSNPL) &
    = (/ &
    6023600.0d0, &
    408523.5d0, &
    328900.55_dp, &
    3098710.0d0, &
    1047.350d0, &
    3498.0d0, &
    22960.d0, &
    19314d0, &
    130000000.0d0 &
    /)

  ! Mass of the planet (in kg, assuming the Sun's mass is
  ! SSM_SUN_KG = 1.9891 x 10^30 kg).
  real(dp), parameter :: SSM_PL_KG(SSNPL) &
    = (/ &
    SSM_SUN_KG / SSM_SUN_M_PL(1), &
    SSM_SUN_KG / SSM_SUN_M_PL(2), &
    SSM_SUN_KG / SSM_SUN_M_PL(3), &
    SSM_SUN_KG / SSM_SUN_M_PL(4), &
    SSM_SUN_KG / SSM_SUN_M_PL(5), &
    SSM_SUN_KG / SSM_SUN_M_PL(6), &
    SSM_SUN_KG / SSM_SUN_M_PL(7), &
    SSM_SUN_KG / SSM_SUN_M_PL(8), &
    SSM_SUN_KG / SSM_SUN_M_PL(9) &
    /)

  ! Average radius of the planet (in m).
  real(dp), parameter :: SSR_PL_M(SSNPL) &
    = (/ &
    2439700.0, &
    6051900.0, &
    6378140.0, &
    3397000.0, &
    71492000.0, &
    60268000.0, &
    25559000.0, &
    24764000.0, &
    1151000.0 &
    /)

  ! Semi-major axis of the planet's orbit (in AU).
  real(dp), parameter :: SSA_PL_AU(SSNPL) &
    = (/ &
    0.38709893d0, &
    0.7233199d0, &
    1.00000011d0, &
    1.52366231d0, &
    5.20336301d0, &
    9.53707032d0, &
    19.19126393d0, &
    30.06896348d0, &
    39.48168677d0 &
    /)

  ! Rate of change of the planet's semi-major axis [in AU cent^(-1)].
  real(dp), parameter :: SSDADT_PL_AUCENT(SSNPL) &
    = (/ &
    0.00000066d0, &
    0.00000092d0, &
    - 0.00000005d0, &
    - 0.00007221d0, &
    0.00060737d0, &
    - 0.00301530d0, &
    0.00152025d0, &
    - 0.00125196d0, &
    - 0.0076912d0 &
    /)

  ! Eccentricity of the planet's orbit.
  real(dp), parameter :: SSE_PL(SSNPL) &
    = (/ &
    0.20563069d0, &
    0.00677323d0, &
    0.01671022d0, &
    0.09341233d0, &
    0.04839266d0, &
    0.05415060d0, &
    0.04716771d0, &
    0.00858587d0, &
    0.24880766d0 &
    /)

  ! Rate of change of the eccentricity of the planet's orbit [in cent^(-1)].
  real(dp), parameter :: SSDEDT_PL_CENT(SSNPL) &
    = (/ &
    0.00002527d0, &
    - 0.00004938d0, &
    - 0.00003804d0, &
    0.00011902d0, &
    - 0.00012880d0, &
    - 0.00036762d0, &
    - 0.00019150d0, &
    0.00002514d0, &
    0.00006465d0 &
    /)

  ! The inclination of the planet's orbit, relative to the ecliptic (in deg).
  real(dp), parameter :: SSINC_PL_DEG(SSNPL) &
    = (/ &
    7.00487d0, &
    3.39471d0, &
    0.00005d0, &
    1.85061d0, &
    1.30530d0, &
    2.48446d0, &
    0.76986d0, &
    1.76917d0, &
    17.14175d0 &
    /)

  ! The rate of change of the inclination of the planet's orbit
  ! [in arcsec cent^(-1)].
  real(dp), parameter :: SSDINCDT_PL_ARCSECCENT(SSNPL) &
    = (/ &
    - 23.51d0, &
    - 2.86d0, &
    - 46.94d0, &
    - 25.47d0, &
    - 4.15d0, &
    6.11d0, &
    - 2.09d0, &
    - 3.64d0, &
    11.07d0 &
    /)

  ! The longitude of the ascending node of the planet's orbit (in deg).
  ! (Note that the ascending node for the Earth is usually given as
  ! - 11.26064 deg, but is given as 348.73936 deg to ensure that
  ! 0 <= ascnode < 2 pi is satisfied.)
  real(dp), parameter :: SSASCNODE_PL_DEG(SSNPL) &
    = (/ &
    48.33167d0, &
    76.68069d0, &
    348.73936d0, &
    49.57854d0, &
    100.55615d0, &
    113.71504d0, &
    74.22988d0, &
    131.72169d0, &
    110.30347d0 &
    /)

  ! The rate of change of the longitude of the planet's orbit
  ! [in arcsec cent^(-1)].
  real(dp), parameter :: SSDASCNODEDT_PL_ARCSECCENT(SSNPL) &
    = (/ &
    - 446.30d0, &
    - 999.89d0, &
    - 18228.25d0, &
    - 1020.19d0, &
    1217.17d0, &
    - 1591.05d0, &
    1681.40d0, &
    - 151.25d0, &
    - 37.33d0 &
    /)

  ! The longitude of pericentre of the planet's orbit (in deg).
  real(dp), parameter :: SSLONPERI_PL_DEG(SSNPL) &
    = (/ &
    77.45645d0, &
    131.53298d0, &
    102.94719d0, &
    336.04084d0, &
    14.75385d0, &
    92.43194d0, &
    170.96424d0, &
    44.97135d0, &
    224.06676d0 &
    /)

  ! The rate of change of the longitude of pericentre of the planet's
  ! orbit [in arcsec cent^(-1)].
  real(dp), parameter :: SSDLONPERIDT_PL_ARCSECCENT(SSNPL) &
    = (/ &
    573.57d0, &
    - 108.80d0, &
    1198.28d0, &
    1560.78d0, &
    839.93d0, &
    - 1948.89d0, &
    1312.56d0, &
    - 844.43d0, &
    - 132.25d0 &
    /)

  ! The mean longitude of the planet's orbit (in deg).
  real(dp), parameter :: SSMEANLON_PL_DEG(SSNPL) &
    = (/ &
    252.25084d0, &
    181.97973d0, &
    100.46435d0, &
    355.45332d0, &
    34.40438d0, &
    49.94432d0, &
    313.23218d0, &
    304.88003d0, &
    238.92881d0 &
    /)

  ! The rate of change of the mean longitude of the planet's orbit,
  ! expressed in two parts, due to the fact that this will have
  ! cycled through many times (for some planets) during the century
  ! over which these rates are measured. The first part is the number
  ! of full revolutions per century, the second the number of additional
  ! arcsec. Thus the full angle through which the mean longitude will
  ! have changed in a century is given by (in radians):
  !
  !   ssdmeanlondt_pl_radcent &
  !     = twopi * SSDMEANLONDT_PL_REVCENT &
  !     + GNARCSEC_RAD * SSDMEANLONDT_PL_ARCSECCENT
  real(dp), parameter :: SSDMEANLONDT_PL_REVCENT(SSNPL) &
    = (/ &
    415.0d0, &
    162.0d0, &
    99.0d0, &
    53.0d0, &
    8.0d0, &
    3.0d0, &
    1.0d0, &
    0.0d0, &
    0.0d0 &
    /)

  real(dp), parameter :: SSDMEANLONDT_PL_ARCSECCENT(SSNPL) &
    = (/ &
    261628.29d0, &
    712136.06d0, &
    1293740.63d0, &
    217103.78d0, &
    557078.35d0, &
    513052.95d0, &
    246547.79d0, &
    786449.21d0, &
    522747.90d0 &
    /)

  ! Constants that give the ecliptic longitude of the Sun.
  real(dp), parameter :: &
    SSSUN_L0_DEG = 280.264_dp, &
    SSSUN_L1_DEG = 36000.770_dp, &
    SSSUN_G0_DEG = 357.528_dp, &
    SSSUN_G1_DEG = 35999.050_dp, &
    SSSUN_F1_DEG = 1.915_dp, &
    SSSUN_F2_DEG = 0.020_dp

  ! Return the index number from the name of a planet.
  interface ss_pl_name
    module procedure ss_pl_name_main, ss_pl_name_solsys
  end interface

  ! Initialises a planet structure.
  interface ss_planet_init
    module procedure ss_planet_init_name, ss_planet_init_pl, ss_planet_init_0
  end interface

contains

  ! Returns the planet of the given name from the Solar system structure.
  function ss_planet_name(name, solsys) result(planet)
    character(len=*), intent(in) :: name
    type(sssolsys), intent(in) :: solsys
    type(ssplanet) :: planet

    integer :: pl

    pl = ss_pl_name(name, solsys)
    call gn_assert (pl>=0, 'ss_planet_name: pl < 0', pl)
    call gn_assert (pl<=solsys%npl, 'ss_planet_name: pl > n_pl', pl)
    planet = solsys%planets(pl)
  end function ss_planet_name

  ! Returns the index pl from the parameter array of main planets
  ! corresponding to the input name.
  function ss_pl_name_main(name) result(pl)
    character(len=*), intent(in) :: name
    integer :: pl

    integer :: pl_temp
    logical :: found

    ! Search through the array of names.
    pl = 0
    found = .false.
    do pl_temp = 1, SSNPL
      if (name == SSNAME_PL(pl_temp)) then
        found = .true.
        pl = pl_temp
      end if
    end do

    ! If the name doesn't correspond to any in the array of main planets
    ! there is a serious error.
    if (.not. found) call gn_warning('ss_pl_name_main: no main planet', name)
  end function ss_pl_name_main

  ! Returns the index pl from the array of planets currently stored in
  ! solsys corresponding to the input name.
  function ss_pl_name_solsys(name, solsys) result(pl)
    character(len=*), intent(in) :: name
    type(sssolsys), intent(in) :: solsys
    integer :: pl

    integer :: pl_temp
    logical :: found

    ! Search through the array of names in the Solar system structure.
    pl = 0
    found = .false.
    do pl_temp = 1, SSNPL
      if (name == solsys%planets(pl_temp)%name) then
        found = .true.
        pl = pl_temp
      end if
    end do

    ! If the name doesn't correspond to any in the array of main planets
    ! there is a serious error.
    if (.not. found) &
      call gn_warning('ss_pl_name_main: no planet in solsys', name)
  end function ss_pl_name_solsys

  ! Initialize the Solar system planets, adjusting their orbits to the
  ! best fit values at J2000.
  function ss_solsys_init() result(solsys)
    type(sssolsys) :: solsys

    integer :: pl
    type(ssplanet) :: planet
    type(gnsec) :: t_0_local

    ! Use the default time of J2000.
    t_0_local = gn_cent2sec(SST_0_CENT)

    write(*, '(/,a,f6.1,a,/,a,i10,a,/)') &
      'Initialising Solar system with reference year ', &
      gn_sec2year(t_0_local), '.', '  (', t_0_local%int, &
      ' seconds since 1970.0.)'

    ! Initially set the number of planets to zero.
    solsys = ss_solsys_init_0()

    ! Initialise the central star.
    solsys%star = ss_star_init(SSM_SUN_KG, SSR_SUN_M, SSRADFORCE_SUN_N)

    ! Then add the main planets in turn, automatically handling the
    ! memory requirements.
    do pl = 1, SSNPL
      planet = ss_planet_init(SSNAME_PL(pl), t_0_local)
      write(*, '(a, i2, a, a, a)') 'Adding planet ', pl, ' (', &
        trim(SSNAME_PL(pl)), ') to Solar system ...'
      call ss_addplanet(solsys, planet)
    end do
  end function ss_solsys_init

  ! Creates the Solar system structure with no memory allocated and no
  ! planets.
  function ss_solsys_init_0() result(solsys)
    type(sssolsys) :: solsys

    solsys%npl = 0
    solsys%planets => null()
  end function ss_solsys_init_0

  ! Free the Solar system structure.
  subroutine ss_solsys_free(solsys)
    type(sssolsys), intent(in out) :: solsys

    ! Free the array of planets and make the pointer null.
    if (associated(solsys%planets)) deallocate(solsys%planets)
    solsys%planets => null()

    ! Set the number of planets to zero.
    solsys%npl = 0
  end subroutine ss_solsys_free

  ! Add a planet to the Solar system, allocating memory as required. The
  ! planet being added must be fully initialised.
  subroutine ss_addplanet(solsys, planet)
    type(sssolsys), intent(in out) :: solsys
    type(ssplanet), intent(in) :: planet

    type(sssolsys) :: solsys_temp
    logical :: matched
    integer :: pl

    call gn_assert(solsys%npl>=0,'ss_addplanet: npl < 0', solsys%npl)
    if (solsys%npl == 0) then

      ! With no planets in the Solar system at present just allocate
      ! memory and add in a new one.
      solsys%npl = 1
      allocate(solsys%planets(solsys%npl))

      solsys%planets(solsys%npl) = planet

    else

      ! Go through the current list of planets and check that the new
      ! planet doesn't share a name with any of them; if this is the
      ! case then don't add the new planet.
      matched = .false.
      do pl = 1, solsys%npl
        if (planet%name == solsys%planets(pl)%name) then
          matched = .true.
          call gn_warning('ss_addplanet: new planet has already used name', &
            planet%name)
        end if
      end do

      if (.not. matched) then

        ! Copy the current Solar system structure over to the temporary copy.
        solsys_temp%npl = solsys%npl
        allocate(solsys_temp%planets(solsys_temp%npl))

        solsys_temp%star = solsys%star

        do pl = 1, solsys%npl
          solsys_temp%planets(pl) = solsys%planets(pl)
        end do

        ! Then reallocate the memory to fit the new planet.
        call ss_solsys_free(solsys)

        solsys%npl = solsys_temp%npl + 1
        allocate(solsys%planets(solsys%npl))

        ! Copy the central star back into the structure.
        solsys%star = solsys_temp%star

        ! Copy the original planets back onto the newly allocated structure.
        do pl = 1, solsys_temp%npl
          solsys%planets(pl) = solsys_temp%planets(pl)
        end do

        ! Finally add in the new planet.
        solsys%planets(solsys%npl) = planet

        ! Deallocate the temporary array.
        call ss_solsys_free(solsys_temp)

      end if

    end if
  end subroutine ss_addplanet

  ! Initialise a major planet given its name.
  function ss_planet_init_name(name, t_0) result(planet)
    character(len=*), intent(in) :: name
    type(gnsec), intent(in) :: t_0
    type(ssplanet) :: planet

    integer :: pl

    ! Simply find the planet's index number and call the relevant routine.
    pl = ss_pl_name(name)
    planet = ss_planet_init(pl, t_0)
  end function ss_planet_init_name

  ! Initialise a major planet from its index.
  function ss_planet_init_pl(pl, t_0) result(planet)
    integer, intent(in) :: pl
    type(gnsec), intent(in) :: t_0
    type(ssplanet) :: planet

    character(len=filenamelen) :: name
    real(dp) :: delta_t_cent, m, r, a, e, inc, ascnode, lonperi, &
      meanlon, m_star

    call gn_assert(pl>=0,'ss_planet_init_pl: pl < 0', pl)
    call gn_assert(pl<=SSNPL,'ss_planet_init_pl: pl > SSNPL', pl)

    ! Convert from the reference time as input to a relative reference
    ! time (i.e., the time at which the orbital elements, etc. are
    ! correct).
    delta_t_cent = gn_sec2cent(t_0) - SST_0_CENT

    ! Extract the basic planetary data from the constant arrays in this
    ! module, adjusting the evolving quantities to their best fit
    ! values at time t_0. Also, the tabulated quantities are converted
    ! from their almanac values to the package's default units: radians
    ! for angles, metres for lengths, kilograms for masses.
    name = SSNAME_PL(pl)

    m = SSM_PL_KG(pl)
    r = SSR_PL_M(pl)

    a = GNAU_M * (SSA_PL_AU(pl) + delta_t_cent * SSDADT_PL_AUCENT(pl))
    e = SSE_PL(pl) + delta_t_cent * SSDEDT_PL_CENT(pl)
    inc = GNDEG_RAD * (SSINC_PL_DEG(pl) &
      + delta_t_cent * GNARCSEC_DEG * SSDINCDT_PL_ARCSECCENT(pl))
    ascnode = GNDEG_RAD * (SSASCNODE_PL_DEG(pl) &
      + delta_t_cent * GNARCSEC_DEG * SSDASCNODEDT_PL_ARCSECCENT(pl))
    lonperi = GNDEG_RAD * (SSLONPERI_PL_DEG(pl) &
      + delta_t_cent * GNARCSEC_DEG * SSDLONPERIDT_PL_ARCSECCENT(pl))
    meanlon = GNDEG_RAD * (SSMEANLON_PL_DEG(pl) &
      + delta_t_cent * GNARCSEC_DEG * SSDMEANLONDT_PL_ARCSECCENT(pl) &
      + delta_t_cent * 360.0 * SSDMEANLONDT_PL_REVCENT(pl))

    ! The mass of the central star of the Solar system.
    m_star = SSM_SUN_KG

    ! Finally, initialise the planet using the above quantities.
    planet = ss_planet_init(name, m, r, t_0, a, e, inc, ascnode, lonperi, &
      meanlon, m_star)
  end function ss_planet_init_pl

  ! Approximate the ecliptic longitude of the central star, in radians, at
  ! time t.
  function ss_lambda_star(t) result(lambda_star)
    type(gnsec), intent(in) :: t
    real(dp) :: lambda_star

    real(dp) :: t_cent, t_rel_cent, lambda, gamma

    ! Time conversion to centuries since J2000.
    t_cent = gn_sec2cent(t)
    t_rel_cent = t_cent + 0.12_dp - GNJ2000_CENT

    ! Calculate the central star's ecliptic longitude.
    lambda = GNDEG_RAD * SSSUN_L0_DEG + GNDEG_RAD * SSSUN_L1_DEG * t_rel_cent
    gamma = GNDEG_RAD * SSSUN_G0_DEG + GNDEG_RAD * SSSUN_G1_DEG * t_rel_cent
    gamma = mod(gamma, twopi)

    lambda_star = lambda + GNDEG_RAD * SSSUN_F1_DEG * sin(gamma) &
      + GNDEG_RAD * SSSUN_F2_DEG * sin(2.0_dp * gamma)

    ! For elegance convert lambda_star to be in the range 0 <= lambda_star
    ! <= 2 pi.
    lambda_star = mod(lambda_star, twopi)
    if (lambda_star < 0.0) lambda_star = lambda_star + twopi
  end function ss_lambda_star

end module solarsystem_solsys
