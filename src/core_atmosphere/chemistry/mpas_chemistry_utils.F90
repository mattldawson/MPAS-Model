! Copyright (C) 2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Chemistry utility functions: solar geometry for TUV-x
!
module mpas_chemistry_utils

   use mpas_kind_types, only : RKIND
   use iso_fortran_env, only : real64

   implicit none

   private
   public :: compute_solar_zenith_angle, compute_earth_sun_distance

   ! Universal physical constants (species-independent)
   real (kind=real64), parameter, public :: AVOGADRO = 6.02214076e23_real64
   real (kind=real64), parameter, public :: MW_AIR = 0.029_real64          ! [kg/mol]
   real (kind=real64), parameter, public :: SCALE_HEIGHT_AIR = 8.01_real64 ! [km]

   real (kind=RKIND), parameter :: pi = 3.14159265358979323846_RKIND
   real (kind=RKIND), parameter :: deg2rad = pi / 180.0_RKIND

contains

   !> Compute solar zenith angle from latitude, longitude, and time
   !! Uses standard solar geometry (Meeus, Astronomical Algorithms)
   !!
   !! @param[in] lat_rad    Latitude in radians
   !! @param[in] lon_rad    Longitude in radians
   !! @param[in] julday     Julian day of year (1-366)
   !! @param[in] ut_hours   UTC hour of day (0-24)
   !! @return    SZA in radians
   function compute_solar_zenith_angle(lat_rad, lon_rad, julday, ut_hours) result(sza)
      real (kind=RKIND), intent(in) :: lat_rad
      real (kind=RKIND), intent(in) :: lon_rad
      integer,           intent(in) :: julday
      real (kind=RKIND), intent(in) :: ut_hours
      real (kind=RKIND) :: sza

      real (kind=RKIND) :: decl, hour_angle, cos_sza
      real (kind=RKIND) :: gamma

      ! Fractional year in radians
      gamma = 2.0_RKIND * pi * real(julday - 1, RKIND) / 365.0_RKIND

      ! Solar declination (Spencer, 1971)
      decl = 0.006918_RKIND &
           - 0.399912_RKIND * cos(gamma) &
           + 0.070257_RKIND * sin(gamma) &
           - 0.006758_RKIND * cos(2.0_RKIND * gamma) &
           + 0.000907_RKIND * sin(2.0_RKIND * gamma) &
           - 0.002697_RKIND * cos(3.0_RKIND * gamma) &
           + 0.00148_RKIND  * sin(3.0_RKIND * gamma)

      ! Hour angle: solar noon at lon=0 when ut_hours=12
      hour_angle = (ut_hours - 12.0_RKIND) * (pi / 12.0_RKIND) + lon_rad

      ! Cosine of SZA
      cos_sza = sin(lat_rad) * sin(decl) &
              + cos(lat_rad) * cos(decl) * cos(hour_angle)

      ! Clamp and convert
      cos_sza = max(-1.0_RKIND, min(1.0_RKIND, cos_sza))
      sza = acos(cos_sza)

   end function compute_solar_zenith_angle

   !> Compute Earth-Sun distance in AU from Julian day
   !! Uses a simple approximation (Spencer, 1971)
   !!
   !! @param[in] julday  Julian day of year (1-366)
   !! @return    Earth-Sun distance in AU
   function compute_earth_sun_distance(julday) result(dist_au)
      integer, intent(in) :: julday
      real (kind=RKIND) :: dist_au

      real (kind=RKIND) :: gamma

      gamma = 2.0_RKIND * pi * real(julday - 1, RKIND) / 365.0_RKIND

      ! Inverse square of distance ratio (Spencer, 1971)
      dist_au = 1.000110_RKIND &
              + 0.034221_RKIND * cos(gamma) &
              + 0.001280_RKIND * sin(gamma) &
              + 0.000719_RKIND * cos(2.0_RKIND * gamma) &
              + 0.000077_RKIND * sin(2.0_RKIND * gamma)

      ! dist_au is actually (r0/r)^2, so r/r0 = 1/sqrt(dist_au)
      dist_au = 1.0_RKIND / sqrt(dist_au)

   end function compute_earth_sun_distance

end module mpas_chemistry_utils
