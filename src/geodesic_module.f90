!*****************************************************************************************
!>
!  Implementation of geodesic routines in modern Fortran
!
!### See also
!  * This module contains modernized versions of the Fortran
!    code from [geographiclib-fortran](https://github.com/geographiclib/geographiclib-fortran).
!    The `geographiclib` subroutines are documented at:
!    [sourceforge](https://geographiclib.sourceforge.io/html/Fortran/)
!    These are Fortran implementation of the geodesic algorithms described in:
!    C. F. F. Karney, [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z),
!    J. Geodesy 87, 43--55 (2013). [Apr 23, 2022 version]
!  * Some of the code was also split off from the
!    [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit)
!    (see `geodesy_module.f90`). [Nov 6, 2022 version]
!
!### geographiclib Notes
! The principal advantages of these algorithms over previous ones
! (e.g., Vincenty, 1975) are
! - accurate to round off for |<i>f</i>| &lt; 1/50;
! - the solution of the inverse problem is always found;
! - differential and integral properties of geodesics are computed.
!
! The shortest path between two points on the ellipsoid at (\e lat1, \e
! lon1) and (\e lat2, \e lon2) is called the geodesic.  Its length is
! \e s12 and the geodesic from point 1 to point 2 has forward azimuths
! \e azi1 and \e azi2 at the two end points.
!
! Traditionally two geodesic problems are considered:
! - the direct problem -- given \e lat1, \e lon1, \e s12, and \e azi1,
!   determine \e lat2, \e lon2, and \e azi2.  This is solved by the
!   subroutine direct().
! - the inverse problem -- given \e lat1, \e lon1, \e lat2, \e lon2,
!   determine \e s12, \e azi1, and \e azi2.  This is solved by the
!   subroutine inverse().
!
! The ellipsoid is specified by its equatorial radius \e a (typically
! in meters) and flattening \e f.  The routines are accurate to round
! off with real(wp) arithmetic provided that |<i>f</i>| &lt;
! 1/50; for the WGS84 ellipsoid, the errors are less than 15
! nanometers.  (Reasonably accurate results are obtained for |<i>f</i>|
! &lt; 1/5.)  For a prolate ellipsoid, specify \e f &lt; 0.
!
! The routines also calculate several other quantities of interest
! - \e SS12 is the area between the geodesic from point 1 to point 2
!   and the equator; i.e., it is the area, measured counter-clockwise,
!   of the geodesic quadrilateral with corners (\e lat1,\e lon1), (0,\e
!   lon1), (0,\e lon2), and (\e lat2,\e lon2).
! - \e m12, the reduced length of the geodesic is defined such that if
!   the initial azimuth is perturbed by \e dazi1 (radians) then the
!   second point is displaced by \e m12 \e dazi1 in the direction
!   perpendicular to the geodesic.  On a curved surface the reduced
!   length obeys a symmetry relation, \e m12 + \e m21 = 0.  On a flat
!   surface, we have \e m12 = \e s12.
! - \e MM12 and \e MM21 are geodesic scales.  If two geodesics are
!   parallel at point 1 and separated by a small distance \e dt, then
!   they are separated by a distance \e MM12 \e dt at point 2.  \e MM21
!   is defined similarly (with the geodesics being parallel to one
!   another at point 2).  On a flat surface, we have \e MM12 = \e MM21
!   = 1.
! - \e a12 is the arc length on the auxiliary sphere.  This is a
!   construct for converting the problem to one in spherical
!   trigonometry.  \e a12 is measured in degrees.  The spherical arc
!   length from one equator crossing to the next is always 180 deg.
!
! If points 1, 2, and 3 lie on a single geodesic, then the following
! addition rules hold:
! - \e s13 = \e s12 + \e s23
! - \e a13 = \e a12 + \e a23
! - \e SS13 = \e SS12 + \e SS23
! - \e m13 = \e m12 \e MM23 + \e m23 \e MM21
! - \e MM13 = \e MM12 \e MM23 - (1 - \e MM12 \e MM21) \e
!   m23 / \e m12
! - \e MM31 = \e MM32 \e MM21 - (1 - \e MM23 \e MM32) \e
!   m12 / \e m23
!
! The shortest distance returned by the solution of the inverse problem
! is (obviously) uniquely defined.  However, in a few special cases
! there are multiple azimuths which yield the same shortest distance.
! Here is a catalog of those cases:
! - \e lat1 = -\e lat2 (with neither point at a pole).  If \e
!   azi1 = \e azi2, the geodesic is unique.  Otherwise there are two
!   geodesics and the second one is obtained by setting [\e azi1, \e
!   azi2] &rarr; [\e azi2, \e azi1], [\e MM12, \e MM21] &rarr; [\e
!   MM21, \e MM12], \e SS12 &rarr; -\e SS12.  (This occurs when
!   the longitude difference is near &plusmn;180 deg for oblate
!   ellipsoids.)
! - \e lon2 = \e lon1 &plusmn; 180 deg (with neither point at a pole).
!   If \e azi1 = 0 deg or &plusmn;180 deg, the geodesic is unique.
!   Otherwise there are two geodesics and the second one is obtained by
!   setting [\e azi1, \e azi2] &rarr; [-\e azi1, -\e azi2],
!   \e SS12 &rarr; -\e SS12.  (This occurs when \e lat2 is near
!   -\e lat1 for prolate ellipsoids.)
! - Points 1 and 2 at opposite poles.  There are infinitely many
!   geodesics which can be generated by setting [\e azi1, \e azi2]
!   &rarr; [\e azi1, \e azi2] + [\e d, -\e d], for arbitrary \e
!   d.  (For spheres, this prescription applies when points 1 and 2 are
!   antipodal.)
! - \e s12 = 0 (coincident points).  There are infinitely many
!   geodesics which can be generated by setting [\e azi1, \e azi2]
!   &rarr; [\e azi1, \e azi2] + [\e d, \e d], for arbitrary \e d.
!
! These routines are a simple transcription of the corresponding C++
! classes in <a href="https://geodesic_module.sourceforge.io">
! geodesic_module</a>.  Because of the limitations of Fortran 77, the
! classes have been replaced by simple subroutines with no attempt to
! save "state" across subroutine calls.  Most of the internal comments
! have been retained.  However, in the process of transcription some
! documentation has been lost and the documentation for the C++
! classes, geodesic_module::Geodesic, geodesic_module::GeodesicLine, and
! geodesic_module::PolygonAreaT, should be consulted.  The C++ code
! remains the "reference implementation".  Think twice about
! restructuring the internals of the Fortran code since this may make
! porting fixes from the C++ code more difficult.
!
!### License
!
!  `geographiclib-fortran` Copyright (c) Charles Karney (2012-2022) <charles@karney.com> and
!  licensed under the MIT/X11 License.  For more information, see
!  https://geodesic_module.sourceforge.io/

module geodesic_module

  use iso_fortran_env, only: wp => real64

  implicit none

  private

  real(wp),parameter :: zero = 0.0_wp
  real(wp),parameter :: one = 1.0_wp
  real(wp),parameter :: two = 2.0_wp
  real(wp),parameter :: three = 3.0_wp

  integer,parameter :: maxit1 = 20
  integer,parameter :: maxit2 = maxit1 + digits(one) + 10

  real(wp),parameter :: dblmin = tiny(one) !! 0.5_wp**1022
  real(wp),parameter :: dbleps = epsilon(one) !! 0.5_wp**(digits-1)

  real(wp),parameter :: pi = atan2(zero, -one)
  real(wp),parameter :: degree = pi/180.0_wp
  real(wp),parameter :: twopi = two * pi
  real(wp),parameter :: halfpi = pi / two

  ! This is about cbrt(dblmin).  With other implementations, sqrt(dblmin)
  ! is used.  The larger value is used here to avoid complaints about a
  ! IEEE_UNDERFLOW_FLAG IEEE_DENORMAL signal.  This is triggered when
  ! inverse is called with points at opposite poles.
  real(wp),parameter :: tiny2 = dblmin**(one/three) !! 0.5_wp**((1022+1)/3)
  real(wp),parameter :: tol0 = dbleps

  ! Increase multiplier in defn of tol1 from 100 to 200 to fix inverse
  ! case 52.784459512564 0 -52.784459512563990912 179.634407464943777557
  ! which otherwise failed for Visual Studio 10 (Release and Debug)
  real(wp),parameter :: tol1 = 200.0_wp * tol0
  real(wp),parameter :: tol2 = sqrt(tol0)

  ! Check on bisection interval
  real(wp),parameter :: tolb = tol0 * tol2
  real(wp),parameter :: xthresh = 1000.0_wp * tol2

  public :: direct
  public :: inverse
  public :: direct_vincenty
  public :: inverse_vincenty
  public :: area
  public :: heikkinen
  public :: olson
  public :: cartesian_to_geodetic_triaxial
  public :: CartesianIntoGeodeticI, CartesianIntoGeodeticII
  public :: cartesian_to_geodetic_triaxial_2
  public :: geodetic_to_cartesian
  public :: geodetic_to_cartesian_triaxial
  public :: geodetic_to_cartesian_triaxial_2
  public :: great_circle_distance
  public :: geocentric_radius

  ! other routines
  public :: AngDiff,AngNormalize,sumx,LatFix,atan2d,AngRound,sincosd

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve the direct geodesic problem.
!
!### Notes
!  If `arcmode` is not set, `s12a12` is `s12` and `a12s12` is
!  `a12`; otherwise, `s12a12` is `a12` and `a12s12` is `s12`.  If
!  `unroll` is not set, the value `lon2` returned is in the range
!  [-180 deg, 180 deg]; if `unroll` is set, the longitude variable
!  is "unrolled" so that `lon2 - lon1` indicates how many
!  times and in what sense the geodesic encircles the ellipsoid.
!
!  If either point is at a pole, the azimuth is defined by keeping the
!  longitude fixed, writing `lat = lat = +/- (90 deg - Epsilon)`,
!  and taking the limit Epsilon --> 0+.  An arc length
!  greater that 180 deg signifies a geodesic which is not a shortest
!  path.  (For a prolate ellipsoid, an additional condition is necessary
!  for a shortest path: the longitudinal extent must not exceed of
!  180 deg.)

subroutine direct(a, f, lat1, lon1, azi1, s12a12, flags, &
                  lat2, lon2, azi2, outmask, a12s12, m12, MM12, MM21, SS12)

  real(wp), intent(in) :: a !! the equatorial radius (meters).
  real(wp), intent(in) :: f !! the flattening of the ellipsoid.  Setting `f = 0` gives
                                    !! a sphere.  Negative `f` gives a prolate ellipsoid.
  real(wp), intent(in) :: lat1 !! lat1 latitude of point 1 (degrees). `lat1` should be in the range [-90 deg, 90 deg].
  real(wp), intent(in) :: lon1 !! lon1 longitude of point 1 (degrees).
  real(wp), intent(in) :: azi1 !! azi1 azimuth at point 1 (degrees).
  real(wp), intent(in) :: s12a12 !! if `arcmode` is not set, this is the distance
                                         !! from point 1 to point 2 (meters); otherwise it is the arc
                                         !! length from point 1 to point 2 (degrees); it can be negative.
  integer, intent(in) :: flags !! a bitor'ed combination of the `arcmode` and `unroll` flags.
                               !!
                               !! `flags` is an integer in [0, 4) whose binary bits are interpreted as follows:
                               !!
                               !!  * 1 the `arcmode` flag
                               !!  * 2 the `unroll` flag
  integer, intent(in) :: outmask !! a bitor'ed combination of mask values
                               !! specifying which of the following parameters should be set.
                               !!
                               !! outmask is an integer in [0, 16) whose binary bits are interpreted
                               !! as follows:
                               !!
                               !!  * 1 return `a12`
                               !!  * 2 return `m12`
                               !!  * 4 return `MM12` and `MM21`
                               !!  * 8 return `SS12`
  real(wp), intent(out) :: lat2 !! latitude of point 2 (degrees).
  real(wp), intent(out) :: lon2 !! longitude of point 2 (degrees).
  real(wp), intent(out) :: azi2 !! (forward) azimuth at point 2 (degrees). The value `azi2` returned is in the range [-180 deg, 180 deg].
  real(wp), intent(out) :: a12s12 !! if `arcmode` is not set, this is the arc length
                                          !! from point 1 to point 2 (degrees); otherwise it is the distance
                                          !! from point 1 to point 2 (meters).
  real(wp), intent(out) :: m12 !! reduced length of geodesic (meters).
  real(wp), intent(out) :: MM12 !! geodesic scale of point 2 relative to point 1 (dimensionless).
  real(wp), intent(out) :: MM21 !! geodesic scale of point 1 relative to point 2 (dimensionless).
  real(wp), intent(out) :: SS12 !! area under the geodesic (\(m^2\)).

  integer,parameter :: ord = 6
  integer,parameter :: nC1 = ord
  integer,parameter :: nC1p = ord
  integer,parameter :: nC2 = ord
  integer,parameter :: nA3 = ord
  integer,parameter :: nA3x = nA3
  integer,parameter :: nC3 = ord
  integer,parameter :: nC3x = (nC3 * (nC3 - 1)) / 2
  integer,parameter :: nC4 = ord
  integer,parameter :: nC4x = (nC4 * (nC4 + 1)) / 2

  real(wp) :: A3x(0:nA3x-1), C3x(0:nC3x-1), C4x(0:nC4x-1), &
              C1a(nC1), C1pa(nC1p), C2a(nC2), C3a(nC3-1), C4a(0:nC4-1)

  logical :: arcmode, unroll, arcp, redlp, scalp, areap
  real(wp) :: e2, f1, ep2, n, b, c2, &
              salp0, calp0, k2, eps, &
              salp1, calp1, ssig1, csig1, cbet1, sbet1, dn1, somg1, comg1, &
              salp2, calp2, ssig2, csig2, sbet2, cbet2, dn2, somg2, comg2, &
              ssig12, csig12, salp12, calp12, omg12, lam12, lon12, &
              sig12, stau1, ctau1, tau12, t, s, c, serr, E, &
              A1m1, A2m1, A3c, A4, AB1, AB2, &
              B11, B12, B21, B22, B31, B41, B42, J12

e2 = f * (2 - f)
ep2 = e2 / (1 - e2)
f1 = 1 - f
n = f / (2 - f)
b = a * f1
c2 = 0

arcmode = mod(flags/1, 2) == 1
unroll = mod(flags/2, 2) == 1

arcp = mod(outmask/1, 2) == 1
redlp = mod(outmask/2, 2) == 1
scalp = mod(outmask/4, 2) == 1
areap = mod(outmask/8, 2) == 1

if (areap) then
  if (e2 == 0) then
    c2 = a**2
  else if (e2 > 0) then
    c2 = (a**2 + b**2 * atanh(sqrt(e2)) / sqrt(e2)) / 2
  else
    c2 = (a**2 + b**2 * atan(sqrt(abs(e2))) / sqrt(abs(e2))) / 2
  end if
end if

call A3coeff(n, A3x)
call C3coeff(n, C3x)
if (areap) call C4coeff(n, C4x)

! Guard against underflow in salp0
call sincosd(AngRound(azi1), salp1, calp1)

call sincosd(AngRound(LatFix(lat1)), sbet1, cbet1)
sbet1 = f1 * sbet1
call norm(sbet1, cbet1)
! Ensure cbet1 = +dbleps at poles
cbet1 = max(tiny2, cbet1)
dn1 = sqrt(1 + ep2 * sbet1**2)

! Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
! alp0 in [0, pi/2 - |bet1|]
salp0 = salp1 * cbet1
! Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
! is slightly better (consider the case salp1 = 0).
calp0 = hypot(calp1, salp1 * sbet1)
! Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
! sig = 0 is nearest northward crossing of equator.
! With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
! With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
! With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
! Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
! With alp0 in (0, pi/2], quadrants for sig and omg coincide.
! No atan2(0,0) ambiguity at poles since cbet1 = +dbleps.
! With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
ssig1 = sbet1
somg1 = salp0 * sbet1
if (sbet1 /= 0 .or. calp1 /= 0) then
  csig1 = cbet1 * calp1
else
  csig1 = 1
end if
comg1 = csig1
! sig1 in (-pi, pi]
call norm(ssig1, csig1)
! norm(somg1, comg1); -- don't need to normalize!

k2 = calp0**2 * ep2
eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2)

A1m1 = A1m1f(eps)
call C1f(eps, C1a)
B11 = SinCosSeries(.true., ssig1, csig1, C1a, nC1)
s = sin(B11)
c = cos(B11)
! tau1 = sig1 + B11
stau1 = ssig1 * c + csig1 * s
ctau1 = csig1 * c - ssig1 * s
! Not necessary because C1pa reverts C1a
!    B11 = -SinCosSeries(true, stau1, ctau1, C1pa, nC1p)

if (.not. arcmode) call C1pf(eps, C1pa)

if (redlp .or. scalp) then
  A2m1 = A2m1f(eps)
  call C2f(eps, C2a)
  B21 = SinCosSeries(.true., ssig1, csig1, C2a, nC2)
else
! Suppress bogus warnings about unitialized variables
  A2m1 = 0
  B21 = 0
end if

call C3f(eps, C3x, C3a)
A3c = -f * salp0 * A3f(eps, A3x)
B31 = SinCosSeries(.true., ssig1, csig1, C3a, nC3-1)

if (areap) then
  call C4f(eps, C4x, C4a)
! Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
  A4 = a**2 * calp0 * salp0 * e2
  B41 = SinCosSeries(.false., ssig1, csig1, C4a, nC4)
else
! Suppress bogus warnings about unitialized variables
  A4 = 0
  B41 = 0
end if

if (arcmode) then
! Interpret s12a12 as spherical arc length
  sig12 = s12a12 * degree
  call sincosd(s12a12, ssig12, csig12)
! Suppress bogus warnings about unitialized variables
  B12 = 0
else
! Interpret s12a12 as distance
  tau12 = s12a12 / (b * (1 + A1m1))
  s = sin(tau12)
  c = cos(tau12)
! tau2 = tau1 + tau12
  B12 = - SinCosSeries(.true., &
      stau1 * c + ctau1 * s, ctau1 * c - stau1 * s, C1pa, nC1p)
  sig12 = tau12 - (B12 - B11)
  ssig12 = sin(sig12)
  csig12 = cos(sig12)
  if (abs(f) > 0.01_wp) then
! Reverted distance series is inaccurate for |f| > 1/100, so correct
! sig12 with 1 Newton iteration.  The following table shows the
! approximate maximum error for a = WGS_a() and various f relative to
! GeodesicExact.
!     erri = the error in the inverse solution (nm)
!     errd = the error in the direct solution (series only) (nm)
!     errda = the error in the direct solution (series + 1 Newton) (nm)
!
!       f     erri  errd errda
!     -1/5    12e6 1.2e9  69e6
!     -1/10  123e3  12e6 765e3
!     -1/20   1110 108e3  7155
!     -1/50  18.63 200.9 27.12
!     -1/100 18.63 23.78 23.37
!     -1/150 18.63 21.05 20.26
!      1/150 22.35 24.73 25.83
!      1/100 22.35 25.03 25.31
!      1/50  29.80 231.9 30.44
!      1/20   5376 146e3  10e3
!      1/10  829e3  22e6 1.5e6
!      1/5   157e6 3.8e9 280e6
    ssig2 = ssig1 * csig12 + csig1 * ssig12
    csig2 = csig1 * csig12 - ssig1 * ssig12
    B12 = SinCosSeries(.true., ssig2, csig2, C1a, nC1)
    serr = (1 + A1m1) * (sig12 + (B12 - B11)) - s12a12 / b
    sig12 = sig12 - serr / sqrt(1 + k2 * ssig2**2)
    ssig12 = sin(sig12)
    csig12 = cos(sig12)
! Update B12 below
  end if
end if

! sig2 = sig1 + sig12
ssig2 = ssig1 * csig12 + csig1 * ssig12
csig2 = csig1 * csig12 - ssig1 * ssig12
dn2 = sqrt(1 + k2 * ssig2**2)
if (arcmode .or. abs(f) > 0.01_wp) &
    B12 = SinCosSeries(.true., ssig2, csig2, C1a, nC1)
AB1 = (1 + A1m1) * (B12 - B11)

! sin(bet2) = cos(alp0) * sin(sig2)
sbet2 = calp0 * ssig2
! Alt: cbet2 = hypot(csig2, salp0 * ssig2)
cbet2 = hypot(salp0, calp0 * csig2)
if (cbet2 == 0) then
! I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
  cbet2 = tiny2
  csig2 = cbet2
end if
! tan(omg2) = sin(alp0) * tan(sig2)
! No need to normalize
somg2 = salp0 * ssig2
comg2 = csig2
! tan(alp0) = cos(sig2)*tan(alp2)
! No need to normalize
salp2 = salp0
calp2 = calp0 * csig2
! East or west going?
E = sign(1.0_wp, salp0)
! omg12 = omg2 - omg1
if (unroll) then
  omg12 = E * (sig12 &
      - (atan2(    ssig2, csig2) - atan2(    ssig1, csig1)) &
      + (atan2(E * somg2, comg2) - atan2(E * somg1, comg1)))
else
  omg12 = atan2(somg2 * comg1 - comg2 * somg1, &
      comg2 * comg1 + somg2 * somg1)
end if

lam12 = omg12 + A3c * &
    ( sig12 + (SinCosSeries(.true., ssig2, csig2, C3a, nC3-1) &
    - B31))
lon12 = lam12 / degree
if (unroll) then
  lon2 = lon1 + lon12
else
  lon2 = AngNormalize(AngNormalize(lon1) + AngNormalize(lon12))
end if
lat2 = atan2d(sbet2, f1 * cbet2)
azi2 = atan2d(salp2, calp2)

if (redlp .or. scalp) then
  B22 = SinCosSeries(.true., ssig2, csig2, C2a, nC2)
  AB2 = (1 + A2m1) * (B22 - B21)
  J12 = (A1m1 - A2m1) * sig12 + (AB1 - AB2)
end if
! Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
! accurate cancellation in the case of coincident points.
if (redlp) m12 = b * ((dn2 * (csig1 * ssig2) - &
    dn1 * (ssig1 * csig2)) - csig1 * csig2 * J12)
if (scalp) then
  t = k2 * (ssig2 - ssig1) * (ssig2 + ssig1) / (dn1 + dn2)
  MM12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1
  MM21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2
end if

if (areap) then
  B42 = SinCosSeries(.false., ssig2, csig2, C4a, nC4)
  if (calp0 == 0 .or. salp0 == 0) then
! alp12 = alp2 - alp1, used in atan2 so no need to normalize
    salp12 = salp2 * calp1 - calp2 * salp1
    calp12 = calp2 * calp1 + salp2 * salp1
  else
! tan(alp) = tan(alp0) * sec(sig)
! tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
! = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
! If csig12 > 0, write
!   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
! else
!   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
! No need to normalize
    if (csig12 <= 0) then
      salp12 = csig1 * (1 - csig12) + ssig12 * ssig1
    else
      salp12 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
    end if
    salp12 = calp0 * salp0 * salp12
    calp12 = salp0**2 + calp0**2 * csig1 * csig2
  end if
  SS12 = c2 * atan2(salp12, calp12) + A4 * (B42 - B41)
end if

if (arcp) then
  if (arcmode) then
    a12s12 = b * ((1 + A1m1) * sig12 + AB1)
  else
    a12s12 = sig12 / degree
  end if
end if

end subroutine direct

!*****************************************************************************************
!>
!  Solve the inverse geodesic problem.
!
!  `lat1` and  `lat2` should be in the range [-90 deg, 90 deg].
!  The values of `azi1` and `azi2` returned are in the range
!  [-180 deg, 180 deg].
!
!## Notes
!  If either point is at a pole, the azimuth is defined by keeping the
!  longitude fixed, writing `lat = +/- (90 deg - Epsilon),
!  and taking the limit Epsilon --> 0+.
!
!  The solution to the inverse problem is found using Newton's method.
!  If this fails to converge (this is very unlikely in geodetic
!  applications but does occur for very eccentric ellipsoids), then the
!  bisection method is used to refine the solution.

subroutine inverse(a, f, lat1, lon1, lat2, lon2, &
                  s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)

real(wp), intent(in)  :: a !! the equatorial radius (meters).
real(wp), intent(in)  :: f !! the flattening of the ellipsoid.  Setting `f = 0` gives
                           !! a sphere.  Negative `f` gives a prolate ellipsoid.
real(wp), intent(in)  :: lat1 !! latitude of point 1 (degrees).
real(wp), intent(in)  :: lon1 !! longitude of point 1 (degrees).
real(wp), intent(in)  :: lat2 !! latitude of point 2 (degrees).
real(wp), intent(in)  :: lon2 !! longitude of point 2 (degrees).
integer, intent(in)   :: outmask !! a bitor'ed combination of mask values
                               !! specifying which of the following parameters should be set.
                               !! `outmask` is an integer in [0, 16) whose binary bits are interpreted
                               !! as follows:
                               !!
                               !!  *  1 return `a12`
                               !!  *  2 return `m12`
                               !!  *  4 return `MM12` and `MM21`
                               !!  *  8 return `SS12`
real(wp), intent(out) :: s12 !! distance from point 1 to point 2 (meters).
real(wp), intent(out) :: azi1 !! azimuth at point 1 (degrees).
real(wp), intent(out) :: azi2 !! (forward) azimuth at point 2 (degrees).
real(wp), intent(out) :: a12 !! arc length from point 1 to point 2 (degrees).
real(wp), intent(out) :: m12 !! reduced length of geodesic (meters).
real(wp), intent(out) :: MM12 !! MM12 geodesic scale of point 2 relative to point 1 (dimensionless).
real(wp), intent(out) :: MM21 !! MM21 geodesic scale of point 1 relative to point 2 (dimensionless).
real(wp), intent(out) :: SS12 !! SS12 area under the geodesic (\(m^2\)).

integer ord, nA3, nA3x, nC3, nC3x, nC4, nC4x, nC
parameter (ord = 6, nA3 = ord, nA3x = nA3, &
    nC3 = ord, nC3x = (nC3 * (nC3 - 1)) / 2, &
    nC4 = ord, nC4x = (nC4 * (nC4 + 1)) / 2, &
    nC = ord)
real(wp) A3x(0:nA3x-1), C3x(0:nC3x-1), C4x(0:nC4x-1), &
    Ca(nC)

integer latsign, lonsign, swapp, numit
logical arcp, redlp, scalp, areap, meridian, tripn, tripb

real(wp) e2, f1, ep2, n, b, c2, &
    lat1x, lat2x, salp0, calp0, k2, eps, &
    salp1, calp1, ssig1, csig1, cbet1, sbet1, dbet1, dn1, &
    salp2, calp2, ssig2, csig2, sbet2, cbet2, dbet2, dn2, &
    slam12, clam12, salp12, calp12, omg12, lam12, lon12, lon12s, &
    salp1a, calp1a, salp1b, calp1b, &
    dalp1, sdalp1, cdalp1, nsalp1, alp12, somg12, comg12, domg12, &
    sig12, v, dv, dnm, dummy, &
    A4, B41, B42, s12x, m12x, a12x, sdomg12, cdomg12

integer :: lmask

f1 = 1 - f
e2 = f * (2 - f)
ep2 = e2 / f1**2
n = f / ( 2 - f)
b = a * f1
c2 = 0

arcp = mod(outmask/1, 2) == 1
redlp = mod(outmask/2, 2) == 1
scalp = mod(outmask/4, 2) == 1
areap = mod(outmask/8, 2) == 1
if (scalp) then
  lmask = 16 + 2 + 4
else
  lmask = 16 + 2
end if

if (areap) then
  if (e2 == 0) then
    c2 = a**2
  else if (e2 > 0) then
    c2 = (a**2 + b**2 * atanh(sqrt(e2)) / sqrt(e2)) / 2
  else
    c2 = (a**2 + b**2 * atan(sqrt(abs(e2))) / sqrt(abs(e2))) / 2
  end if
end if

call A3coeff(n, A3x)
call C3coeff(n, C3x)
if (areap) call C4coeff(n, C4x)

! Compute longitude difference (AngDiff does this carefully).  Result is
! in [-180, 180] but -180 is only for west-going geodesics.  180 is for
! east-going and meridional geodesics.
! If very close to being on the same half-meridian, then make it so.
lon12 = AngDiff(lon1, lon2, lon12s)
! Make longitude difference positive.
lonsign = int(sign(1.0_wp, lon12))
lon12 = lonsign * lon12
lon12s = lonsign * lon12s
lam12 = lon12 * degree
! Calculate sincos of lon12 + error (this applies AngRound internally).
call sincosde(lon12, lon12s, slam12, clam12)
! the supplementary longitude difference
lon12s = (180 - lon12) - lon12s

! If really close to the equator, treat as on equator.
lat1x = AngRound(LatFix(lat1))
lat2x = AngRound(LatFix(lat2))
! Swap points so that point with higher (abs) latitude is point 1
! If one latitude is a nan, then it becomes lat1.
if (abs(lat1x) < abs(lat2x) .or. lat2x /= lat2x) then
  swapp = -1
else
  swapp = 1
end if
if (swapp < 0) then
  lonsign = -lonsign
  call swap(lat1x, lat2x)
end if
! Make lat1 <= 0
latsign = int(sign(1.0_wp, -lat1x))
lat1x = lat1x * latsign
lat2x = lat2x * latsign
! Now we have
!
!     0 <= lon12 <= 180
!     -90 <= lat1 <= 0
!     lat1 <= lat2 <= -lat1
!
! longsign, swapp, latsign register the transformation to bring the
! coordinates to this canonical form.  In all cases, 1 means no change
! was made.  We make these transformations so that there are few cases
! to check, e.g., on verifying quadrants in atan2.  In addition, this
! enforces some symmetries in the results returned.

call sincosd(lat1x, sbet1, cbet1)
sbet1 = f1 * sbet1
call norm(sbet1, cbet1)
! Ensure cbet1 = +dbleps at poles
cbet1 = max(tiny2, cbet1)

call sincosd(lat2x, sbet2, cbet2)
sbet2 = f1 * sbet2
call norm(sbet2, cbet2)
! Ensure cbet2 = +dbleps at poles
cbet2 = max(tiny2, cbet2)

! If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
! |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1
! is a better measure.  This logic is used in assigning calp2 in
! Lambda12.  Sometimes these quantities vanish and in that case we force
! bet2 = +/- bet1 exactly.  An example where is is necessary is the
! inverse problem 48.522876735459 0 -48.52287673545898293
! 179.599720456223079643 which failed with Visual Studio 10 (Release and
! Debug)

if (cbet1 < -sbet1) then
  if (cbet2 == cbet1) sbet2 = sign(sbet1, sbet2)
else
  if (abs(sbet2) == -sbet1) cbet2 = cbet1
end if

dn1 = sqrt(1 + ep2 * sbet1**2)
dn2 = sqrt(1 + ep2 * sbet2**2)

! Suppress bogus warnings about unitialized variables
a12x = 0
meridian = lat1x == -90 .or. slam12 == 0

if (meridian) then

! Endpoints are on a single full meridian, so the geodesic might lie on
! a meridian.

! Head to the target longitude
  calp1 = clam12
  salp1 = slam12
! At the target we're heading north
  calp2 = 1
  salp2 = 0

! tan(bet) = tan(sig) * cos(alp)
  ssig1 = sbet1
  csig1 = calp1 * cbet1
  ssig2 = sbet2
  csig2 = calp2 * cbet2

! sig12 = sig2 - sig1
  sig12 = atan2(max(0d0, csig1 * ssig2 - ssig1 * csig2) + 0d0, &
                         csig1 * csig2 + ssig1 * ssig2)
  call Lengths(n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, &
      cbet1, cbet2, lmask, &
      s12x, m12x, dummy, MM12, MM21, ep2, Ca)

! Add the check for sig12 since zero length geodesics might yield m12 <
! 0.  Test case was
!
!    echo 20.001 0 20.001 0 | GeodSolve -i
!
! In fact, we will have sig12 > pi/2 for meridional geodesic which is
! not a shortest path.
  if (sig12 < 1 .or. m12x >= 0) then
    if (sig12 < 3 * tiny2 .or. &
        (sig12 < tol0 .and. &
        (s12x < 0 .or. m12x < 0))) then
! Prevent negative s12 or m12 for short lines
      sig12 = 0
      m12x = 0
      s12x = 0
    end if
    m12x = m12x * b
    s12x = s12x * b
    a12x = sig12 / degree
  else
! m12 < 0, i.e., prolate and too close to anti-podal
    meridian = .false.
  end if
end if

omg12 = 0
! somg12 = 2 marks that it needs to be calculated
somg12 = 2
comg12 = 0
if (.not. meridian .and. sbet1 == 0 .and. &
    (f <= 0 .or. lon12s >= f * 180)) then

! Geodesic runs along equator
  calp1 = 0
  calp2 = 0
  salp1 = 1
  salp2 = 1
  s12x = a * lam12
  sig12 = lam12 / f1
  omg12 = sig12
  m12x = b * sin(sig12)
  if (scalp) then
    MM12 = cos(sig12)
    MM21 = MM12
  end if
  a12x = lon12 / f1
else if (.not. meridian) then
! Now point1 and point2 belong within a hemisphere bounded by a
! meridian and geodesic is neither meridional or equatorial.

! Figure a starting point for Newton's method
  sig12 = InverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, &
      slam12, clam12, f, A3x, salp1, calp1, salp2, calp2, dnm, Ca)

  if (sig12 >= 0) then
! Short lines (InverseStart sets salp2, calp2, dnm)
    s12x = sig12 * b * dnm
    m12x = dnm**2 * b * sin(sig12 / dnm)
    if (scalp) then
      MM12 = cos(sig12 / dnm)
      MM21 = MM12
    end if
    a12x = sig12 / degree
    omg12 = lam12 / (f1 * dnm)
  else

! Newton's method.  This is a straightforward solution of f(alp1) =
! lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
! root in the interval (0, pi) and its derivative is positive at the
! root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
! alp1.  During the course of the iteration, a range (alp1a, alp1b) is
! maintained which brackets the root and with each evaluation of
! f(alp) the range is shrunk, if possible.  Newton's method is
! restarted whenever the derivative of f is negative (because the new
! value of alp1 is then further from the solution) or if the new
! estimate of alp1 lies outside (0,pi); in this case, the new starting
! guess is taken to be (alp1a + alp1b) / 2.

! Bracketing range
    salp1a = tiny2
    calp1a = 1
    salp1b = tiny2
    calp1b = -1
    tripn = .false.
    tripb = .false.
    do numit = 0, maxit2-1
      ! the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
      ! WGS84 and random input: mean = 2.85, sd = 0.60
      v = Lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, &
          salp1, calp1, slam12, clam12, f, A3x, C3x, salp2, calp2, &
          sig12, ssig1, csig1, ssig2, csig2, &
          eps, domg12, numit < maxit1, dv, Ca)
      ! Reversed test to allow escape with NaNs
      if (tripn) then
        dummy = 8
      else
        dummy = 1
      end if
      if (tripb .or. .not. (abs(v) >= dummy * tol0)) exit
      ! Update bracketing values
      if (v > 0 .and. (numit > maxit1 .or. &
          calp1/salp1 > calp1b/salp1b)) then
        salp1b = salp1
        calp1b = calp1
      else if (v < 0 .and. (numit > maxit1 .or. &
            calp1/salp1 < calp1a/salp1a)) then
        salp1a = salp1
        calp1a = calp1
      end if
      if (numit < maxit1 .and. dv > 0) then
        dalp1 = -v/dv
        sdalp1 = sin(dalp1)
        cdalp1 = cos(dalp1)
        nsalp1 = salp1 * cdalp1 + calp1 * sdalp1
        if (nsalp1 > 0 .and. abs(dalp1) < pi) then
          calp1 = calp1 * cdalp1 - salp1 * sdalp1
          salp1 = nsalp1
          call norm(salp1, calp1)
          ! In some regimes we don't get quadratic convergence because
          ! slope -> 0.  So use convergence conditions based on dbleps
          ! instead of sqrt(dbleps).
          tripn = abs(v) <= 16 * tol0
          cycle
        end if
      end if
      ! Either dv was not positive or updated value was outside legal
      ! range.  Use the midpoint of the bracket as the next estimate.
      ! This mechanism is not needed for the WGS84 ellipsoid, but it does
      ! catch problems with more eccentric ellipsoids.  Its efficacy is
      ! such for the WGS84 test set with the starting guess set to alp1 =
      ! 90deg:
      ! the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
      ! WGS84 and random input: mean = 4.74, sd = 0.99
      salp1 = (salp1a + salp1b)/2
      calp1 = (calp1a + calp1b)/2
      call norm(salp1, calp1)
      tripn = .false.
      tripb = abs(salp1a - salp1) + (calp1a - calp1) < tolb &
          .or. abs(salp1 - salp1b) + (calp1 - calp1b) < tolb
    end do
    call Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, &
        cbet1, cbet2, lmask, &
        s12x, m12x, dummy, MM12, MM21, ep2, Ca)
    m12x = m12x * b
    s12x = s12x * b
    a12x = sig12 / degree
    if (areap) then
      sdomg12 = sin(domg12)
      cdomg12 = cos(domg12)
      somg12 = slam12 * cdomg12 - clam12 * sdomg12
      comg12 = clam12 * cdomg12 + slam12 * sdomg12
    end if
  end if
end if

! Convert -0 to 0
s12 = 0 + s12x
if (redlp) m12 = 0 + m12x

if (areap) then
! From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
  salp0 = salp1 * cbet1
  calp0 = hypot(calp1, salp1 * sbet1)
  if (calp0 /= 0 .and. salp0 /= 0) then
! From Lambda12: tan(bet) = tan(sig) * cos(alp)
    ssig1 = sbet1
    csig1 = calp1 * cbet1
    ssig2 = sbet2
    csig2 = calp2 * cbet2
    k2 = calp0**2 * ep2
    eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2)
! Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
    A4 = a**2 * calp0 * salp0 * e2
    call norm(ssig1, csig1)
    call norm(ssig2, csig2)
    call C4f(eps, C4x, Ca)
    B41 = SinCosSeries(.false., ssig1, csig1, Ca, nC4)
    B42 = SinCosSeries(.false., ssig2, csig2, Ca, nC4)
    SS12 = A4 * (B42 - B41)
  else
! Avoid problems with indeterminate sig1, sig2 on equator
    SS12 = 0
  end if

  if (.not. meridian .and. somg12 == 2) then
    somg12 = sin(omg12)
    comg12 = cos(omg12)
  end if

  if (.not. meridian .and. comg12 >= 0.7071d0 &
      .and. sbet2 - sbet1 < 1.75d0) then
! Use tan(Gamma/2) = tan(omg12/2)
! * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
! with tan(x/2) = sin(x)/(1+cos(x))
    domg12 = 1 + comg12
    dbet1 = 1 + cbet1
    dbet2 = 1 + cbet2
    alp12 = 2 * atan2(somg12 * (sbet1 * dbet2 + sbet2 * dbet1), &
        domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) )
  else
! alp12 = alp2 - alp1, used in atan2 so no need to normalize
    salp12 = salp2 * calp1 - calp2 * salp1
    calp12 = calp2 * calp1 + salp2 * salp1
! The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
! salp12 = -0 and alp12 = -180.  However this depends on the sign
! being attached to 0 correctly.  The following ensures the correct
! behavior.
    if (salp12 == 0 .and. calp12 < 0) then
      salp12 = tiny2 * calp1
      calp12 = -1
    end if
    alp12 = atan2(salp12, calp12)
  end if
  SS12 = SS12 + c2 * alp12
  SS12 = SS12 * swapp * lonsign * latsign
! Convert -0 to 0
  SS12 = 0 + SS12
end if

! Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
if (swapp < 0) then
  call swap(salp1, salp2)
  call swap(calp1, calp2)
  if (scalp) call swap(MM12, MM21)
end if

salp1 = salp1 * swapp * lonsign
calp1 = calp1 * swapp * latsign
salp2 = salp2 * swapp * lonsign
calp2 = calp2 * swapp * latsign

azi1 = atan2d(salp1, calp1)
azi2 = atan2d(salp2, calp2)

if (arcp) a12 = a12x

end subroutine inverse

!*****************************************************************************************
!>
!  Determine the area of a geodesic polygon
!
!  Arbitrarily complex polygons are allowed.  In the case of
!  self-intersecting polygons the area is accumulated "algebraically",
!  e.g., the areas of the 2 loops in a figure-8 polygon will partially
!  cancel.  There's no need to "close" the polygon by repeating the
!  first vertex.  The area returned is signed with counter-clockwise
!  traversal being treated as positive.

subroutine area(a, f, lats, lons, n, AA, PP)

  integer, intent(in) :: n !! number of points
  real(wp), intent(in) :: a !! the equatorial radius (meters).
  real(wp), intent(in) :: f !! the flattening of the ellipsoid.  Setting `f = 0` gives
                            !! a sphere.  Negative `f` gives a prolate ellipsoid.
  real(wp), intent(in) :: lats(n) !! an array of the latitudes of the vertices (degrees).
                                  !! lats should be in the range [-90 deg, 90 deg].
  real(wp), intent(in) :: lons(n) !! an array of the longitudes of the vertices (degrees).
  real(wp), intent(out) :: AA !! the (signed) area of the polygon (\(m^2\)).
  real(wp), intent(out) :: PP !! the perimeter of the polygon.

integer i, outmask, cross
real(wp) s12, azi1, azi2, dummy, SS12, b, e2, c2, area0, &
                 Aacc(2), Pacc(2)

outmask = 8
Aacc = 0.0_wp
Pacc = 0.0_wp
cross = 0
do i = 0, n-1
  call inverse(a, f, lats(i+1), lons(i+1), &
      lats(mod(i+1,n)+1), lons(mod(i+1,n)+1), &
      s12, azi1, azi2, outmask, dummy, dummy, dummy, dummy, SS12)
  call accadd(Pacc, s12)
  call accadd(Aacc, -SS12)
  cross = cross + transit(lons(i+1), lons(mod(i+1,n)+1))
end do
PP = Pacc(1)
b = a * (1 - f)
e2 = f * (2 - f)
if (e2 == 0) then
  c2 = a**2
else if (e2 > 0) then
  c2 = (a**2 + b**2 * atanh(sqrt(e2)) / sqrt(e2)) / 2
else
  c2 = (a**2 + b**2 * atan(sqrt(abs(e2))) / sqrt(abs(e2))) / 2
end if
area0 = 4 * pi * c2
call accrem(Aacc, area0)
if (mod(abs(cross), 2) == 1) then
  if (Aacc(1) < 0) then
    call accadd(Aacc, +area0/2)
  else
    call accadd(Aacc, -area0/2)
  end if
end if
if (Aacc(1) > area0/2) then
  call accadd(Aacc, -area0)
else if (Aacc(1) <= -area0/2) then
  call accadd(Aacc, +area0)
end if
AA = Aacc(1)

end subroutine area

!*****************************************************************************************
!>
!

subroutine Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, &
    cbet1, cbet2, outmask, &
    s12b, m12b, m0, MM12, MM21, ep2, Ca)

real(wp),intent(in) :: eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, &
                       cbet1, cbet2, ep2
integer,intent(in) :: outmask
real(wp),intent(out) :: s12b, m12b, m0, MM12, MM21
real(wp) :: Ca(*) !! temporary storage

integer ord, nC1, nC2
parameter (ord = 6, nC1 = ord, nC2 = ord)

real(wp) m0x, J12, A1, A2, B1, B2, csig12, t, Cb(nC2)
logical distp, redlp, scalp
integer l

! Return m12b = (reduced length)/b; also calculate s12b = distance/b,
! and m0 = coefficient of secular term in expression for reduced length.

distp = (mod(outmask/16, 2) == 1)
redlp = (mod(outmask/2, 2) == 1)
scalp = (mod(outmask/4, 2) == 1)

! Suppress compiler warnings
m0x = 0
J12 = 0
A1 = 0
A2 = 0
if (distp .or. redlp .or. scalp) then
  A1 = A1m1f(eps)
  call C1f(eps, Ca)
  if (redlp .or. scalp) then
    A2 = A2m1f(eps)
    call C2f(eps, Cb)
    m0x = A1 - A2
    A2 = 1 + A2
  end if
  A1 = 1 + A1
end if
if (distp) then
  B1 = SinCosSeries(.true., ssig2, csig2, Ca, nC1) - &
      SinCosSeries(.true., ssig1, csig1, Ca, nC1)
! Missing a factor of b
  s12b = A1 * (sig12 + B1)
  if (redlp .or. scalp) then
    B2 = SinCosSeries(.true., ssig2, csig2, Cb, nC2) - &
        SinCosSeries(.true., ssig1, csig1, Cb, nC2)
    J12 = m0x * sig12 + (A1 * B1 - A2 * B2)
  end if
else if (redlp .or. scalp) then
! Assume here that nC1 >= nC2
  do l = 1, nC2
    Cb(l) = A1 * Ca(l) - A2 * Cb(l)
  end do
  J12 = m0x * sig12 + (SinCosSeries(.true., ssig2, csig2, Cb, nC2) - &
      SinCosSeries(.true., ssig1, csig1, Cb, nC2))
end if
if (redlp) then
  m0 = m0x
! Missing a factor of b.
! Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
! accurate cancellation in the case of coincident points.
  m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - &
      csig1 * csig2 * J12
end if
if (scalp) then
  csig12 = csig1 * csig2 + ssig1 * ssig2
  t = ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
  MM12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1
  MM21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2
end if

end subroutine Lengths

!*****************************************************************************************
!>
!  Solve `k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0` for positive root `k`.
!  This solution is adapted from `Geocentric::Reverse`.

real(wp) function Astroid(x, y)

real(wp),intent(in) :: x, y

real(wp) :: k, p, q, r, S, r2, r3, disc, u, &
            T3, T, ang, v, uv, w

p = x**2
q = y**2
r = (p + q - 1) / 6
if ( .not. (q == 0 .and. r < 0) ) then
! Avoid possible division by zero when r = 0 by multiplying equations
! for s and t by r^3 and r, resp.
! S = r^3 * s
  S = p * q / 4
  r2 = r**2
  r3 = r * r2
! The discriminant of the quadratic equation for T3.  This is zero on
! the evolute curve p^(1/3)+q^(1/3) = 1
  disc = S * (S + 2 * r3)
  u = r
  if (disc >= 0) then
    T3 = S + r3
! Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
! of precision due to cancellation.  The result is unchanged because
! of the way the T is used in definition of u.
! T3 = (r * t)^3
    if (T3 < 0) then
      disc = -sqrt(disc)
    else
      disc = sqrt(disc)
    end if
    T3 = T3 + disc
! N.B. cbrt always returns the real root.  cbrt(-8) = -2.
! T = r * t
    T = cbrt(T3)
! T can be zero; but then r2 / T -> 0.
    if (T /= 0) u = u + T + r2 / T
  else
! T is complex, but the way u is defined the result is real.
    ang = atan2(sqrt(-disc), -(S + r3))
! There are three possible cube roots.  We choose the root which
! avoids cancellation.  Note that disc < 0 implies that r < 0.
    u = u + 2 * r * cos(ang / 3)
  end if
! guaranteed positive
  v = sqrt(u**2 + q)
! Avoid loss of accuracy when u < 0.
! u+v, guaranteed positive
  if (u < 0) then
    uv = q / (v - u)
  else
    uv = u + v
  end if
! positive?
  w = (uv - q) / (2 * v)
! Rearrange expression for k to avoid loss of accuracy due to
! subtraction.  Division by 0 not possible because uv > 0, w >= 0.
! guaranteed positive
  k = uv / (sqrt(uv + w**2) + w)
else
! q == 0 && r <= 0
! y = 0 with |x| <= 1.  Handle this case directly.
! for y small, positive root is k = abs(y)/sqrt(1-x^2)
  k = 0
end if
Astroid = k

end function Astroid

!*****************************************************************************************
!>
!  Return a starting point for Newton's method in salp1 and calp1
!  (function value is -1).  If Newton's method doesn't need to be used,
!  return also salp2, calp2, and dnm and function value is sig12.

real(wp) function InverseStart(sbet1, cbet1, dn1, &
                         sbet2, cbet2, dn2, lam12, slam12, clam12, f, A3x, &
                         salp1, calp1, salp2, calp2, dnm, &
                         Ca)

real(wp),intent(in) :: sbet1, cbet1, dn1, sbet2, cbet2, dn2, &
                       lam12, slam12, clam12, f
real(wp),intent(inout) :: A3x(*)
real(wp),intent(out) :: salp1, calp1, salp2, calp2, dnm
real(wp) :: Ca(*) !! temporary

logical :: shortline
real(wp) :: f1, e2, ep2, n, etol2, k2, eps, sig12, &
            sbet12, cbet12, sbet12a, omg12, somg12, comg12, ssig12, csig12, &
            x, y, lamscale, betscale, cbet12a, bt12a, m12b, m0, dummy, &
            k, omg12a, sbetm2, lam12x

f1 = 1 - f
e2 = f * (2 - f)
ep2 = e2 / (1 - e2)
n = f / (2 - f)
! The sig12 threshold for "really short".  Using the auxiliary sphere
! solution with dnm computed at (bet1 + bet2) / 2, the relative error in
! the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
! (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a
! given f and sig12, the max error occurs for lines near the pole.  If
! the old rule for computing dnm = (dn1 + dn2)/2 is used, then the error
! increases by a factor of 2.)  Setting this equal to epsilon gives
! sig12 = etol2.  Here 0.1 is a safety factor (error decreased by 100)
! and max(0.001, abs(f)) stops etol2 getting too large in the nearly
! spherical case.
etol2 = 0.1_wp * tol2 / &
    sqrt( max(0.001_wp, abs(f)) * min(1.0_wp, 1 - f/2) / 2 )

! Return value
sig12 = -1
! bet12 = bet2 - bet1 in [0, pi); bt12a = bet2 + bet1 in (-pi, 0]
sbet12 = sbet2 * cbet1 - cbet2 * sbet1
cbet12 = cbet2 * cbet1 + sbet2 * sbet1
sbet12a = sbet2 * cbet1 + cbet2 * sbet1

shortline = cbet12 >= 0 .and. sbet12 < 0.5_wp .and. &
    cbet2 * lam12 < 0.5_wp

if (shortline) then
  sbetm2 = (sbet1 + sbet2)**2
! sin((bet1+bet2)/2)^2
!  =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
  sbetm2 = sbetm2 / (sbetm2 + (cbet1 + cbet2)**2)
  dnm = sqrt(1 + ep2 * sbetm2)
  omg12 = lam12 / (f1 * dnm)
  somg12 = sin(omg12)
  comg12 = cos(omg12)
else
  somg12 = slam12
  comg12 = clam12
end if

salp1 = cbet2 * somg12
if (comg12 >= 0) then
  calp1 = sbet12 + cbet2 * sbet1 * somg12**2 / (1 + comg12)
else
  calp1 = sbet12a - cbet2 * sbet1 * somg12**2 / (1 - comg12)
end if

ssig12 = hypot(salp1, calp1)
csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12

if (shortline .and. ssig12 < etol2) then
! really short lines
  salp2 = cbet1 * somg12
  if (comg12 >= 0) then
    calp2 = somg12**2 / (1 + comg12)
  else
    calp2 = 1 - comg12
  end if
  calp2 = sbet12 - cbet1 * sbet2 * calp2
  call norm(salp2, calp2)
! Set return value
  sig12 = atan2(ssig12, csig12)
else if (abs(n) > 0.1_wp .or. csig12 >= 0 .or. &
      ssig12 >= 6 * abs(n) * pi * cbet1**2) then
! Nothing to do, zeroth order spherical approximation is OK
else
! lam12 - pi
  lam12x = atan2(-slam12, -clam12)
! Scale lam12 and bet2 to x, y coordinate system where antipodal point
! is at origin and singular point is at y = 0, x = -1.
  if (f >= 0) then
! x = dlong, y = dlat
    k2 = sbet1**2 * ep2
    eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2)
    lamscale = f * cbet1 * A3f(eps, A3x) * pi
    betscale = lamscale * cbet1
    x = lam12x / lamscale
    y = sbet12a / betscale
  else
! f < 0: x = dlat, y = dlong
    cbet12a = cbet2 * cbet1 - sbet2 * sbet1
    bt12a = atan2(sbet12a, cbet12a)
! In the case of lon12 = 180, this repeats a calculation made in
! Inverse.
    call Lengths(n, pi + bt12a, &
        sbet1, -cbet1, dn1, sbet2, cbet2, dn2, cbet1, cbet2, 2, &
        dummy, m12b, m0, dummy, dummy, ep2, Ca)
    x = -1 + m12b / (cbet1 * cbet2 * m0 * pi)
    if (x < -0.01_wp) then
      betscale = sbet12a / x
    else
      betscale = -f * cbet1**2 * pi
    end if
    lamscale = betscale / cbet1
    y = lam12x / lamscale
  end if

  if (y > -tol1 .and. x > -1 - xthresh) then
! strip near cut
    if (f >= 0) then
      salp1 = min(1.0_wp, -x)
      calp1 = - sqrt(1 - salp1**2)
    else
      if (x > -tol1) then
        calp1 = 0
      else
        calp1 = 1
      end if
      calp1 = max(calp1, x)
      salp1 = sqrt(1 - calp1**2)
    end if
  else
! Estimate alp1, by solving the astroid problem.
!
! Could estimate alpha1 = theta + pi/2, directly, i.e.,
!   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
!   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
!
! However, it's better to estimate omg12 from astroid and use
! spherical formula to compute alp1.  This reduces the mean number of
! Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
! (min 0 max 5).  The changes in the number of iterations are as
! follows:
!
! change percent
!    1       5
!    0      78
!   -1      16
!   -2       0.6
!   -3       0.04
!   -4       0.002
!
! The histogram of iterations is (m = number of iterations estimating
! alp1 directly, n = number of iterations estimating via omg12, total
! number of trials = 148605):
!
!  iter    m      n
!    0   148    186
!    1 13046  13845
!    2 93315 102225
!    3 36189  32341
!    4  5396      7
!    5   455      1
!    6    56      0
!
! Because omg12 is near pi, estimate work with omg12a = pi - omg12
    k = Astroid(x, y)
    if (f >= 0) then
      omg12a = -x * k/(1 + k)
    else
      omg12a = -y * (1 + k)/k
    end if
    omg12a = lamscale * omg12a
    somg12 = sin(omg12a)
    comg12 = -cos(omg12a)
! Update spherical estimate of alp1 using omg12 instead of lam12
    salp1 = cbet2 * somg12
    calp1 = sbet12a - cbet2 * sbet1 * somg12**2 / (1 - comg12)
  end if
end if
! Sanity check on starting guess.  Backwards check allows NaN through.
if (.not. (salp1 <= 0)) then
  call norm(salp1, calp1)
else
  salp1 = 1
  calp1 = 0
end if
InverseStart = sig12

end function InverseStart

!*****************************************************************************************
!>
!

real(wp) function Lambda12(sbet1, cbet1, dn1, &
                         sbet2, cbet2, dn2, salp1, calp1, slm120, clm120, f, A3x, C3x, &
                         salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, &
                         domg12, diffp, dlam12, Ca)

real(wp),intent(in) :: sbet1, cbet1, dn1, sbet2, cbet2, dn2, &
                       salp1, slm120, clm120, f, C3x(*)
real(wp),intent(inout) :: A3x(*)
real(wp),intent(inout) :: calp1
logical,intent(in) :: diffp
real(wp),intent(out) :: salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, &
                        eps, domg12
real(wp),intent(out) :: dlam12
real(wp) :: Ca(*)

integer,parameter :: ord = 6
integer,parameter :: nC3 = ord

real(wp) :: f1, e2, ep2, salp0, calp0, &
            somg1, comg1, somg2, comg2, somg12, comg12, &
            lam12, eta, B312, k2, dummy

f1 = 1 - f
e2 = f * (2 - f)
ep2 = e2 / (1 - e2)
! Break degeneracy of equatorial line.  This case has already been
! handled.
if (sbet1 == 0 .and. calp1 == 0) calp1 = -tiny2

! sin(alp1) * cos(bet1) = sin(alp0)
salp0 = salp1 * cbet1
! calp0 > 0
calp0 = hypot(calp1, salp1 * sbet1)

! tan(bet1) = tan(sig1) * cos(alp1)
! tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
ssig1 = sbet1
somg1 = salp0 * sbet1
csig1 = calp1 * cbet1
comg1 = csig1
call norm(ssig1, csig1)
! norm(somg1, comg1); -- don't need to normalize!

! Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
! about this case, since this can yield singularities in the Newton
! iteration.
! sin(alp2) * cos(bet2) = sin(alp0)
if (cbet2 /= cbet1) then
  salp2 = salp0 / cbet2
else
  salp2 = salp1
end if
! calp2 = sqrt(1 - sq(salp2))
!       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
! and subst for calp0 and rearrange to give (choose positive sqrt
! to give alp2 in [0, pi/2]).
if (cbet2 /= cbet1 .or. abs(sbet2) /= -sbet1) then
  if (cbet1 < -sbet1) then
    calp2 = (cbet2 - cbet1) * (cbet1 + cbet2)
  else
    calp2 = (sbet1 - sbet2) * (sbet1 + sbet2)
  end if
  calp2 = sqrt((calp1 * cbet1)**2 + calp2) / cbet2
else
  calp2 = abs(calp1)
end if
! tan(bet2) = tan(sig2) * cos(alp2)
! tan(omg2) = sin(alp0) * tan(sig2).
ssig2 = sbet2
somg2 = salp0 * sbet2
csig2 = calp2 * cbet2
comg2 = csig2
call norm(ssig2, csig2)
! norm(somg2, comg2); -- don't need to normalize!

! sig12 = sig2 - sig1, limit to [0, pi]
sig12 = atan2(max(0d0, csig1 * ssig2 - ssig1 * csig2) + 0d0, &
                       csig1 * csig2 + ssig1 * ssig2)

! omg12 = omg2 - omg1, limit to [0, pi]
somg12 = max(0d0, comg1 * somg2 - somg1 * comg2) + 0d0
comg12 =          comg1 * comg2 + somg1 * somg2
! eta = omg12 - lam120
eta = atan2(somg12 * clm120 - comg12 * slm120, &
    comg12 * clm120 + somg12 * slm120)
k2 = calp0**2 * ep2
eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2)
call C3f(eps, C3x, Ca)
B312 = (SinCosSeries(.true., ssig2, csig2, Ca, nC3-1) - &
    SinCosSeries(.true., ssig1, csig1, Ca, nC3-1))
domg12 = -f * A3f(eps, A3x) * salp0 * (sig12 + B312)
lam12 = eta + domg12

if (diffp) then
  if (calp2 == 0) then
    dlam12 = - 2 * f1 * dn1 / sbet1
  else
    call Lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, &
        cbet1, cbet2, 2, &
        dummy, dlam12, dummy, dummy, dummy, ep2, Ca)
    dlam12 = dlam12 * f1 / (calp2 * cbet2)
  end if
end if
Lambda12 = lam12

end function Lambda12

!*****************************************************************************************
!>
!  Evaluate A3

    real(wp) function A3f(eps, A3x)

    integer,parameter :: ord = 6
    integer,parameter :: nA3 = ord
    integer,parameter :: nA3x = nA3

    real(wp),intent(in) :: eps
    real(wp),intent(out) :: A3x(0: nA3x-1)

    A3f = polyval(nA3 - 1, A3x, eps)

    end function A3f
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate C3 coeffs.
!  Elements c[1] thru c[nC3-1] are set

    subroutine C3f(eps, C3x, c)

    integer,parameter :: ord = 6
    integer,parameter :: nC3 = ord
    integer,parameter :: nC3x = (nC3 * (nC3 - 1)) / 2

    real(wp),intent(in) :: eps, C3x(0:nC3x-1)
    real(wp),intent(out) :: c(nC3-1)

    integer :: o, m, l
    real(wp) :: mult

    mult = 1
    o = 0
    do l = 1, nC3 - 1
        m = nC3 - l - 1
        mult = mult * eps
        c(l) = mult * polyval(m, C3x(o), eps)
        o = o + m + 1
    end do

    end subroutine C3f
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate C4
!  Elements c[0] thru c[nC4-1] are set

    subroutine C4f(eps, C4x, c)

    integer,parameter :: ord = 6
    integer,parameter :: nC4 = ord
    integer,parameter :: nC4x = (nC4 * (nC4 + 1)) / 2

    real(wp),intent(in) :: eps, C4x(0:nC4x-1)
    real(wp),intent(out) :: c(0:nC4-1)

    integer :: o, m, l
    real(wp) :: mult

    mult = 1
    o = 0
    do l = 0, nC4 - 1
    m = nC4 - l - 1
    c(l) = mult * polyval(m, C4x(o), eps)
    o = o + m + 1
    mult = mult * eps
    end do

    end subroutine C4f
!*****************************************************************************************

!*****************************************************************************************
!>
!  The scale factor A1-1 = mean value of (d/dsigma)I1 - 1

    real(wp) function A1m1f(eps)

    integer, parameter :: ord = 6
    integer, parameter :: nA1 = ord

    real(wp),intent(in) :: eps
    real(wp),dimension(nA1/2 + 2),parameter :: coeff = [1, 4, 64, 0, 256]

    real(wp) :: t
    integer :: o, m

    o = 1
    m = nA1/2
    t = polyval(m, coeff(o), eps**2) / coeff(o + m + 1)
    A1m1f = (t + eps) / (1 - eps)

    end function A1m1f
!*****************************************************************************************

!*****************************************************************************************
!>
!  The coefficients C1[l] in the Fourier expansion of B1

    subroutine C1f(eps, c)

    integer,parameter :: ord = 6
    integer,parameter :: nC1 = ord
    integer,parameter :: ncoeff = (nC1**2 + 7*nC1 - 2*(nC1/2))/4

    real(wp),intent(in) :: eps
    real(wp),intent(out) :: c(nC1)

    real(wp) :: eps2, d
    integer :: o, m, l

    real(wp),dimension(ncoeff),parameter :: coeff = [ -1, 6, -16, 32, &
                                                      -9, 64, -128, 2048, &
                                                      9, -16, 768, &
                                                      3, -5, 512, &
                                                      -7, 1280, &
                                                      -7, 2048 ]

    eps2 = eps**2
    d = eps
    o = 1
    do l = 1, nC1
        m = (nC1 - l) / 2
        c(l) = d * polyval(m, coeff(o), eps2) / coeff(o + m + 1)
        o = o + m + 2
        d = d * eps
    end do

    end subroutine C1f
!*****************************************************************************************

!*****************************************************************************************
!>
!  The coefficients C1p[l] in the Fourier expansion of B1p

    subroutine C1pf(eps, c)

    integer,parameter :: ord = 6
    integer,parameter :: nC1p = ord
    integer,parameter :: ncoeff = (nC1p**2 + 7*nC1p - 2*(nC1p/2))/4

    real(wp),intent(in) :: eps
    real(wp),intent(out) :: c(nC1p)

    real(wp) :: eps2, d
    integer :: o, m, l

    real(wp),dimension(ncoeff),parameter :: coeff = [ 205, -432, 768, 1536, &
                                                      4005, -4736, 3840, 12288, &
                                                      -225, 116, 384, &
                                                      -7173, 2695, 7680, &
                                                      3467, 7680, &
                                                      38081, 61440 ]

    eps2 = eps**2
    d = eps
    o = 1
    do l = 1, nC1p
        m = (nC1p - l) / 2
        c(l) = d * polyval(m, coeff(o), eps2) / coeff(o + m + 1)
        o = o + m + 2
        d = d * eps
    end do

    end subroutine C1pf
!*****************************************************************************************

!*****************************************************************************************
!>
!  The scale factor A2-1 = mean value of (d/dsigma)I2 - 1

    real(wp) function A2m1f(eps)

    real(wp),intent(in) :: eps

    real(wp) :: t
    integer :: o, m

    integer,parameter :: ord = 6
    integer,parameter :: nA2 = ord
    integer,parameter :: ncoeff = nA2/2 + 2

    real(wp),dimension(ncoeff),parameter :: coeff = [-11, -28, -192, 0, 256]

    o = 1
    m = nA2/2
    t = polyval(m, coeff(o), eps**2) / coeff(o + m + 1)
    A2m1f = (t - eps) / (1 + eps)

    end function A2m1f
!*****************************************************************************************

!*****************************************************************************************
!>
!  The coefficients C2[l] in the Fourier expansion of B2

    subroutine C2f(eps, c)

    integer,parameter :: ord = 6
    integer,parameter :: nC2 = ord
    integer,parameter :: ncoeff = (nC2**2 + 7*nC2 - 2*(nC2/2))/4

    real(wp),intent(in) :: eps
    real(wp),intent(out) :: c(nC2)

    real(wp) :: eps2, d
    integer :: o, m, l

    real(wp),dimension(ncoeff),parameter :: coeff = [ 1, 2, 16, 32, &
                                                      35, 64, 384, 2048, &
                                                      15, 80, 768, &
                                                      7, 35, 512, &
                                                      63, 1280, &
                                                      77, 2048 ]

    eps2 = eps**2
    d = eps
    o = 1
    do l = 1, nC2
        m = (nC2 - l) / 2
        c(l) = d * polyval(m, coeff(o), eps2) / coeff(o + m + 1)
        o = o + m + 2
        d = d * eps
    end do

    end subroutine C2f
!*****************************************************************************************

!*****************************************************************************************
!>
!  The scale factor A3 = mean value of (d/dsigma)I3

    subroutine A3coeff(n, A3x)

    integer,parameter :: ord = 6
    integer,parameter :: nA3 = ord
    integer,parameter :: nA3x = nA3
    integer,parameter :: ncoeff = (nA3**2 + 7*nA3 - 2*(nA3/2))/4

    real(wp),intent(in) :: n
    real(wp),intent(out) :: A3x(0:nA3x-1)

    integer :: o, m, k, j
    real(wp),dimension(ncoeff),parameter :: coeff = [ -3, 128, &
                                                      -2, -3, 64, &
                                                      -1, -3, -1, 16, &
                                                      3, -1, -2, 8, &
                                                      1, -1, 2, &
                                                      1, 1 ]

    o = 1
    k = 0
    do j = nA3 - 1, 0, -1
        m = min(nA3 - j - 1, j)
        A3x(k) = polyval(m, coeff(o), n) / coeff(o + m + 1)
        k = k + 1
        o = o + m + 2
    end do

    end subroutine A3coeff
!*****************************************************************************************

!*****************************************************************************************
!>
!  The coefficients C3[l] in the Fourier expansion of B3

    subroutine C3coeff(n, C3x)

    integer,parameter :: ord = 6
    integer,parameter :: nC3 = ord
    integer,parameter :: nC3x = (nC3 * (nC3 - 1)) / 2
    integer,parameter :: ncoeff = ((nC3-1)*(nC3**2 + 7*nC3 - 2*(nC3/2)))/8

    real(wp),intent(in) :: n
    real(wp),intent(out) :: C3x(0:nC3x-1)

    integer :: o, m, l, j, k
    real(wp),dimension(ncoeff),parameter :: coeff = [ 3, 128, &
                                                    2, 5, 128, &
                                                    -1, 3, 3, 64, &
                                                    -1, 0, 1, 8, &
                                                    -1, 1, 4, &
                                                    5, 256, &
                                                    1, 3, 128, &
                                                    -3, -2, 3, 64, &
                                                    1, -3, 2, 32, &
                                                    7, 512, &
                                                    -10, 9, 384, &
                                                    5, -9, 5, 192, &
                                                    7, 512, &
                                                    -14, 7, 512, &
                                                    21, 2560 ]

    o = 1
    k = 0
    do l = 1, nC3 - 1
        do j = nC3 - 1, l, -1
            m = min(nC3 - j - 1, j)
            C3x(k) = polyval(m, coeff(o), n) / coeff(o + m + 1)
            k = k + 1
            o = o + m + 2
        end do
    end do

    end subroutine C3coeff
!*****************************************************************************************

!*****************************************************************************************
!>
!  The coefficients C4[l] in the Fourier expansion of I4

    subroutine C4coeff(n, C4x)

    integer,parameter :: ord = 6
    integer,parameter :: nC4 = ord
    integer,parameter :: nC4x = (nC4 * (nC4 + 1)) / 2
    integer,parameter :: ncoeff = (nC4 * (nC4 + 1) * (nC4 + 5)) / 6

    real(wp),intent(in) :: n
    real(wp),intent(out) :: C4x(0:nC4x-1)

    integer :: o, m, l, j, k
    real(wp),dimension(ncoeff),parameter :: coeff = [ &
        97, 15015,   1088, 156, 45045,   -224, -4784, 1573, 45045, &
        -10656, 14144, -4576, -858, 45045, &
        64, 624, -4576, 6864, -3003, 15015, &
        100, 208, 572, 3432, -12012, 30030, 45045, &
        1, 9009,   -2944, 468, 135135,   5792, 1040, -1287, 135135, &
        5952, -11648, 9152, -2574, 135135, &
        -64, -624, 4576, -6864, 3003, 135135, &
        8, 10725,   1856, -936, 225225,   -8448, 4992, -1144, 225225, &
        -1440, 4160, -4576, 1716, 225225, &
        -136, 63063,   1024, -208, 105105, &
        3584, -3328, 1144, 315315, &
        -128, 135135,   -2560, 832, 405405,   128, 99099 ]

    o = 1
    k = 0
    do l = 0, nC4 - 1
    do j = nC4 - 1, l, -1
        m = nC4 - j - 1
        C4x(k) = polyval(m, coeff(o), n) / coeff(o + m + 1)
        k = k + 1
        o = o + m + 2
    end do
    end do

    end subroutine C4coeff
!*****************************************************************************************

!*****************************************************************************************
!>
!
    real(wp) function sumx(u, v, t)

    real(wp),intent(in) :: u, v
    real(wp),intent(out) :: t

    real(wp) :: up, vpp

    sumx = u + v
    up = sumx - v
    vpp = sumx - up
    up = up - u
    vpp = vpp - v
    if (sumx == 0.0_wp) then
        t = sumx
    else
        t = 0.0_wp - (up + vpp)
    end if

    end function sumx
!*****************************************************************************************

!*****************************************************************************************
!>
!  the remainder function but not worrying how ties are handled
!  y must be positive

    real(wp) function remx(x, y)

    real(wp),intent(in) :: x, y

    remx = mod(x, y)
    if (remx < -y/2.0_wp) then
        remx = remx + y
    else if (remx > +y/2.0_wp) then
        remx = remx - y
    end if

    end function remx
!*****************************************************************************************

!*****************************************************************************************
!>
!
    real(wp) function AngNormalize(x)

    real(wp),intent(in) :: x

    AngNormalize = remx(x, 360.0_wp)
    if (abs(AngNormalize) == 180.0_wp) then
        AngNormalize = sign(180.0_wp, x)
    end if

    end function AngNormalize
!*****************************************************************************************

!*****************************************************************************************
!>
    real(wp) function LatFix(x)

    real(wp),intent(in) :: x

    LatFix = x
    if (.not. (abs(x) > 90.0_wp)) return
    ! concoct a NaN
    LatFix = sqrt(90.0_wp - abs(x))

    end function LatFix
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute y - x.  x and y must both lie in [-180, 180].  The result is
!  equivalent to computing the difference exactly, reducing it to (-180,
!  180] and rounding the result.  Note that this prescription allows -180
!  to be returned (e.g., if x is tiny and negative and y = 180).  The
!  error in the difference is returned in e

    real(wp) function AngDiff(x, y, e)

    real(wp),intent(in) :: x, y
    real(wp),intent(out) :: e

    real(wp) d, t

    d = sumx(remx(-x, 360.0_wp), remx(y, 360.0_wp), t)
    d = sumx(remx(d, 360.0_wp), t, e)
    if (d == 0 .or. abs(d) == 180) then
        if (e == 0) then
            d = sign(d, y - x)
        else
            d = sign(d, -e)
        end if
    end if
    AngDiff = d

    end function AngDiff
!*****************************************************************************************

!*****************************************************************************************
!>
!  The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
!  for reals = 0.7 pm on the earth if x is an angle in degrees.  (This is
!  about 1000 times more resolution than we get with angles around 90
!  degrees.)  We use this to avoid having to deal with near singular
!  cases when x is non-zero but tiny (e.g., 1.0e-200).

    real(wp) function AngRound(x)

    real(wp),intent(in) :: x

    real(wp) :: y, z

    z = 1.0_wp/16.0_wp
    y = abs(x)
    ! The compiler mustn't "simplify" z - (z - y) to y
    if (y < z) y = z - (z - y)
    AngRound = sign(y, x)

    end function AngRound
!*****************************************************************************************

!*****************************************************************************************
!>
!  Swap two real values

    subroutine swap(x, y)

    real(wp),intent(inout) :: x, y

    real(wp) :: z

    z = x
    x = y
    y = z

    end subroutine swap
!*****************************************************************************************

!*****************************************************************************************
!>

    subroutine norm(x, y)

    real(wp),intent(inout) :: x, y

    real(wp) :: r

    r = hypot(x, y)
    x = x/r
    y = y/r

    end subroutine norm
!*****************************************************************************************

!*****************************************************************************************
!>

    real(wp) function log1px(x)

    real(wp),intent(in) :: x

    real(wp) :: y, z

    y = 1 + x
    z = y - 1
    if (z == 0) then
        log1px = x
    else
        log1px = x * log(y) / z
    end if

    end function log1px
!*****************************************************************************************

! real(wp) function atanhx(x)
! ! input
! real(wp) x
!
! ! With Fortran 2008, this becomes: atanhx = atanh(x)
! real(wp) y
! y = abs(x)
! y = log1px(2 * y/(1 - y))/2
! atanhx = sign(y, x)
!
! end function atanhx

!*****************************************************************************************
!>
!  Cube root function

    real(wp) function cbrt(x)

    real(wp),intent(in) :: x

    cbrt = sign(abs(x)**(1.0_wp/3.0_wp), x)

    end function cbrt
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate
!```
!  y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
!             sum(c[i] * cos((2*i-1) * x), i, 1, n)
!```
!  using Clenshaw summation.
!  Approx operation count = (n + 5) mult and (2 * n + 2) add

real(wp) function SinCosSeries(sinp, sinx, cosx, c, n)

logical,intent(in) :: sinp
integer,intent(in) :: n
real(wp),intent(in) :: sinx, cosx, c(n)

real(wp) :: ar, y0, y1
integer :: n2, k

! 2 * cos(2 * x)
ar = 2 * (cosx - sinx) * (cosx + sinx)
! accumulators for sum
if (mod(n, 2) == 1) then
  y0 = c(n)
  n2 = n - 1
else
  y0 = 0
  n2 = n
end if
y1 = 0
! Now n2 is even
do k = n2, 2, -2
  ! Unroll loop x 2, so accumulators return to their original role
  y1 = ar * y0 - y1 + c(k)
  y0 = ar * y1 - y0 + c(k-1)
end do
if (sinp) then
! sin(2 * x) * y0
  SinCosSeries = 2 * sinx * cosx * y0
else
! cos(x) * (y0 - y1)
  SinCosSeries = cosx * (y0 - y1)
end if

end function SinCosSeries
!*****************************************************************************************

!*****************************************************************************************
!>

    integer function transit(lon1, lon2)

    real(wp),intent(in) :: lon1, lon2

    real(wp) :: lon1x, lon2x, lon12, e

    lon12 = AngDiff(lon1, lon2, e)
    lon1x = AngNormalize(lon1)
    lon2x = AngNormalize(lon2)
    if (lon12 > 0 .and. ((lon1x < 0 .and. lon2x >= 0) .or. &
                         (lon1x > 0 .and. lon2x == 0))) then
        transit = 1
    else if (lon12 < 0 .and. lon1x >= 0 .and. lon2x < 0) then
        transit = -1
    else
        transit = 0
    end if

    end function transit
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add y to an accumulator.

    subroutine accadd(s, y)

    real(wp),intent(in) :: y
    real(wp),intent(inout) :: s(2)

    real(wp) :: z, u

    z = sumx(y, s(2), u)
    s(1) = sumx(z, s(1), s(2))
    if (s(1) == 0) then
        s(1) = u
    else
        s(2) = s(2) + u
    end if

    end subroutine accadd
!*****************************************************************************************

!*****************************************************************************************
!>
!  Reduce s to [-y/2, y/2].

    subroutine accrem(s, y)

    real(wp),intent(in) :: y
    real(wp),intent(inout) :: s(2)

    s(1) = remx(s(1), y)
    call accadd(s, 0.0_wp)

    end subroutine accrem
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute `sin(x)` and `cos(x)` with `x` in degrees

    pure subroutine sincosd(x, sinx, cosx)

    real(wp),intent(in) :: x
    real(wp),intent(out) :: sinx, cosx

    real(wp) :: r, s, c
    integer :: q

    r = mod(x, 360.0_wp)
    q = nint(r / 90.0_wp)
    r = (r - 90.0_wp * q) * degree
    s = sin(r)
    c = cos(r)
    q = mod(q + 4, 4)
    select case (q)
    case(0)
        sinx =  s
        cosx =  c
    case(1)
        sinx =  c
        cosx = -s
    case(2)
        sinx = -s
        cosx = -c
    case(3)
        sinx = -c
        cosx =  s
    end select

    if (sinx == 0.0_wp) then
        sinx = sign(sinx, x)
    end if
    cosx = 0.0_wp + cosx

    end subroutine sincosd
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute sin(x+t) and cos(x+t) with x in degrees

    subroutine sincosde(x, t, sinx, cosx)

    real(wp),intent(in) :: x, t
    real(wp),intent(inout) :: sinx, cosx

    real(wp) :: r, s, c
    integer :: q

    q = nint(x / 90.0_wp)
    r = x - 90.0_wp * q
    r = AngRound(r + t) * degree
    s = sin(r)
    c = cos(r)
    q = mod(q + 4, 4)
    select case (q)
    case (0)
        sinx =  s
        cosx =  c
    case (1)
        sinx =  c
        cosx = -s
    case (2)
        sinx = -s
        cosx = -c
    case (3)
        sinx = -c
        cosx =  s
    end select

    if (sinx == 0.0_wp) then
    sinx = sign(sinx, x)
    end if
    cosx = 0.0_wp + cosx

    end subroutine sincosde
!*****************************************************************************************

!*****************************************************************************************
!>

    real(wp) function atan2d(y, x)

    real(wp),intent(in) :: x, y

    real(wp) :: xx, yy
    integer :: q

    if (abs(y) > abs(x)) then
        xx = y
        yy = x
        q = 2
    else
        xx = x
        yy = y
        q = 0
    end if
    if (xx < 0) then
        xx = -xx
        q = q + 1
    end if
    atan2d = atan2(yy, xx) / degree
    if (q == 1) then
        atan2d = sign(180.0_wp, y) - atan2d
    else if (q == 2) then
        atan2d =       90       - atan2d
    else if (q == 3) then
        atan2d =      -90       + atan2d
    end if

    end function atan2d
!*****************************************************************************************

!*****************************************************************************************
!>

    real(wp) function polyval(N, p, x)

    integer,intent(in) :: N
    real(wp),intent(in) :: p(0:N), x

    integer i
    if (N < 0) then
        polyval = 0
    else
        polyval = p(0)
    end if
    do i = 1, N
        polyval = polyval * x + p(i)
    end do

    end function polyval
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Heikkinen routine for cartesian to geodetic transformation
!
!### References
!  1. M. Heikkinen, "Geschlossene formeln zur berechnung raumlicher
!     geodatischer koordinaten aus rechtwinkligen Koordinaten".
!     Z. Ermess., 107 (1982), 207-211 (in German).
!  2. E. D. Kaplan, "Understanding GPS: Principles and Applications",
!     Artech House, 1996.

    pure subroutine heikkinen(rvec, a, b, h, lon, lat)

    implicit none

    real(wp),dimension(3),intent(in) :: rvec  !! position vector [km]
    real(wp),intent(in)  :: a                 !! geoid semimajor axis [km]
    real(wp),intent(in)  :: b                 !! geoid semiminor axis [km]
    real(wp),intent(out) :: h                 !! geodetic altitude [km]
    real(wp),intent(out) :: lon               !! longitude [rad]
    real(wp),intent(out) :: lat               !! geodetic latitude [rad]

    real(wp) :: f,e_2,ep,r,e2,ff,g,c,s,pp,q,r0,u,v,z0,x,y,z,z2,r2,tmp,a2,b2

    x   = rvec(1)
    y   = rvec(2)
    z   = rvec(3)
    a2  = a*a
    b2  = b*b
    f   = (a-b)/a
    e_2 = (2.0_wp*f-f*f)
    ep  = sqrt(a2/b2 - 1.0_wp)
    z2  = z*z
    r   = sqrt(x**2 + y**2)
    r2  = r*r
    e2  = a2 - b2
    ff  = 54.0_wp * b2 * z2
    g   = r2 + (1.0_wp - e_2)*z2 - e_2*e2
    c   = e_2**2 * ff * r2 / g**3
    s   = (1.0_wp + c + sqrt(c**2 + 2.0_wp*c))**(1.0_wp/3.0_wp)
    pp  = ff / ( 3.0_wp*(s + 1.0_wp/s + 1.0_wp)**2 * g**2 )
    q   = sqrt( 1.0_wp + 2.0_wp*e_2**2 * pp )
    r0  = -pp*e_2*r/(1.0_wp+q) + &
            sqrt( max(0.0_wp, 1.0_wp/2.0_wp * a2 * (1.0_wp + 1.0_wp/q) - &
                ( pp*(1.0_wp-e_2)*z2 )/(q*(1.0_wp+q)) - &
                1.0_wp/2.0_wp * pp * r2) )
    u   = sqrt( (r - e_2*r0)**2 + z2 )
    v   = sqrt( (r - e_2*r0)**2 + (1.0_wp - e_2)*z2 )
    z0  = b**2 * z / (a*v)

    h   = u*(1.0_wp - b2/(a*v) )
    lat = atan2( (z + ep**2*z0), r )
    lon = atan2( y, x )

    end subroutine heikkinen
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Olson routine for cartesian to geodetic transformation.
!
!# References
!  1. Olson, D. K., Converting Earth-Centered, Earth-Fixed Coordinates to
!     Geodetic Coordinates, IEEE Transactions on Aerospace and Electronic
!     Systems, 32 (1996) 473-476.

    pure subroutine olson(rvec, a, b, h, long, lat)

    implicit none

    real(wp),dimension(3),intent(in) :: rvec !!position vector [km]
    real(wp),intent(in)  :: a                !!geoid semimajor axis [km]
    real(wp),intent(in)  :: b                !!geoid semiminor axis [km]
    real(wp),intent(out) :: h                !!geodetic altitude [km]
    real(wp),intent(out) :: long             !!longitude [rad]
    real(wp),intent(out) :: lat              !!geodetic latitude [rad]

    real(wp) :: f,x,y,z,e2,a1,a2,a3,a4,a5,a6,w,zp,&
                w2,r2,r,s2,c2,u,v,s,ss,c,g,rg,rf,m,p,z2

    x  = rvec(1)
    y  = rvec(2)
    z  = rvec(3)
    f  = (a-b)/a
    e2 = f * (2.0_wp - f)
    a1 = a * e2
    a2 = a1 * a1
    a3 = a1 * e2 / 2.0_wp
    a4 = 2.5_wp * a2
    a5 = a1 + a3
    a6 = 1.0_wp - e2
    zp = abs(z)
    w2 = x*x + y*y
    w  = sqrt(w2)
    z2 = z * z
    r2 = z2 + w2
    r  = sqrt(r2)

    if (r < 100.0_wp) then

        lat = 0.0_wp
        long = 0.0_wp
        h = -1.0e7_wp

    else

        s2 = z2 / r2
        c2 = w2 / r2
        u  = a2 / r
        v  = a3 - a4 / r

        if (c2 > 0.3_wp) then
            s = (zp / r) * (1.0_wp + c2 * (a1 + u + s2 * v) / r)
            lat = asin(s)
            ss = s * s
            c = sqrt(1.0_wp - ss)
        else
            c = (w / r) * (1.0_wp - s2 * (a5 - u - c2 * v) / r)
            lat = acos(c)
            ss = 1.0_wp - c * c
            s = sqrt(ss)
        end if

        g   = 1.0_wp - e2 * ss
        rg  = a / sqrt(g)
        rf  = a6 * rg
        u   = w - rg * c
        v   = zp - rf * s
        f   = c * u + s * v
        m   = c * v - s * u
        p   = m / (rf / g + f)
        lat = lat + p
        if (z < 0.0_wp) lat = -lat
        h = f + m * p / 2.0_wp
        long = atan2( y, x )

    end if

    end subroutine olson
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Solve the "direct" geodetic problem: given the latitude and longitude of one
!  point and the azimuth and distance to a second point, determine the latitude
!  and longitude of that second point.  The solution is obtained using the
!  algorithm by Vincenty.
!
!# References
!  1. T. Vincenty, "[Direct and Inverse Solutions of Geodesics on the
!     Ellipsoid with Application of Nested Equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)",
!     Survey Review XXII. 176, April 1975.
!  2. [PC Software Download - INVERSE and FORWARD](http://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/),
!     National Geodetic Survey. Version 3.0 (November, 2012).

    subroutine direct_vincenty(a,f,glat1,glon1,faz,s,glat2,glon2,baz)

    implicit none

    real(wp),intent(in)  :: a        !! semimajor axis of ellipsoid [m]
    real(wp),intent(in)  :: f        !! flattening of ellipsoid [-]
    real(wp),intent(in)  :: glat1    !! latitude of 1 [rad]
    real(wp),intent(in)  :: glon1    !! longitude of 1 [rad]
    real(wp),intent(in)  :: faz      !! forward azimuth 1->2 [rad]
    real(wp),intent(in)  :: s        !! distance from 1->2 [m]
    real(wp),intent(out) :: glat2    !! latitude of 2 [rad]
    real(wp),intent(out) :: glon2    !! longitude of 2 [rad]
    real(wp),intent(out) :: baz      !! back azimuth 2->1 [rad]

    real(wp) :: r,tu,sf,cf,cu,su,sa,c2a,x,c,d,y,sy,cy,cz,e

    real(wp),parameter :: eps = 0.5e-13_wp

    r  = 1.0_wp - f
    tu = r*sin(glat1)/cos(glat1)
    sf = sin(faz)
    cf = cos(faz)
    if ( cf/=0.0_wp ) then
        baz = atan2(tu,cf)*2.0_wp
    else
        baz = 0.0_wp
    end if
    cu  = 1.0_wp/sqrt(tu*tu+1.0_wp)
    su  = tu*cu
    sa  = cu*sf
    c2a = -sa*sa + 1.0_wp
    x   = sqrt((1.0_wp/r/r-1.0_wp)*c2a+1.0_wp) + 1.0_wp
    x   = (x-2.0_wp)/x
    c   = 1.0_wp - x
    c   = (x*x/4.0_wp+1.0_wp)/c
    d   = (0.375_wp*x*x-1.0_wp)*x
    tu  = s/r/a/c
    y   = tu
    do
        sy = sin(y)
        cy = cos(y)
        cz = cos(baz+y)
        e  = cz*cz*2.0_wp - 1.0_wp
        c  = y
        x  = e*cy
        y  = e + e - 1.0_wp
        y  = (((sy*sy*4.0_wp-3.0_wp)*y*cz*d/6.0_wp+x)*d/4.0_wp-cz)*sy*d + tu
        if ( abs(y-c)<=eps ) exit
    end do
    baz   = cu*cy*cf - su*sy
    c     = r*sqrt(sa*sa+baz*baz)
    d     = su*cy + cu*sy*cf
    glat2 = atan2(d,c)
    c     = cu*cy - su*sy*cf
    x     = atan2(sy*sf,c)
    c     = ((-3.0_wp*c2a+4.0_wp)*f+4.0_wp)*c2a*f/16.0_wp
    d     = ((e*cy*c+cz)*sy*c+y)*sa
    glon2 = glon1 + x - (1.0_wp-c)*d*f
    baz   = atan2(sa,baz) + pi

    end subroutine direct_vincenty
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Geodetic latitude, longitude, and height to Cartesian position vector.
!
!# References
!  1. E. D. Kaplan, "Understanding GPS: Principles and Applications",
!     Artech House, 1996.

    subroutine geodetic_to_cartesian(a,b,glat,lon,h,r)

    implicit none

    real(wp),intent(in) :: a                !! geoid semimajor axis [km]
    real(wp),intent(in) :: b                !! geoid semiminor axis [km]
    real(wp),intent(in) :: glat             !! geodetic latitude [rad]
    real(wp),intent(in) :: lon              !! longitude [rad]
    real(wp),intent(in) :: h                !! geodetic altitude [km]
    real(wp),dimension(3),intent(out) :: r  !! Cartesian position vector [x,y,z]

    real(wp) :: e2,slat,clat,slon,clon,tlat,ome2,d,q,aod

    slat    = sin(glat)
    clat    = cos(glat)
    tlat    = tan(glat)
    slon    = sin(lon)
    clon    = cos(lon)
    e2      = 1.0_wp - (b*b)/(a*a)
    ome2    = 1.0_wp - e2
    d       = sqrt( 1.0_wp + ome2*tlat*tlat )
    q       = sqrt( 1.0_wp - e2*slat*slat   )
    aod     = a/d

    r(1) = aod*clon + h*clon*clat
    r(2) = aod*slon + h*slon*clat
    r(3) = a*ome2*slat/q + h*slat

    end subroutine geodetic_to_cartesian
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 7/13/2014
!
!  Great circle distance on a spherical body, using the Vincenty algorithm.
!
!# References
!  * T. Vincenty, "[Direct and Inverse Solutions of Geodesics on the Ellipsoid
!    with Application of Nested Equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)",
!    Survey Review XXII. 176, April 1975.

    pure function great_circle_distance(r,long1,lat1,long2,lat2) result(d)

    implicit none

    real(wp)            :: d        !! great circle distance from 1 to 2 [km]
    real(wp),intent(in) :: r        !! radius of the body [km]
    real(wp),intent(in) :: long1    !! longitude of first site [rad]
    real(wp),intent(in) :: lat1     !! latitude of the first site [rad]
    real(wp),intent(in) :: long2    !! longitude of the second site [rad]
    real(wp),intent(in) :: lat2     !! latitude of the second site [rad]

    real(wp) :: c1,s1,c2,s2,dlon,clon,slon

    !Compute aux variables:
    c1    = cos(lat1)
    s1    = sin(lat1)
    c2    = cos(lat2)
    s2    = sin(lat2)
    dlon  = long1-long2
    clon  = cos(dlon)
    slon  = sin(dlon)

    d = r*atan2( sqrt((c2*slon)**2+(c1*s2-s1*c2*clon)**2), (s1*s2+c1*c2*clon) )

    end function great_circle_distance
!*****************************************************************************************

!*****************************************************************************************
!>
!  The distance from the center of a celestial body (e.g., the Earth) to a point
!  on the spheroid surface at a specified geodetic latitude.
!
!### Reference
!  * [Geocentric radius](https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius)

    pure function geocentric_radius(a,b,lat) result(r)

    implicit none

    real(wp),intent(in) :: a    !! equatorial radius (km)
    real(wp),intent(in) :: b    !! polar radius of point (km)
    real(wp),intent(in) :: lat  !! geodetic latitude of point (rad)
    real(wp)            :: r    !! distance from center of body to point (km)

    !local variables:
    real(wp) :: num,den,cl2,sl2,a2,b2

    if (a==zero .and. b==zero) then
        r = zero
    else
        cl2 = cos(lat)**2
        sl2 = sin(lat)**2
        a2  = a*a
        b2  = b*b
        num = cl2 * a2**2 + sl2 * b2**2
        den = cl2 * a2    + sl2 * b2
        r   = sqrt(num/den)
    end if

    end function geocentric_radius
!*****************************************************************************************

!*****************************************************************************************
!>
!  INVERSE computes the geodetic azimuth and distance between two points,
!  given their geographic positions.
!
!  Version for long-line and antipodal cases.
!  Latitudes may be 90 degrees exactly.
!
!### Reference
!  * T. Vincenty, "[Direct and Inverse Solutions of Geodesics on the Ellipsoid
!    with Application of Nested Equations](http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)",
!    Survey Review XXII. 176, April 1975.
!  * [inverse.for](http://www.ngs.noaa.gov/PC_PROD/Inv_Fwd/source/inverse.for)
!    Version 3.0 (November, 2012).
!
!### History
!  * Original programmed by thaddeus vincenty, 1975, 1976
!  * Removed back side solution option, debugged, revised -- 2011may01 -- dgm
!    this version of code is interim -- antipodal boundary needs work
!  * Jacob Williams, 1/25/2016 : refactored into modern Fortran.

    subroutine inverse_vincenty(a,rf,b1,l1,b2,l2,faz,baz,s,it,sig,lam,kind)

    implicit none

    real(wp),intent(in)  :: a     !! Equatorial semimajor axis
    real(wp),intent(in)  :: rf    !! reciprocal flattening (1/f)
    real(wp),intent(in)  :: b1    !! latitude of point 1 (rad, positive north)
    real(wp),intent(in)  :: l1    !! longitude of point 1 (rad, positive east)
    real(wp),intent(in)  :: b2    !! latitude of point 2 (rad, positive north)
    real(wp),intent(in)  :: l2    !! longitude of point 2 (rad, positive east)
    real(wp),intent(out) :: faz   !! Forward azimuth (rad, clockwise from north)
    real(wp),intent(out) :: baz   !! Back azimuth (rad, clockwise from north)
    real(wp),intent(out) :: s     !! Ellipsoidal distance
    integer,intent(out)  :: it    !! iteration count
    real(wp),intent(out) :: sig   !! spherical distance on auxiliary sphere
    real(wp),intent(out) :: lam   !! longitude difference on auxiliary sphere
    integer,intent(out)  :: kind  !! solution flag: kind=1, long-line; kind=2, antipodal

    real(wp) :: beta1,beta2,biga,bigb,bige,bigf,boa,c,cosal2,coslam,&
                cossig,costm,costm2,cosu1,cosu2,d,dsig,ep2,l,prev,&
                sinal,sinlam,sinsig,sinu1,sinu2,tem1,tem2,temp,test,z

    real(wp),parameter :: tol = 1.0e-14_wp   !! convergence tolerance
    real(wp),parameter :: eps = 1.0e-15_wp   !! tolerance for zero

    boa = 1.0_wp - 1.0_wp/rf   ! b/a

    beta1 = atan(boa*tan(b1))  ! better reduced latitude
    sinu1 = sin(beta1)
    cosu1 = cos(beta1)
    beta2 = atan(boa*tan(b2))  ! better reduced latitude
    sinu2 = sin(beta2)
    cosu2 = cos(beta2)

    l = l2 - l1  ! longitude difference [-pi,pi]
    if ( l>pi ) l = l - twopi
    if ( l<-pi ) l = l + twopi
    prev = l
    test = l
    it   = 0
    kind = 1
    lam  = l

    longline : do  ! long-line loop (kind=1)

        sinlam = sin(lam)
        coslam = cos(lam)
        temp   = cosu1*sinu2 - sinu1*cosu2*coslam
        sinsig = sqrt((cosu2*sinlam)**2+temp**2)
        cossig = sinu1*sinu2 + cosu1*cosu2*coslam
        sig    = atan2(sinsig,cossig)
        if ( abs(sinsig)<eps ) then
            sinal = cosu1*cosu2*sinlam/sign(eps,sinsig)
        else
            sinal = cosu1*cosu2*sinlam/sinsig
        endif
        cosal2 = -sinal**2 + 1.0_wp
        if ( abs(cosal2)<eps ) then
            costm = -2.0_wp*(sinu1*sinu2/sign(eps,cosal2)) + cossig
        else
            costm = -2.0_wp*(sinu1*sinu2/cosal2) + cossig
        endif
        costm2 = costm*costm
        c = ((-3.0_wp*cosal2+4.0_wp)/rf+4.0_wp)*cosal2/rf/16.0_wp

        antipodal : do  ! antipodal loop (kind=2)

            it = it + 1
            d = (((2.0_wp*costm2-1.0_wp)*cossig*c+costm)*sinsig*c+sig)*(1.0_wp-c)/rf
            if ( kind==1 ) then
                lam = l + d*sinal
                if ( abs(lam-test)>=tol ) then
                    if ( abs(lam)>pi ) then
                        kind   = 2
                        lam    = pi
                        if ( l<0.0_wp ) lam = -lam
                        sinal  = 0.0_wp
                        cosal2 = 1.0_wp
                        test   = 2.0_wp
                        prev   = test
                        sig    = pi - abs(atan(sinu1/cosu1)+atan(sinu2/cosu2))
                        sinsig = sin(sig)
                        cossig = cos(sig)
                        c      = ((-3.0_wp*cosal2+4.0_wp)/rf+4.0_wp)*cosal2/rf/16.0_wp
                        if ( abs(sinal-prev)<tol ) exit longline
                        if ( abs(cosal2)<eps ) then
                            costm = -2.0_wp*(sinu1*sinu2/sign(eps,cosal2)) + cossig
                        else
                            costm = -2.0_wp*(sinu1*sinu2/cosal2) + cossig
                        endif
                        costm2 = costm*costm
                        cycle antipodal
                    endif
                    if ( ((lam-test)*(test-prev))<0.0_wp .and. it>5 ) &
                         lam = (2.0_wp*lam+3.0_wp*test+prev)/6.0_wp  ! refined converge.
                    prev = test
                    test = lam
                    cycle longline
                endif
            else
                sinal  = (lam-l)/d
                if ( ((sinal-test)*(test-prev))<0.0_wp .and. it>5 ) &
                       sinal = (2.0_wp*sinal+3.0_wp*test+prev)/6.0_wp  ! refined converge.
                prev   = test
                test   = sinal
                cosal2 = -sinal**2 + 1.0_wp
                sinlam = sinal*sinsig/(cosu1*cosu2)
                coslam = -sqrt(abs(-sinlam**2+1.0_wp))
                lam    = atan2(sinlam,coslam)
                temp   = cosu1*sinu2 - sinu1*cosu2*coslam
                sinsig = sqrt((cosu2*sinlam)**2+temp**2)
                cossig = sinu1*sinu2 + cosu1*cosu2*coslam
                sig    = atan2(sinsig,cossig)
                c      = ((-3.0_wp*cosal2+4.0_wp)/rf+4.0_wp)*cosal2/rf/16.0_wp
                if ( abs(sinal-prev)>=tol ) then
                    if ( abs(cosal2)<eps ) then
                        costm = -2.0_wp*(sinu1*sinu2/sign(eps,cosal2)) + cossig
                    else
                        costm = -2.0_wp*(sinu1*sinu2/cosal2) + cossig
                    endif
                    costm2 = costm*costm
                    cycle antipodal
                endif
            endif

            exit longline  !finished

        end do antipodal

    end do longline

    ! convergence

    if ( kind==2 ) then  ! antipodal
        faz  = sinal/cosu1
        baz  = sqrt(-faz**2+1.0_wp)
        if ( temp<0.0_wp ) baz = -baz
        faz  = atan2(faz,baz)
        tem1 = -sinal
        tem2 = sinu1*sinsig - cosu1*cossig*baz
        baz  = atan2(tem1,tem2)
    else  ! long-line
        tem1 = cosu2*sinlam
        tem2 = cosu1*sinu2 - sinu1*cosu2*coslam
        faz  = atan2(tem1,tem2)
        tem1 = -cosu1*sinlam
        tem2 = sinu1*cosu2 - cosu1*sinu2*coslam
        baz  = atan2(tem1,tem2)
    endif
    if ( faz<0.0_wp ) faz = faz + twopi
    if ( baz<0.0_wp ) baz = baz + twopi

    ! helmert 1880 from vincenty "geodetic inverse solution between antipodal points"

    ep2  = 1.0_wp/(boa*boa) - 1.0_wp
    bige = sqrt(1.0_wp+ep2*cosal2)
    bigf = (bige-1.0_wp)/(bige+1.0_wp)
    biga = (1.0_wp+bigf*bigf/4.0_wp)/(1.0_wp-bigf)
    bigb = bigf*(1.0_wp-0.375_wp*bigf*bigf)
    z    = bigb/6.0_wp*costm*(-3.0_wp+4.0_wp*sinsig**2)*(-3.0_wp+4.0_wp*costm2)
    dsig = bigb*sinsig*(costm+bigb/4.0_wp*(cossig*(-1.0_wp+2.0_wp*costm2)-z))
    s    = (boa*a)*biga*(sig-dsig)

    end subroutine inverse_vincenty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Function computes the Cartesian coordinates given the
!  geodetic latitude (phi), longitude (lambda) and
!  height (h) of a point related to an ellipsoid
!  defined by its three semiaxes ax, ay and b
!
!### History
!  * Jacob Williams, 10/29/2022 : Fortran verison of this algorithm,
!    based on the Matlab (v1.0 01/03/2019) code.

subroutine geodetic_to_cartesian_triaxial(ax, ay, b, phi, lambda, h, r)

    real(wp),intent(in) :: ax !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: ay !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: b  !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: phi !! geodetic latitude (radians)
    real(wp),intent(in) :: lambda !! geodetic longitude (radians)
    real(wp),intent(in) :: h !! geodetic height
    real(wp),dimension(3),intent(out) :: r  !! Cartesian position vector [x,y,z]

    real(wp) :: ee2,ex2,N,cp,sp,cl,sl

    cp  = cos(phi)
    sp  = sin(phi)
    cl  = cos(lambda)
    sl  = sin(lambda)
    ee2 = (ax*ax-ay*ay)/(ax*ax)
    ex2 = (ax*ax-b*b)/(ax*ax)
    N   = ax/sqrt(one - ex2*sp*sp - ee2*cp*cp*sl*sl)

    r = [(N+h)*cp*cl, &
         (N*(one-ee2)+h)*cp*sl, &
         (N*(one-ex2)+h)*sp ]

end subroutine geodetic_to_cartesian_triaxial

!********************************************************************************
!>
!  Geodetic to Cartesian for Triaxial Ellipsoid.
!
!### References
!  * S. Bektas, "Geodetic Computations on Triaxial Ellipsoid",
!    International Journal of Mining Science (IJMS),
!    Volume 1, Issue 1, June 2015, p 25-34

    pure subroutine geodetic_to_cartesian_triaxial_2(a,b,c,lat,long,h,r)

    implicit none

    real(wp),intent(in) :: a    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: b    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: c    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: lat  !! latitude (rad)
    real(wp),intent(in) :: long !! longitude (rad)
    real(wp),intent(in) :: h    !! altitude
    real(wp),dimension(3),intent(out) :: r  !! Cartesian coordinates (x,y,z)

    real(wp) :: ex2,ee2,v,a2,clat,slat,clon,slon,omee2,omex2

    a2    = a * a
    ex2   = (a2-c**2)/a2
    ee2   = (a2-b**2)/a2
    clat  = cos(lat)
    slat  = sin(lat)
    clon  = cos(long)
    slon  = sin(long)
    omee2 = 1.0_wp-ee2
    omex2 = 1.0_wp-ex2
    v     = a/sqrt(1.0_wp-ex2*slat**2-ee2*clat**2*slon**2)

    r = [(v+h)*clon*clat, &
         (v*omee2+h)*slon*clat, &
         (v*omex2+h)*slat ]

    end subroutine geodetic_to_cartesian_triaxial_2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Function computes the geodetic latitude (phi), longitude (lambda) and
!  height (h) of a point related to an ellipsoid
!  defined by its three semiaxes ax, ay and b (0 < b <= ay <= ax)
!  given Cartesian coordinates Xi, Yi, Zi and tolerance (tol).
!  Latitude and longitude are returned in radians.
!
!### Reference
!  * G. Panou and R. Korakitis, "Cartesian to geodetic coordinates conversion
!    on an ellipsoid using the bisection method".
!    Journal of Geodesy volume 96, Article number: 66 (2022).
!    [(link)](https://link.springer.com/article/10.1007/s00190-022-01650-9)
!  * [C++ code](https://www.researchgate.net/publication/353739609_PK-code)
!  * [MATLAB code](https://www.researchgate.net/publication/333904614_Cartesian2Geodetic_General_Panou_Korakitis)
!
!### History
!  * Jacob Williams, 10/29/2022 : Fortran verison of this algorithm.

subroutine cartesian_to_geodetic_triaxial(ax, ay, b, r, tol, phi, lambda, h)

    real(wp),intent(in) :: ax !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: ay !! semiaxes (0 < b <= ay <= ax)
    real(wp),intent(in) :: b  !! semiaxes (0 < b <= ay <= ax)
    real(wp),dimension(3),intent(in) :: r !! Cartesian coordinates (x,y,z)
    real(wp),intent(in) :: tol !! tolerance (may be set to zero)
    real(wp),intent(out) :: phi !! geodetic latitude (radians)
    real(wp),intent(out) :: lambda !! geodetic longitude (radians)
    real(wp),intent(out) :: h !! geodetic height

    real(wp) :: kx,ky,cx,cy,cz,XX,YY,ZZ,x,y,z,Xo,Yo,Zo,m,Mm,axax,ayay,b2
    integer :: n

    if (ax < ay .or. ay < b) error stop 'error in cartesian_to_geodetic_triaxial: invalid ax,ay,b'

    axax = ax*ax
    ayay = ay*ay
    b2   = b*b
    kx   = (axax-b2)/ax
    ky   = (ayay-b2)/ay
    cx   = (axax)/(b2)
    cy   = (ayay)/(b2)
    cz   = (axax)/(ayay)

    XX = abs(r(1))
    YY = abs(r(2))
    ZZ = abs(r(3))

    ! Compute geodetic latitude/longitude
    if (ZZ == zero) then
        if (XX == zero .and. YY == zero) then
            x = zero
            y = zero
            z = b
        else if (ky*XX*ky*XX+kx*YY*kx*YY < kx*ky*kx*ky) then
            x = ax*XX/kx
            y = ay*YY/ky
            z = b*sqrt(one-((x*x)/(axax))-((y*y)/(ayay)))
        else if (XX == zero) then
            x = zero
            y = ay
            z = zero
        else if (YY == zero) then
            x = ax
            y = zero
            z = zero
        else
            Xo = XX/ax
            Yo = YY/ay
            call bisection_special_2(cz, Xo, Yo, tol, n, m, Mm)
            x = cz*XX/(cz+m)
            y = YY/(one+m)
            z = zero
        end if
    else
        if (XX == zero .and. YY == zero) then
            x = zero
            y = zero
            z = b
        else
            Xo = XX/ax
            Yo = YY/ay
            Zo = ZZ/b
            call bisection_special_3(cx, cy, Xo, Yo, Zo, tol, n, m, Mm)
            x = cx*XX/(cx+m)
            y = cy*YY/(cy+m)
            if (m < zero .and. ky*XX*ky*XX + kx*YY*kx*YY < kx*ky*kx*ky) then
                z = b*sqrt(one-((x*x)/(axax))-((y*y)/(ayay)))
            else
                z = ZZ/(one+m)
            end if
        end if
    end if

    call xyz2fl(ax, ay, b, x, y, z, phi, lambda)        ! analytic method used for initial guess
    call xyz2philambda(ax, ay, b, x, y, z, phi, lambda) ! iterative method
    call philambda_quadrant(r(1), r(2), r(3), phi, lambda)

    ! Compute geodetic height
    h = norm2([XX-x, YY-y, ZZ-z])
    if ((XX+YY+ZZ) < (x+y+z)) h = -h

end subroutine cartesian_to_geodetic_triaxial

!*****************************************************************************************
!>

subroutine bisection_special_2(cz, Xo, Yo, tol, n, m, Gm)

    real(wp),intent(in) :: cz, Xo, Yo, tol
    integer,intent(out) :: n
    real(wp),intent(out) :: m, Gm

    real(wp) :: d1,Gd1,d2,d,MM

    d1 = -one+Yo
    Gd1 = (cz*Xo*cz*Xo)/((cz+d1)*(cz+d1))
    d2 = -one+sqrt(cz*Xo*cz*Xo+Yo*Yo)
    d = (d2 - d1)/two
    n = 0
    m = -two

    do while (d > tol)
        n = n + 1
        MM = m
        m = d1 + d
        Gm = ((cz*Xo*cz*Xo)/((cz+m)*(cz+m)))+((Yo*Yo)/((one+m)**2))-one
        if (MM == m + tol .or. Gm == zero) return
        if (sign(one,Gm) == sign(one,Gd1)) then
            d1 = m
            Gd1 = Gm
        else
            d2 = m
        end if
        d = (d2 - d1)/two
    end do

    n = n + 1
    m = d1 + d
    Gm = ((cz*Xo*cz*Xo)/((cz+m)*(cz+m)))+((Yo*Yo)/((one+m)**2))-one

end subroutine bisection_special_2

!*****************************************************************************************
!>

subroutine bisection_special_3(cx, cy, Xo, Yo, Zo, tol, n, m, Hm)

    real(wp),intent(in) :: cx, cy, Xo, Yo, Zo, tol
    integer,intent(out) :: n
    real(wp),intent(out) :: m, Hm

    real(wp) :: d1,Hd1,d2,d,MM

    d1 = -one+Zo
    Hd1 = ((cx*Xo*cx*Xo)/((cx+d1)*(cx+d1)))+((cy*Yo*cy*Yo)/((cy+d1)*(cy+d1)))
    d2 = -one+sqrt(cx*Xo*cx*Xo+cy*Yo*cy*Yo+Zo*Zo)
    d = (d2 - d1)/two
    n = 0
    m = -two

    do while (d > tol)
        n = n + 1
        MM = m
        m = d1 + d
        Hm = ((cx*Xo*cx*Xo)/((cx+m)*(cx+m)))+((cy*Yo*cy*Yo)/&
             ((cy+m)*(cy+m)))+((Zo*Zo)/((one+m)**2))-one
        if (MM == m + tol .or. Hm == zero) return
        if (sign(one,Hm) == sign(one,Hd1)) then
            d1 = m
            Hd1 = Hm
        else
            d2 = m
        end if
        d = (d2 - d1)/two
    end do

    n = n + 1
    m = d1 + d
    Hm = ((cx*Xo*cx*Xo)/((cx+m)*(cx+m)))+((cy*Yo*cy*Yo)/&
         ((cy+m)*(cy+m)))+((Zo*Zo)/((one+m)**2))-one

end subroutine bisection_special_3

!*****************************************************************************************
!>

subroutine philambda_quadrant(x, y, z, phi, lambda)

    real(wp),intent(in) :: x, y, z
    real(wp),intent(inout) :: phi, lambda

    if (z < zero) then
        phi = -phi
    end if

    if (x >= zero) then
        if (y >= zero) then
            lambda = lambda
        else
            lambda = -lambda
        end if
    else
        if (y >= zero) then
            lambda = pi-lambda
        else
            lambda = lambda-pi
        end if
    end if

end subroutine philambda_quadrant

!******************************************************************************
!>
!  Determination of the geodetic latitude and longitude
!
!@note This one has a different stopping criterion than the reference.

    subroutine xyz2philambda(ax, ay, b, x, y, z, phi, lambda)

     real(wp),intent(in) :: ax, ay, b, x, y, z
     real(wp),intent(inout) :: phi, lambda !! input: initial guess, output: refined values

     real(wp) :: ee2,ex2,Sphi,Cphi,Slambda,Clambda,&
                 Den,NN,onemee2,onemex2,dndphi,dxdphi,&
                 dydphi,dzdphi,dndlam,dxdlam,dydlam,dzdlam
     integer :: n
     real(wp),dimension(3,2) :: J
     real(wp),dimension(2,3) :: Jt !! transpose of J
     real(wp),dimension(3,1) :: dl
     real(wp),dimension(2,2) :: Nmat, Ninv
     real(wp),dimension(2,1) :: dx
     real(wp),dimension(3) :: r0

    ! real(wp) :: s0, SS0
    ! real(wp),dimension(3,1) :: UU
    ! real(wp),dimension(1,1) :: tmp

     integer,parameter :: maxiter = 100 !! maximum number of iterations
     real(wp),parameter :: stop_tol = 10.0_wp * epsilon(one) !! stopping tol for corrections

     ee2 = (ax*ax - ay*ay)/(ax*ax) ! eqn. 5
     ex2 = (ax*ax - b*b)/(ax*ax)   !
     onemee2 = one - ee2
     onemex2 = one - ex2

     !s0 = zero
     do n = 1, maxiter
        !SS0 = s0

        ! Design Matrix J
        Sphi = sin(phi)
        Cphi = cos(phi)
        Slambda = sin(lambda)
        Clambda = cos(lambda)

        NN  = ax/sqrt(one-ex2*Sphi*Sphi-ee2*Cphi*Cphi*Slambda*Slambda) ! eqn. 4
        Den = two*(one-ex2*Sphi**2-ee2*Cphi**2*Slambda**2)**(three/two)
        dndphi = -ax*sin(two*phi)*(ex2 - ee2*Slambda**2) / Den
        dxdphi = (dndphi*Cphi - NN*Sphi) * Clambda
        dydphi = onemee2*(dndphi*Cphi - NN*Sphi) * Slambda
        dzdphi = onemex2*(dndphi*Sphi + NN*Cphi)
        dndlam = -ax*ee2*Cphi**2*sin(two*lambda) / Den
        dxdlam = (dndlam*Clambda - NN*Slambda)*Cphi
        dydlam = onemee2*(dndlam*Slambda + NN*Clambda)*Cphi
        dzdlam = onemex2*dndlam*Sphi
        J = reshape([dxdphi,dydphi,dzdphi,dxdlam,dydlam,dzdlam],[3,2])

        ! Vector dl
        call geodetic_to_cartesian_triaxial_2(ax,ay,b,phi,lambda,0.0_wp,r0) ! just use the main one with alt=0
        dl(:,1) = [x,y,z] - r0 ! eqn. 51

        ! Solution
        Jt      = transpose(J)
        Nmat    = matmul(Jt,J) ! eqn. 53
        Ninv    = (one / (Nmat(1,1)*Nmat(2,2) - Nmat(1,2)*Nmat(2,1))) * &
                  reshape([Nmat(2,2),-Nmat(2,1),-Nmat(1,2),Nmat(1,1)], [2,2]) ! eqn. 54
        dx      = matmul(Ninv, matmul(Jt,dl)) ! eqn. 52
        phi     = phi    + dx(1,1)   ! corrections. eqn. 55
        lambda  = lambda + dx(2,1)   !

        ! ! original:
        ! UU      = matmul(J,dx) - dl
        ! tmp     = sqrt(matmul(transpose(UU),UU))
        ! s0      = tmp(1,1)
        ! if (s0 == SS0) exit

        ! JW: I think this is a better stopping criterion:
        if (all(abs(dx) <= stop_tol)) exit

     end do

    end subroutine xyz2philambda
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the transformation of Cartesian to geodetic coordinates on the surface of the ellipsoid
!  assuming x,y,z are all non-negative
!  Angular coordinates in radians
!
!  This is based on the [C++ version](https://www.researchgate.net/publication/353739609_PK-code)

subroutine xyz2fl(ax, ay, b, x, y, z, latitude, longitude)

    real(wp),intent(in) :: ax, ay, b, x, y, z
    real(wp),intent(out) :: latitude, longitude

    real(wp) :: nom,den,dex,xme,rot
    real(wp) :: ax2,ay2,b2,Ex2,Ee2,lex2,lee2,mex,mee

    ! note: these could be precomputed:
    ax2  = ax*ax
    ay2  = ay*ay
    b2   = b*b
    Ex2  = ax2-b2
    Ee2  = ax2-ay2
    lex2 = Ex2/ax2
    lee2 = Ee2/ax2
    mex  = one-lex2
    mee  = one-lee2

    nom = mee*z
    xme = mee*x
    dex = xme*xme+y*y
    den = mex*sqrt(dex)
    rot = sqrt(dex)

    if (den==zero)  then
        latitude = halfpi
        longitude = zero
    else
        if (nom<=den) then
            latitude = atan(nom/den)
        else
            latitude = halfpi-atan(den/nom)
        end if
        if (y<=xme) then
            den = xme+rot
            longitude = two*atan(y/den)
        else
            den = y+rot
            longitude = halfpi - two*atan(xme/den)
        end if
    end if

end subroutine xyz2fl

!****************************************************************
!>
!  Numerical solution to polynomial equation using Newton-Raphson method

pure function solve_polynomial(B, x0, error) result(x)

    real(wp),dimension(0:6),intent(in) :: B !! Polynomial `B = B(6) x^6 + B(5) x^5 + ... + B(0)`
    real(wp),intent(in) :: x0 !! Initial point
    real(wp),intent(in) :: error !! Maximum error
    real(wp) :: x !! root found after applying Newton-Raphson method to `B`
                  !! The function returns the value when the correction
                  !! is smaller than error.

    real(wp) :: f,fp,corr
    integer :: i, j !! counter

    integer,parameter :: maxiter = 100 !! maximum number of iterations

    x = x0
    do i = 1, maxiter
        f  = B(6)
        do j = 5, 0, -1
            if (j==5) then
                fp = f
            else
                fp = x*fp + f
            end if
            f  = x*f + B(j)
        end do
        if (fp==zero) exit ! singular point
        corr = f/fp
        x = x - corr
        if (abs(corr)<=error) exit
    end do

end function solve_polynomial

!****************************************************************
!>
!  Horner's method to compute `B(x-c)` in terms of `B(x)`.

pure subroutine horner(B, c, BB)

  real(wp),dimension(0:6),intent(in) :: B !! Polynomial `B = B(6) x^6 + B(5) x^5 + ... + B(0)`
  real(wp),intent(in) :: c
  real(wp),dimension(0:6),intent(out) :: BB !! Polynomial `BB` such that `B(x-c) = BB(x)`

  integer :: i,j !! counters

  BB = B

  do i = 0,6
    do j = 5,i,-1
      BB(j) = BB(j) - BB(j+1)*c
    end do
  end do

end subroutine horner

!******************************************************************
!>
!  Cartesian to Geodetic I
!
!### See also
!  * [[CartesianIntoGeodeticII]]
!
!### Reference
!  * Gema Maria Diaz-Toca, Leandro Marin, Ioana Necula,
!    "Direct transformation from Cartesian into geodetic coordinates on a triaxial ellipsoid"
!    Computers & Geosciences, Volume 142, September 2020, 104551.
!    [link](https://www.sciencedirect.com/science/article/pii/S0098300420305410?via%3Dihub),
!  * [C++ code](https://data.mendeley.com/datasets/s5f6sww86x/2) [CC BY 4.0 License]

subroutine CartesianIntoGeodeticI(ax, ay, az, r, latitude, longitude, altitude, error)

    real(wp),intent(in) :: ax, ay, az !! semiaxes of the celestial body: ax>ay>az
    real(wp),dimension(3),intent(in) :: r !! cartesian coordinates of the considered point
                                          !! in the first octant: xG, yG, zG with (xG,yG,zG)<>(0,0,0)
    real(wp),intent(out) :: latitude, longitude, altitude !! geodetic coordinates of the considered point
    real(wp),intent(in) :: error !! Values smaller than error treated as 0.0

    real(wp) :: ax2,ay2,az2,ax4,ay4,az4,b5,b4,b3,b3x,b3y,b3z,&
                b2,b2x,b2y,b2z,b1,b1x,b1y,b1z,b0,b0x,b0y,b0z,eec,exc
    real(wp) :: xg2,yg2,zg2,aux,xG,yG,zG
    real(wp) :: xE,yE,zE,k,B(0:6),BB(0:6)
    logical :: done

    call special_cases(ax,ay,az,r(1),r(2),r(3),latitude,longitude,altitude,done)
    if (done) return

    ! Computations independent of xG,yG,zG. They can be precomputed, if necessary.
    ax2 = ax*ax
    ay2 = ay*ay
    az2 = az*az
    ax4 = ax2*ax2
    ay4 = ay2*ay2
    az4 = az2*az2
    b5  = 2.0_wp*(ax2+ay2+az2)
    b4  = ax4 + 4.0_wp*ax2*ay2 + ay4 + 4.0_wp*ax2*az2 + 4.0_wp*ay2*az2 + az4
    b3  = 2.0_wp*ax4*ay2 + 2.0_wp*ax2*ay4 + 2.0_wp*ax4*az2 + 8.0_wp*ax2*ay2*az2 + 2.0_wp*ay4*az2 + 2.0_wp*ax2*az4 + 2.0_wp*ay2*az4
    b3x = - 2.0_wp*ax2*ay2 - 2.0_wp*ax2*az2
    b3y = - 2.0_wp*ax2*ay2 - 2.0_wp*ay2*az2
    b3z = - 2.0_wp*ay2*az2 - 2.0_wp*ax2*az2
    b2  = 4.0_wp*ax4*ay2*az2 + 4.0_wp*ax2*ay4*az2 + ax4*az4 + 4.0_wp*ax2*ay2*az4 + ax4*ay4 + ay4*az4
    b2x = -ax2*ay4 -4.0_wp*ax2*ay2*az2 -ax2*az4
    b2y = -ax4*ay2 -4.0_wp*ax2*ay2*az2 -ay2*az4
    b2z = -ax4*az2 -4.0_wp*ax2*ay2*az2 -ay4*az2
    b1  = 2.0_wp*ax4*ay4*az2 + 2.0_wp*ax4*ay2*az4 + 2.0_wp*ax2*ay4*az4
    b1x = - 2.0_wp*ax2*ay4*az2 - 2.0_wp*ax2*ay2*az4
    b1y = - 2.0_wp*ax4*ay2*az2 - 2.0_wp*ax2*ay2*az4
    b1z = - 2.0_wp*ax4*ay2*az2 - 2.0_wp*ax2*ay4*az2
    b0  = ax4*ay4*az4
    b0x = - ax2*ay4*az4
    b0y = - ax4*ay2*az4
    b0z = - ax4*ay4*az2
    eec = (ax2-ay2)/ax2
    exc = (ax2-az2)/ax2

    ! Computations dependant of xG, yG, zG
    xG = abs(r(1))
    yG = abs(r(2))
    zG = abs(r(3))
    xg2 = xG*xG
    yg2 = yG*yG
    zg2 = zG*zG
    aux = xg2/ax2+yg2/ay2+zg2/az2

    B = [b0+b0x*xg2+b0y*yg2+b0z*zg2, &
         b1+b1x*xg2+b1y*yg2+b1z*zg2, &
         b2+b2x*xg2+b2y*yg2+b2z*zg2, &
         b3+b3x*xg2+b3y*yg2+b3z*zg2, &
         b4-(ax2*xg2+ay2*yg2+az2*zg2), &
         b5, &
         1.0_wp ]

  if (abs(aux-1.0_wp) < error) then ! The point is on the ellipsoid

    xE = xG
    yE = yG
    zE = zG

  else if (aux > 1.0_wp) then ! The point is outside the ellipsoid

    k = solve_polynomial(B,(xg2+yg2+zg2)/3.0_wp,error)
    xE = ax2*xG/(ax2+k)
    yE = ay2*yG/(ay2+k)
    zE = az2*zG/(az2+k)

  else if (zG > 0.0_wp) then ! The point  is inside the ellipsoid and zG>0

    call horner(B,az2,BB) ! B(x-az2) = BB(x)
    k = solve_polynomial(BB,(xg2+yg2+zg2)/3.0_wp+az2,error) - az2
    xE = ax2*xG/(ax2+k)
    yE = ay2*yG/(ay2+k)
    zE = az2*zG/(az2+k)

  else if (xG > 0.0_wp .and. yG > 0.0_wp) then ! The point is inside the ellipsoid and zG=0, yG > 0, xG > 0

    call horner(B,ay2,BB)
    k = solve_polynomial(BB,(xg2+yg2+zg2)/3.0_wp+ay2,error) - ay2
    xE = ax2*xG/(ax2+k)
    yE = ay2*yG/(ay2+k)
    zE = 0.0_wp

  else if (xG < error .and. yG > 0.0_wp) then

    xE = 0.0_wp
    yE = ay
    zE = 0.0_wp

  else if (xG > 0.0_wp .and. yG < error) then

    xE = ax
    yE = 0.0_wp
    zE = 0.0_wp

  end if

  ! Computing longitude
  if (xG > 0.0_wp) then
    longitude = atan(yE/((1.0_wp-eec)*xE))
  else if (yG > 0.0_wp) then
    longitude = halfpi
  else
    longitude = huge(1.0_wp)  ! undefined
  end if

  ! Computing latitude
  if (xE > 0.0_wp .or. yE > 0.0_wp) then
    latitude = atan((1.0_wp-eec)/(1.0_wp-exc)*zE/norm2([xE*(1.0_wp-eec),yE]))
  else
    latitude = halfpi
  end if

  ! Computing altitude
  if (aux>=1.0_wp) then
    altitude = norm2([xE-xG,yE-yG,zE-zG])
  else
    altitude = -norm2([xE-xG,yE-yG,zE-zG])
  end if

  call philambda_quadrant(r(1), r(2), r(3), latitude, longitude)

  end subroutine CartesianIntoGeodeticI

!******************************************************************
!>
!  Cartesian into Geodetic II
!
!### See also
!  * [[CartesianIntoGeodeticI]]

subroutine CartesianIntoGeodeticII(ax, ay, az, r, latitude, longitude, altitude, error)

    real(wp),intent(in) :: ax, ay, az !! semiaxes of the celestial body: ax>ay>az
    real(wp),dimension(3),intent(in) :: r !! cartesian coordinates of the considered point
                                          !! in the first octant: xG, yG, zG with (xG,yG,zG)<>(0,0,0)
    real(wp),intent(out) :: latitude, longitude, altitude !! geodetic coordinates of the considered point
    real(wp),intent(in) :: error !! Values smaller than error treated as 0.0

    real(wp) :: aymaz,aypaz,axmaz,axpaz,axpaz2,ax2,ay2,az2,ax4,ay4,az4,az6,&
                az8,temp0,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,&
                temp9,tempa,az6ax2,az6ay2,tempb,maz10,excc,eecc
    real(wp) :: xg2,yg2,zg2,zgxg2,zgyg2,zg3,zg4,aux,xG,yG,zG
    real(wp) :: xE,yE,zE,k,B(0:6)
    logical :: done

    call special_cases(ax,ay,az,r(1),r(2),r(3),latitude,longitude,altitude,done)
    if (done) return

    ! Computations independent of xG,yG,zG. They can be precomputed, if necessary.
    aymaz = ay-az
    aypaz = ay+az
    axmaz = ax-az
    axpaz = ax+az
    axpaz2 = axpaz*axpaz
    ax2 = ax*ax
    ay2 = ay*ay
    az2 = az*az
    ax4 = ax2*ax2
    ay4 = ay2*ay2
    az4 = az2*az2
    az6 = az4*az2
    az8 = az4*az4
    temp0 = aymaz*aymaz*aypaz*aypaz*axmaz*axmaz*axpaz2
    temp1 = 2*az2*aymaz*aypaz*axmaz*axpaz*(ax2+ay2-2*az2)
    temp2 = -az2*(ax4*ay4-2*ax4*ay2*az2+ax4*az4-2*ax2*ay4*az2+4*ax2*ay2*az4-2*ay2*az6+az8-2*ax2*az6+ay4*az4)
    temp3 = -az2*(-ax2*ay4+2*ax2*ay2*az2-ax2*az4)
    temp4 = -az2*(-ax4*ay2+2*ax2*ay2*az2-ay2*az4)
    temp5 = -az2*(-ax4*az2-4*ax2*ay2*az2+6*ax2*az4-ay4*az2+6*ay2*az4-6*az6)
    temp6 =  -2*az4*(ax4*ay2-ax4*az2+ax2*ay4-4*ax2*ay2*az2+3*ax2*az4-ay4*az2+3*ay2*az4-2*az6)
    temp7 = -2*az4*(-ax2*ay2+ax2*az2)
    temp8 = -2*az4*(-ax2*ay2+ay2*az2)
    temp9 = -2*az4*(-ax2*az2-ay2*az2+2*az4)
    tempa = -az6*(ax4+4*ax2*ay2-6*ax2*az2+ay4-6*ay2*az2+6*az4)
    az6ax2 = az6*ax2
    az6ay2 = az6*ay2
    tempb = -2*az8*(ax2+ay2-2*az2)
    maz10 = -az6*az4
    excc = (ax2-az2)/(ax2)
    eecc = (ax2-ay2)/(ax2)

    xG = abs(r(1))
    yG = abs(r(2))
    zG = abs(r(3))
    xg2 = xG*xG
    yg2 = yG*yG
    zg2 = zG*zG
    zgxg2 = zG*xg2
    zgyg2 = zG*yg2
    zg3 = zg2*zG
    zg4 = zg2*zg2
    aux = xg2/ax2+yg2/ay2+zg2/az2

  if (abs(aux-1.0_wp) < error) then ! The point is on the ellipsoid

    xE = xG
    yE = yG
    zE = zG

  else if (zG > error) then ! The point is inside or outside the ellipsoid with zG != 0

    B(6) = temp0
    B(5) = temp1*zG
    B(4) = temp2+temp3*xg2+temp4*yg2+temp5*zg2
    B(3) = zG*temp6+temp7*zgxg2+temp8*zgyg2+temp9*zg3
    B(2) = zg2*(tempa+az6ax2*xg2+az6ay2*yg2+az8*zg2)
    B(1) = tempb*zg3
    B(0) = maz10*zg4

    k = solve_polynomial(B,az*zG/norm2([xG,yG,zG]),error)
    xE = ax2*xG*k/(ax2*k-az2*k+az2*zG)
    yE = ay2*yG*k/(ay2*k-az2*k+az2*zG)
    zE = k

  else if (yG > error) then

    B = [-ay4*ay2*zg2, &
         -2*ay4*(ax2-ay2)*zG, &
         -ay2*(ax4-2*ax2*ay2-ax2*yg2+ay4-ay2*zg2), &
         2*ay2*(ax2-ay2)*zG, &
         ax4+ay4-2*ax2*ay2, &
         0.0_wp, &
         0.0_wp ]

    k = solve_polynomial(B,ay*yG/norm2([xG,yG,zG]),error)
    xE = k*ax2*xG/(ax2*k-ay2*k+ay2*yG)
    yE = k
    zE = 0.0_wp

  else

    xE = ax
    yE = 0.0_wp
    zE = 0.0_wp

  end if

  ! Computing longitude

  if (xG > 0.0_wp) then
    longitude = atan(yE/((1.0_wp-eecc)*xE))
  else if (yG > 0.0_wp) then
    longitude = halfpi
  else
    longitude = huge(1.0_wp) ! undefined
  end if

  ! Computing latitude

  if (xE > 0.0_wp .or. yE > 0.0_wp) then
    latitude = atan((1.0_wp-eecc)/(1.0_wp-excc)*zE/norm2([xE*(1.0_wp-eecc),yE]))
  else
    latitude = halfpi
  end if

  ! Computing altitude

  if (aux>=1.0) then
    altitude = norm2([xE-xG,yE-yG,zE-zG])
  else
    altitude = -norm2([xE-xG,yE-yG,zE-zG])
  end if

  call philambda_quadrant(r(1), r(2), r(3), latitude, longitude)

end subroutine CartesianIntoGeodeticII

!********************************************************************************
!>
!  Cartesian to geodetic for Triaxial Ellipsoid.
!
!### References
!  * S. Bektas, "Geodetic Computations on Triaxial Ellipsoid",
!    International Journal of Mining Science (IJMS),
!    Volume 1, Issue 1, June 2015, p 25-34

    subroutine cartesian_to_geodetic_triaxial_2(a,b,c,r,eps,phi,lambda,h)

    implicit none

    real(wp),intent(in) :: a    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: b    !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in) :: c    !! ellipsoid radii `a >= b >= c`
    real(wp),dimension(3),intent(in) :: r !! Cartesian coordinates (x,y,z)
    real(wp),intent(in)  :: eps  !! convergence tolerance
    real(wp),intent(out) :: phi !! latitude (rad)
    real(wp),intent(out) :: lambda !! longitude (rad)
    real(wp),intent(out) :: h !! altitude

    integer,parameter :: maxiter = 20 !! maximum number of iterations

    integer :: i  !! iteration counter
    real(wp),dimension(3,3) :: AA
    real(wp),dimension(3) :: bvec, xvec
    real(wp) :: a2,b2,c2,x,y,z,ex2,ee2,e,f,g,xo,yo,zo,j11,j12,j21,j23,rmag,omee2
    logical :: success

    x = r(1)
    y = r(2)
    z = r(3)

    if (a<b .or. b<c) error stop 'error in cartesian_to_geodetic_triaxial_2: invalid a,b,c'
    call special_cases(a,b,c,x,y,z,phi,lambda,h,success)
    if (success) return

    rmag  = norm2(r)
    a2    = a*a
    b2    = b*b
    c2    = c*c
    ex2   = (a2-c2)/a2
    ee2   = (a2-b2)/a2
    omee2 = one-ee2
    E     = one/a2
    F     = one/b2
    G     = one/c2
    xo    = a*x/rmag
    yo    = b*y/rmag
    zo    = c*z/rmag

    do i = 1, maxiter

        j11 = F*yo-(yo-y)*E
        j12 = (xo-x)*F-E*xo
        j21 = G*zo-(zo-z)*E
        j23 = (xo-x)*G-E*xo

        ! solve the linear system:
        AA = reshape(-[j11,j21,two*E*xo,&
                       j12,zero,two*F*yo,&
                       zero,j23,two*G*zo], [3,3])
        bvec = [ (xo-x)*F*yo-(yo-y)*E*xo, &
                 (xo-x)*G*zo-(zo-z)*E*xo, &
                 E*xo**2+F*yo**2+G*zo**2-one ]
        call linear_solver(AA,bvec,xvec,success)
        if (.not. success) then
            write(*,*) 'error in cartesian_to_geodetic_triaxial_2: matrix is singular'
            phi    = zero
            lambda = zero
            h      = zero
            return
        end if
        xo = xo + xvec(1)
        yo = yo + xvec(2)
        zo = zo + xvec(3)

        if (maxval(abs(xvec))<eps) exit

    end do

    ! outputs:
    phi = atan(zo*omee2/(one-ex2)/sqrt(omee2**2*xo**2+yo**2))
    lambda = atan2(yo, omee2*xo)
    h = sign(one,z-zo)*sign(one,zo)*sqrt((x-xo)**2+(y-yo)**2+(z-zo)**2)

    contains

    subroutine linear_solver(a,b,x,success)

        !!  Solve the 3x3 system: `A * x = b`
        !!  Reference: https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt

        implicit none

        real(wp),dimension(3,3),intent(in) :: a
        real(wp),dimension(3),intent(in)   :: b
        real(wp),dimension(3),intent(out)  :: x
        logical,intent(out)                :: success

        real(wp) :: det !! determinant of a
        real(wp),dimension(3,3) :: adj !! adjoint of a
        real(wp),dimension(3,3) :: ainv !! inverse of a

        det =   a(1,1)*a(2,2)*a(3,3)  &
              - a(1,1)*a(2,3)*a(3,2)  &
              - a(1,2)*a(2,1)*a(3,3)  &
              + a(1,2)*a(2,3)*a(3,1)  &
              + a(1,3)*a(2,1)*a(3,2)  &
              - a(1,3)*a(2,2)*a(3,1)

        success = abs(det) > dblmin ! check for singularity

        if (success) then
            adj(:,1) = [a(2,2)*a(3,3)-a(2,3)*a(3,2),&
                        a(2,3)*a(3,1)-a(2,1)*a(3,3),&
                        a(2,1)*a(3,2)-a(2,2)*a(3,1)]
            adj(:,2) = [a(1,3)*a(3,2)-a(1,2)*a(3,3),&
                        a(1,1)*a(3,3)-a(1,3)*a(3,1),&
                        a(1,2)*a(3,1)-a(1,1)*a(3,2)]
            adj(:,3) = [a(1,2)*a(2,3)-a(1,3)*a(2,2),&
                        a(1,3)*a(2,1)-a(1,1)*a(2,3),&
                        a(1,1)*a(2,2)-a(1,2)*a(2,1)]
            ainv = adj/det
            x = matmul(ainv,b)
        else
            x = zero
        end if

        end subroutine linear_solver

    end subroutine cartesian_to_geodetic_triaxial_2
!********************************************************************************

!********************************************************************************
!>
!  Special cases for lat/lon/altitude

    subroutine special_cases(a,b,c,x,y,z,phi,lambda,h,done)

    real(wp),intent(in)  :: a      !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in)  :: b      !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in)  :: c      !! ellipsoid radii `a >= b >= c`
    real(wp),intent(in)  :: x      !! Cartesian x coordinate
    real(wp),intent(in)  :: y      !! Cartesian y coordinate
    real(wp),intent(in)  :: z      !! Cartesian z coordinate
    real(wp),intent(out) :: phi    !! latitude (rad)
    real(wp),intent(out) :: lambda !! longitude (rad)
    real(wp),intent(out) :: h      !! altitude
    logical,intent(out)  :: done   !! true if one of the special cases was computed

    logical :: x0, y0, z0

    real(wp),parameter :: zero_tol  = 10.0_wp * epsilon(1.0_wp)  !! zero tolerance for singularities

    x0 = abs(x) <= zero_tol
    y0 = abs(y) <= zero_tol
    z0 = abs(z) <= zero_tol

    if (x0 .and. y0 .and. z0) then ! center of the body
        phi    = zero
        lambda = zero
        h      = -c    ! just pick this value
        done = .true.
        return
    else if (x0 .and. y0) then ! (on the z-axis)
        if (z>=zero) then
            phi    = halfpi
            lambda = zero
            h      = z-c
        else
            phi    = -halfpi
            lambda = zero
            h      = -(z+c)
        end if
        done = .true.
        return
    else if (x0 .and. z0) then  ! on the y-axis
        if (y>=zero) then
            phi    = zero
            lambda = halfpi
            h      = y-b
        else
            phi    = zero
            lambda = -halfpi
            h      = -(y+b)
        end if
        done = .true.
        return
    else if (y0 .and. z0) then  ! on the x-axis
        if (x>=zero) then
            phi    = zero
            lambda = zero
            h      = x-a
        else
            phi    = zero
            lambda = pi
            h      = -(x+a)
        end if
        done = .true.
        return
    end if

    phi    = zero
    lambda = zero
    h      = zero
    done = .false.

    end subroutine special_cases
!*****************************************************************************************

  end module geodesic_module
