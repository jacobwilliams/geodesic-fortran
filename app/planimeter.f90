!> @file planimeter.for
!! @brief A test program for area()

!> A simple program to compute the area of a geodesic polygon.
!!
!! This program reads in up to 10000 lines with lat, lon for each vertex
!! of a polygon.  At the end of input, the program prints the number of
!! vertices, the perimeter of the polygon and its area (for the WGS84
!! ellipsoid).

program planimeter
    use geodesic_module, wp => geodesic_wp
implicit none

integer maxpts
parameter (maxpts = 10000)
real(wp) a, f, lats(maxpts), lons(maxpts), S, P
integer n

! WGS84 values
a = 6378137.0_wp
f = 1.0_wp/298.257223563_wp

n = 0
10 continue
if (n >= maxpts) go to 20
read(*, *, end=20, err=20) lats(n+1), lons(n+1)
n = n+1
go to 10
20 continue
call area(a, f, lats, lons, n, S, P)
print 30, n, P, S
30 format(i6, 1x, f20.8, 1x, f20.3)
stop
end
