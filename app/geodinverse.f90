!> @file geodinverse.for
!! @brief A test program for inverse()

!> A simple program to solve the inverse geodesic problem.
!!
!! This program reads in lines with lat1, lon1, lon2, lat2 and prints
!! out lines with azi1, azi2, s12 (for the WGS84 ellipsoid).

program geodinverse
    use geodesic_module, wp => geodesic_wp
implicit none

real(wp) a, f, lat1, lon1, azi1, lat2, lon2, azi2, s12, &
    dummy1, dummy2, dummy3, dummy4, dummy5
integer outmask

! WGS84 values
a = 6378137.0_wp
f = 1.0_wp/298.257223563_wp

outmask = 0

10 continue
read(*, *, end=90, err=90) lat1, lon1, lat2, lon2
call inverse(a, f, lat1, lon1, lat2, lon2, &
    s12, azi1, azi2, outmask, &
    dummy1, dummy2, dummy3, dummy4, dummy5)
print 20, azi1, azi2, s12
20 format(1x, f20.15, 1x, f20.15, 1x, f19.10)
go to 10
90 continue

stop
end
