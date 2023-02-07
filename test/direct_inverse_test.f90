
!*****************************************************************************************
!>
!  Unit test for the [[direct_vincenty]] and [[inverse_vincenty]] geodetic routines.
!
!@note This was originally from the
!      [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit)

    program direct_inverse_test

    use geodesic_module, wp => geodesic_wp

    implicit none

    !Ellipsoid : GRS80 / WGS84 (NAD83)
    real(wp),parameter :: a  = 6378137.0000_wp      !! Equatorial radius
    real(wp),parameter :: rf = 298.25722210088_wp   !! Inverse flattening
    real(wp),parameter :: f  = 1.0_wp / rf          !! flattening
    real(wp),parameter :: pi = acos(-1.0_wp)
    real(wp),parameter :: rad2deg = 180.0_wp / pi
    integer,parameter :: n_repeat = 1000  !! number of times to repeat speed test
    real(wp),parameter :: fail_tol = 0.1_wp !! large failure tol
    real(wp),parameter :: tol = 1.0e-8_wp !! test success tol

    real(wp) :: glat1,glon1,glat2,glon2
    integer,dimension(:),allocatable :: iseed !! for random number generator
    real(wp),dimension(3) :: errors1, errors2, maxerrors1, maxerrors2
    integer :: isize, i

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' direct_inverse_test'
    write(*,*) '---------------'
    write(*,*) ''

    call random_seed(size=isize)
    allocate(iseed(isize)); iseed = 42
    maxerrors1 = 0.0_wp
    maxerrors2 = 0.0_wp

    do i = 1, n_repeat

        ! specify two points (rad):
        glat1 = get_random_number(-pi/2,pi/2)
        glon1 = get_random_number(-pi,pi)
        glat2 = get_random_number(-pi/2,pi/2)
        glon2 = get_random_number(-pi,pi)

        ! call both routines:
        call vincenty_errors(glat1,glon1,glat2,glon2,errors1)
        call geodeticlib_errors(glat1,glon1,glat2,glon2,errors2)

        ! quick check for failures:
        if (any(abs(errors1) > fail_tol)) then
            write(*,*) 'glat1 = ', glat1
            write(*,*) 'glon1 = ', glon1
            write(*,*) 'glat2 = ', glat2
            write(*,*) 'glon2 = ', glon2
            write(*,*) 'errors = ', errors1
            error stop 'vincenty fail'
        end if
        if (any(abs(errors2) > fail_tol)) then
            write(*,*) 'glat1 = ', glat1
            write(*,*) 'glon1 = ', glon1
            write(*,*) 'glat2 = ', glat2
            write(*,*) 'glon2 = ', glon2
            write(*,*) 'errors = ', errors2
            error stop 'geodeticlib fail'
        end if

        ! save the max error
        where(abs(errors1)>maxerrors1) maxerrors1 = abs(errors1)
        where(abs(errors2)>maxerrors2) maxerrors2 = abs(errors2)

    end do

    write(*,*) ' max vincenty errors:    ', maxerrors1
    write(*,*) ' max geodeticlib errors: ', maxerrors2

    if (any(maxerrors1>tol .or. maxerrors2>tol)) error stop 'test failed'

    contains

    subroutine vincenty_errors(lat1,lon1,lat2,lon2,errors)
        real(wp),intent(in) :: lat1,lon1,lat2,lon2 ! in radians
        real(wp),dimension(3),intent(out) :: errors

        real(wp) :: faz,baz,s,sig,lam,glat2_,glon2_,baz_,e
        integer :: it, kind

        call inverse_vincenty(a,rf,lat1,lon1,lat2,lon2,faz,baz,s,it,sig,lam,kind)
        call direct_vincenty(a,f,lat1,lon1,faz,s,glat2_,glon2_,baz_)

        errors(1) = AngDiff(glat2_*rad2deg , lat2*rad2deg , e)
        errors(2) = AngDiff(glon2_*rad2deg , lon2*rad2deg , e)
        errors(3) = AngDiff(baz_  *rad2deg , baz*rad2deg  , e)

    end subroutine vincenty_errors

    subroutine geodeticlib_errors(lat1,lon1,lat2,lon2,errors)
        real(wp),intent(in) :: lat1,lon1,lat2,lon2 ! in radians
        real(wp),dimension(3),intent(out) :: errors

        real(wp) :: faz,baz,s,sig,lam,glat2_,glon2_,baz_
        real(wp) :: a12, m12, MM12, MM21, SS12, a12s12, e
        integer :: flags, outmask

        outmask = 0
        flags = 0
        call inverse(a, f, lat1*rad2deg, lon1*rad2deg, lat2*rad2deg, lon2*rad2deg, &
                     s, faz, baz, outmask, a12, m12, MM12, MM21, SS12)
        call direct(a, f, lat1*rad2deg, lon1*rad2deg, faz, s, flags, &
                     glat2_,glon2_,baz_, outmask, a12s12, m12, MM12, MM21, SS12)

        errors(1) = AngDiff(glat2_ , lat2*rad2deg, e)
        errors(2) = AngDiff(glon2_ , lon2*rad2deg, e)
        errors(3) = AngDiff(baz_ , baz, e)

    end subroutine geodeticlib_errors

    !*****************************************************************************************
    !> author: Jacob Williams
    !
    !  Returns a uniform random number `x`, such that: `a <= x < b`.

    function get_random_number(a,b) result(x)

        implicit none

        real(wp)            :: x
        real(wp),intent(in) :: a
        real(wp),intent(in) :: b

        call random_number(x)

        x = a + (b-a)*x

        end function get_random_number
    !*****************************************************************************************

    end program direct_inverse_test
!*****************************************************************************************
