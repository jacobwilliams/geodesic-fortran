!> @file geodtest.for
!! @brief Test suite for the geodesic routines in Fortran
!!
!! Run these tests by configuring with cmake and running "make test".
!!
!! Copyright (c) Charles Karney (2015-2022) <charles@karney.com> and
!! licensed under the MIT/X11 License.  For more information, see
!! https://geographiclib.sourceforge.io/

program geodtest
    use tests_module
    use geodesic_module, wp => geodesic_wp

    implicit none

    integer :: n, i

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' geodtest'
    write(*,*) '---------------'
    write(*,*) ''

    n = 0
    i = tstinv()
    if (i > 0) then
        n = n + 1
        print *, 'tstinv fail:', i
    end if
    i = tstdir()
    if (i > 0) then
        n = n + 1
        print *, 'tstdir fail:', i
    end if
    i = tstarc()
    if (i > 0) then
        n = n + 1
        print *, 'tstarc fail:', i
    end if
    i = tstg0()
    if (i > 0) then
        n = n + 1
        print *, 'tstg0 fail:', i
    end if
    i = tstg1()
    if (i > 0) then
        n = n + 1
        print *, 'tstg1 fail:', i
    end if
    i = tstg2()
    if (i > 0) then
        n = n + 1
        print *, 'tstg2 fail:', i
    end if
    i = tstg5()
    if (i > 0) then
        n = n + 1
        print *, 'tstg5 fail:', i
    end if
    i = tstg6()
    if (i > 0) then
        n = n + 1
        print *, 'tstg6 fail:', i
    end if
    i = tstg9()
    if (i > 0) then
        n = n + 1
        print *, 'tstg9 fail:', i
    end if
    i = tstg10()
    if (i > 0) then
        n = n + 1
        print *, 'tstg10 fail:', i
    end if
    i = tstg11()
    if (i > 0) then
        n = n + 1
        print *, 'tstg11 fail:', i
    end if
    i = tstg12()
    if (i > 0) then
        n = n + 1
        print *, 'tstg12 fail:', i
    end if
    i = tstg14()
    if (i > 0) then
        n = n + 1
        print *, 'tstg14 fail:', i
    end if
    i = tstg15()
    if (i > 0) then
        n = n + 1
        print *, 'tstg15 fail:', i
    end if
    i = tstg17()
    if (i > 0) then
        n = n + 1
        print *, 'tstg17 fail:', i
    end if
    i = tstg26()
    if (i > 0) then
        n = n + 1
        print *, 'tstg26 fail:', i
    end if
    i = tstg28()
    if (i > 0) then
        n = n + 1
        print *, 'tstg28 fail:', i
    end if
    i = tstg33()
    if (i > 0) then
        n = n + 1
        print *, 'tstg33 fail:', i
    end if
    i = tstg55()
    if (i > 0) then
        n = n + 1
        print *, 'tstg55 fail:', i
    end if
    i = tstg59()
    if (i > 0) then
        n = n + 1
        print *, 'tstg59 fail:', i
    end if
    i = tstg61()
    if (i > 0) then
        n = n + 1
        print *, 'tstg61 fail:', i
    end if
    i = tstg73()
    if (i > 0) then
        n = n + 1
        print *, 'tstg73 fail:', i
    end if
    i = tstg74()
    if (i > 0) then
        n = n + 1
        print *, 'tstg74 fail:', i
    end if
    i = tstg76()
    if (i > 0) then
        n = n + 1
        print *, 'tstg76 fail:', i
    end if
    i = tstg78()
    if (i > 0) then
        n = n + 1
        print *, 'tstg78 fail:', i
    end if
    i = tstg80()
    if (i > 0) then
        n = n + 1
        print *, 'tstg80 fail:', i
    end if
    i = tstg84()
    if (i > 0) then
        n = n + 1
        print *, 'tstg84 fail:', i
    end if
    i = tstg92()
    if (i > 0) then
        n = n + 1
        print *, 'tstg92 fail:', i
    end if
    i = tstg94()
    if (i > 0) then
        n = n + 1
        print *, 'tstg94 fail:', i
    end if
    i = tstg96()
    if (i > 0) then
        n = n + 1
        print *, 'tstg96 fail:', i
    end if
    i = tstp0()
    if (i > 0) then
        n = n + 1
        print *, 'tstp0 fail:', i
    end if
    i = tstp5()
    if (i > 0) then
        n = n + 1
        print *, 'tstp5 fail:', i
    end if
    i = tstp6()
    if (i > 0) then
        n = n + 1
        print *, 'tstp6 fail:', i
    end if
    i = tstp12()
    if (i > 0) then
        n = n + 1
        print *, 'tstp12 fail:', i
    end if
    i = tstp12r()
    if (i > 0) then
        n = n + 1
        print *, 'tstp12r fail:', i
    end if
    i = tstp13()
    if (i > 0) then
        n = n + 1
        print *, 'tstp13 fail:', i
    end if
    i = tstp15()
    if (i > 0) then
        n = n + 1
        print *, 'tstp15 fail:', i
    end if
    i = tstp19()
    if (i > 0) then
        n = n + 1
        print *, 'tstp19 fail:', i
    end if
    i = tstp21()
    if (i > 0) then
        n = n + 1
        print *, 'tstp21 fail:', i
    end if

    if (n > 0) then
        error stop
    else
        write(*,*) 'All tests passed'
    end if

contains

    ! integer function assert(x, y, d)
    ! real(wp),intent(in) :: x, y, d

    ! if (abs(x - y) <= d) then
    !     assert = 0
    ! else
    !     assert = 1
    !     write(*,'(1x,a,g14.7,a, g14.7,a,g10.3)') &
    !             'assert fails: ', x, ' != ', y, ' +/- ', d
    ! end if

    ! end function assert

    ! integer function chknan(x)
    ! real(wp),intent(in) :: x

    ! if (x /= x) then
    !     chknan = 0
    ! else
    !     chknan = 1
    ! end if

    ! end function chknan

    integer function tstinv()

    real(wp) lat1, lon1, azi1, lat2, lon2, azi2, &
        s12, a12, m12, MM12, MM21, SS12
    real(wp) azi1a, azi2a, s12a, a12a, &
        m12a, MM12a, MM21a, SS12a
    real(wp) a, f
    integer r, i, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 1 + 2 + 4 + 8
    r = 0

    do i = 1,20
        lat1 = tstdat(i, 1)
        lon1 = tstdat(i, 2)
        azi1 = tstdat(i, 3)
        lat2 = tstdat(i, 4)
        lon2 = tstdat(i, 5)
        azi2 = tstdat(i, 6)
        s12 = tstdat(i, 7)
        a12 = tstdat(i, 8)
        m12 = tstdat(i, 9)
        MM12 = tstdat(i, 10)
        MM21 = tstdat(i, 11)
        SS12 = tstdat(i, 12)
        call inverse(a, f, lat1, lon1, lat2, lon2, &
            s12a, azi1a, azi2a, outmask, a12a, m12a, MM12a, MM21a, SS12a)
        r = r + assert(azi1, azi1a, 1.0e-13_wp)
        r = r + assert(azi2, azi2a, 1.0e-13_wp)
        r = r + assert(s12,  s12a,  1.0e-8_wp)
        r = r + assert(a12,  a12a,  1.0e-13_wp)
        r = r + assert(m12,  m12a,  1.0e-8_wp)
        r = r + assert(MM12, MM12a, 1.0e-15_wp)
        r = r + assert(MM21, MM21a, 1.0e-15_wp)
        r = r + assert(SS12, SS12a, 0.1_wp)
    end do

    tstinv = r

    end function tstinv

    integer function tstdir()

    real(wp) lat1, lon1, azi1, lat2, lon2, azi2, &
        s12, a12, m12, MM12, MM21, SS12
    real(wp) lat2a, lon2a, azi2a, a12a, &
        m12a, MM12a, MM21a, SS12a
    real(wp) a, f
    integer r, i, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 1 + 2 + 4 + 8
    flags = 2
    r = 0

    do i = 1,20
        lat1 = tstdat(i, 1)
        lon1 = tstdat(i, 2)
        azi1 = tstdat(i, 3)
        lat2 = tstdat(i, 4)
        lon2 = tstdat(i, 5)
        azi2 = tstdat(i, 6)
        s12 = tstdat(i, 7)
        a12 = tstdat(i, 8)
        m12 = tstdat(i, 9)
        MM12 = tstdat(i, 10)
        MM21 = tstdat(i, 11)
        SS12 = tstdat(i, 12)
        call direct(a, f, lat1, lon1, azi1, s12, flags, &
            lat2a, lon2a, azi2a, outmask, a12a, m12a, MM12a, MM21a, SS12a)
        r = r + assert(lat2, lat2a, 1.0e-13_wp)
        r = r + assert(lon2, lon2a, 1.0e-13_wp)
        r = r + assert(azi2, azi2a, 1.0e-13_wp)
        r = r + assert(a12,  a12a,  1.0e-13_wp)
        r = r + assert(m12,  m12a,  1.0e-8_wp )
        r = r + assert(MM12, MM12a, 1.0e-15_wp)
        r = r + assert(MM21, MM21a, 1.0e-15_wp)
        r = r + assert(SS12, SS12a, 0.1_wp)
    end do

    tstdir = r

    end function tstdir

    integer function tstarc()

    real(wp) lat1, lon1, azi1, lat2, lon2, azi2, &
        s12, a12, m12, MM12, MM21, SS12
    real(wp) lat2a, lon2a, azi2a, s12a, &
        m12a, MM12a, MM21a, SS12a
    real(wp) a, f
    integer r, i, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 1 + 2 + 4 + 8
    flags = 1 + 2
    r = 0

    do i = 1,20
        lat1 = tstdat(i, 1)
        lon1 = tstdat(i, 2)
        azi1 = tstdat(i, 3)
        lat2 = tstdat(i, 4)
        lon2 = tstdat(i, 5)
        azi2 = tstdat(i, 6)
        s12 = tstdat(i, 7)
        a12 = tstdat(i, 8)
        m12 = tstdat(i, 9)
        MM12 = tstdat(i, 10)
        MM21 = tstdat(i, 11)
        SS12 = tstdat(i, 12)
        call direct(a, f, lat1, lon1, azi1, a12, flags, &
            lat2a, lon2a, azi2a, outmask, s12a, m12a, MM12a, MM21a, SS12a)
        r = r + assert(lat2, lat2a, 1.0e-13_wp)
        r = r + assert(lon2, lon2a, 1.0e-13_wp)
        r = r + assert(azi2, azi2a, 1.0e-13_wp)
        r = r + assert(s12,  s12a,  1.0e-8_wp )
        r = r + assert(m12,  m12a,  1.0e-8_wp )
        r = r + assert(MM12, MM12a, 1.0e-15_wp)
        r = r + assert(MM21, MM21a, 1.0e-15_wp)
        r = r + assert(SS12, SS12a, 0.1_wp)
    end do

    tstarc = r

    end function tstarc

    integer function tstg0()
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 40.6_wp, -73.8_wp, 49.01666667_wp, 2.55_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 53.47022_wp,  0.5e-5_wp)
    r = r + assert(azi2, 111.59367_wp, 0.5e-5_wp)
    r = r + assert(s12,  5853226.0_wp, 0.5_wp)

    tstg0 = r

    end function tstg0

    integer function tstg1()
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    flags = 0
    r = 0
    call direct(a, f, 40.63972222_wp, -73.77888889_wp, 53.5_wp, 5850.0e3_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(lat2, 49.01467_wp,  0.5e-5_wp)
    r = r + assert(lon2, 2.56106_wp,   0.5e-5_wp)
    r = r + assert(azi2, 111.62947_wp, 0.5e-5_wp)

    tstg1 = r

    end function tstg1

    integer function tstg2()
    ! Check fix for antipodal prolate bug found 2010-09-04
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    a = 6.4e6_wp
    f = -1.0_wp/150.0_wp
    outmask = 0
    r = 0
    call inverse(a, f, 0.07476_wp, 0.0_wp, -0.07476_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 90.00078_wp,   0.5e-5_wp)
    r = r + assert(azi2, 90.00078_wp,   0.5e-5_wp)
    r = r + assert(s12,  20106193.0_wp, 0.5_wp)
    call inverse(a, f, 0.1_wp, 0.0_wp, -0.1_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 90.00105_wp, 0.5e-5_wp)
    r = r + assert(azi2, 90.00105_wp, 0.5e-5_wp)
    r = r + assert(s12,  20106193.0_wp, 0.5_wp)

    tstg2 = r

    end function tstg2

    integer function tstg4()
    ! Check fix for short line bug found 2010-05-21
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, &
        36.493349428792_wp, 0.0_wp, 36.49334942879201_wp, 0.0000008_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(s12, 0.072_wp, 0.5e-3_wp)

    tstg4 = r

    end function tstg4

    integer function tstg5()
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    flags = 0
    r = 0
    call direct(a, f, 0.01777745589997_wp, 30.0_wp, 0.0_wp, 10.0e6_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    if (lon2 < 0) then
        r = r + assert(lon2, -150.0_wp, 0.5e-5_wp)
        r = r + assert(abs(azi2), 180.0_wp, 0.5e-5_wp)
    else
        r = r + assert(lon2, 30.0_wp, 0.5e-5_wp)
        r = r + assert(azi2, 0.0_wp, 0.5e-5_wp)
    end if

    tstg5 = r

    end function tstg5

    integer function tstg6()
    ! Check fix for volatile sbet12a bug found 2011-06-25 (gcc 4.4d0.4d0
    ! x86 -O3).  Found again on 2012-03-27 with tdm-mingw32 (g++ 4.6d0.1d0).
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 88.202499451857_wp, 0.0_wp, &
        -88.202499451857_wp, 179.981022032992859592_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(s12, 20003898.214_wp, 0.5e-3_wp)
    call inverse(a, f, 89.262080389218_wp, 0.0_wp, &
        -89.262080389218_wp, 179.992207982775375662_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(s12, 20003925.854_wp, 0.5e-3_wp)
    call inverse(a, f, 89.333123580033_wp, 0.0_wp, &
        -89.333123580032997687_wp, 179.99295812360148422_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(s12, 20003926.881_wp, 0.5e-3_wp)

    tstg6 = r

    end function tstg6

    integer function tstg9()
    ! Check fix for volatile x bug found 2011-06-25 (gcc 4.4d0.4d0 x86 -O3)
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 56.320923501171_wp, 0.0_wp, &
        -56.320923501171_wp, 179.664747671772880215_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(s12, 19993558.287_wp, 0.5e-3_wp)

    tstg9 = r

    end function tstg9

    integer function tstg10()
    ! Check fix for adjust tol1_ bug found 2011-06-25 (Visual Studio 10 rel
    ! + debug)
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 52.784459512564_wp, 0.0_wp, &
        -52.784459512563990912_wp, 179.634407464943777557_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(s12, 19991596.095_wp, 0.5e-3_wp)

    tstg10 = r

    end function tstg10

    integer function tstg11()
    ! Check fix for bet2 = -bet1 bug found 2011-06-25 (Visual Studio 10 rel
    ! + debug)
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 48.522876735459_wp, 0.0_wp, &
        -48.52287673545898293_wp, 179.599720456223079643_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(s12, 19989144.774_wp, 0.5e-3_wp)

    tstg11 = r

    end function tstg11

    integer function tstg12()
    ! Check fix for inverse geodesics on extreme prolate/oblate ellipsoids
    ! Reported 2012-08-29 Stefan Guenther <stefan.gunther@embl.de>; fixed
    ! 2012-10-07
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    a = 89.8_wp
    f = -1.83_wp
    outmask = 0
    r = 0
    call inverse(a, f,0.0_wp, 0.0_wp, -10.0_wp, 160.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 120.27_wp, 1.0e-2_wp)
    r = r + assert(azi2, 105.15_wp, 1.0e-2_wp)
    r = r + assert(s12,  266.7_wp,  1.0e-1_wp)

    tstg12 = r

    end function tstg12

    integer function tstg14()
    ! Check fix for inverse ignoring lon12 = nan
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f,0.0_wp, 0.0_wp, 1.0_wp, LatFix(91.0_wp), &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(azi1)
    r = r + chknan(azi2)
    r = r + chknan(s12)

    tstg14 = r

    end function tstg14

    integer function tstg15()
    ! Initial implementation of Math::eatanhe was wrong for e^2 < 0.  This
    ! checks that this is fixed.
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask, flags

    a = 6.4e6_wp
    f = -1.0_wp/150.0_wp
    outmask = 8
    flags = 0
    r = 0
    call direct(a, f, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(SS12, 23700.0_wp, 0.5_wp)

    tstg15 = r

    end function tstg15

    integer function tstg17()
    ! Check fix for LONG_UNROLL bug found on 2015-05-07
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    flags = 2
    r = 0
    call direct(a, f, 40.0_wp, -75.0_wp, -10.0_wp, 2.0e7_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(lat2, -39.0_wp, 1.0_wp)
    r = r + assert(lon2, -254.0_wp, 1.0_wp)
    r = r + assert(azi2, -170.0_wp, 1.0_wp)
    flags = 0
    call direct(a, f, 40.0_wp, -75.0_wp, -10.0_wp, 2.0e7_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(lat2, -39.0_wp, 1.0_wp)
    r = r + assert(lon2, 105.0_wp, 1.0_wp)
    r = r + assert(azi2, -170.0_wp, 1.0_wp)

    tstg17 = r

    end function tstg17

    integer function tstg26()
    ! Check 0/0 problem with area calculation on sphere 2015-09-08
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    a = 6.4d6
    f = 0
    outmask = 8
    r = 0
    call inverse(a, f, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(SS12, 49911046115.0_wp, 0.5_wp)

    tstg26 = r

    end function tstg26

    integer function tstg28()
    ! Check fix for LONG_UNROLL bug found on 2015-05-07
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask, flags

    a = 6.4e6_wp
    f = 0.1e0_wp
    outmask = 1 + 2 + 4 + 8
    flags = 0
    r = 0
    call direct(a, f, 1.0_wp, 2.0_wp, 10.0_wp, 5.0e6_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(a12, 48.55570690_wp, 0.5e-8_wp)

    tstg28 = r

    end function tstg28

    integer function tstg33()
    ! Check max(-0.0,+0.0) issues 2015-08-22 (triggered by bugs in Octave --
    ! sind(-0.0) = +0.0 -- and in some version of Visual Studio --
    ! fmod(-0.0, 360.0) = +0.0.
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 179.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(azi2, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 19926189.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 179.5_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 55.96650_wp, 0.5e-5_wp)
    r = r + assert(azi2, 124.03350_wp, 0.5e-5_wp)
    r = r + assert(s12, 19980862.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 0.00000_wp, 0.5e-5_wp)
    r = r + assert(abs(azi2), 180.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 20003931.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 1.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 0.00000_wp, 0.5e-5_wp)
    r = r + assert(abs(azi2), 180.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 19893357.0_wp, 0.5_wp)
    a = 6.4d6
    f = 0
    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 179.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(azi2, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 19994492.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 0.00000_wp, 0.5e-5_wp)
    r = r + assert(abs(azi2), 180.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 20106193.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 1.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 0.00000_wp, 0.5e-5_wp)
    r = r + assert(abs(azi2), 180.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 19994492.0_wp, 0.5_wp)
    a = 6.4d6
    f = -1/300.0_wp
    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 179.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(azi2, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 19994492.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(azi2, 90.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 20106193.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 0.5_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 33.02493_wp, 0.5e-5_wp)
    r = r + assert(azi2, 146.97364_wp, 0.5e-5_wp)
    r = r + assert(s12, 20082617.0_wp, 0.5_wp)
    call inverse(a, f,0.0_wp, 0.0_wp, 1.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 0.00000_wp, 0.5e-5_wp)
    r = r + assert(abs(azi2), 180.00000_wp, 0.5e-5_wp)
    r = r + assert(s12, 20027270.0_wp, 0.5_wp)

    tstg33 = r

    end function tstg33

    integer function tstg55()
    ! Check fix for nan + point on equator or pole not returning all nans in
    ! Geodesic::Inverse, found 2015-09-23.
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 91.0_wp, 0.0_wp, 0.0_wp, 90.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(azi1)
    r = r + chknan(azi2)
    r = r + chknan(s12)
    call inverse(a, f, 91.0_wp, 0.0_wp, 90.0_wp, 9.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(azi1)
    r = r + chknan(azi2)
    r = r + chknan(s12)

    tstg55 = r

    end function tstg55

    integer function tstg59()
    ! Check for points close with longitudes close to 180 deg apart.
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 5.0_wp, 0.00000000000001_wp, 10.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 0.000000000000035_wp,  1.5e-14_wp)
    r = r + assert(azi2, 179.99999999999996_wp, 1.5e-14_wp)
    r = r + assert(s12,  18345191.174332713_wp, 5.0e-9_wp)

    tstg59 = r

    end function tstg59

    integer function tstg61()
    ! Make sure small negative azimuths are west-going
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    flags = 2
    r = 0
    call direct(a, f, 45.0_wp, 0.0_wp, -0.000000000000000003_wp, 1.0e7_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(lat2, 45.30632_wp, 0.5e-5_wp)
    r = r + assert(lon2, -180.0_wp, 0.5e-5_wp)
    r = r + assert(abs(azi2), 180.0_wp, 0.5e-5_wp)

    tstg61 = r

    end function tstg61

    integer function tstg73()
    ! Check for backwards from the pole bug reported by Anon on 2016-02-13.
    ! This only affected the Java implementation.  It was introduced in Java
    ! version 1.44 and fixed in 1.46-SNAPSHOT on 2016-01-17.
    ! Also the + sign on azi2 is a check on the normalizing of azimuths
    ! (converting -0.0 to +0.0).
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    flags = 0
    r = 0
    call direct(a, f, 90.0_wp, 10.0_wp, 180.0_wp, -1.0e6_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(lat2, 81.04623_wp, 0.5e-5_wp)
    r = r + assert(lon2, -170.0_wp,   0.5e-5_wp)
    r = r + assert(azi2, 0.0_wp,      0.0_wp)
    r = r + assert(sign(1.0_wp, azi2), 1.0_wp, 0.0_wp)

    tstg73 = r

    end function tstg73

    integer function tstg74()
    ! Check fix for inaccurate areas, bug introduced in v1.46, fixed
    ! 2015-10-16.
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 1 + 2 + 4 + 8
    r = 0
    call inverse(a, f, 54.1589_wp, 15.3872_wp, 54.1591_wp, 15.3877_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 55.723110355_wp, 5e-9_wp)
    r = r + assert(azi2, 55.723515675_wp, 5e-9_wp)
    r = r + assert(s12,  39.527686385_wp, 5e-9_wp)
    r = r + assert(a12,   0.000355495_wp, 5e-9_wp)
    r = r + assert(m12,  39.527686385_wp, 5e-9_wp)
    r = r + assert(MM12,  0.999999995_wp, 5e-9_wp)
    r = r + assert(MM21,  0.999999995_wp, 5e-9_wp)
    r = r + assert(SS12, 286698586.30197_wp, 5.0e-4_wp)

    tstg74 = r

    end function tstg74

    integer function tstg76()
    ! The distance from Wellington and Salamanca (a classic failure of
    ! Vincenty
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, &
        -(41+19/60.0_wp), 174+49/60.0_wp, 40+58/60.0_wp, -(5+30/60.0_wp), &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 160.39137649664_wp, 0.5e-11_wp)
    r = r + assert(azi2,  19.50042925176_wp, 0.5e-11_wp)
    r = r + assert(s12,  19960543.857179_wp, 0.5e-6_wp)

    tstg76 = r

    end function tstg76

    integer function tstg78()
    ! An example where the NGS calculator fails to converge
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, 27.2_wp, 0.0_wp, -27.1_wp, 179.5_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1,  45.82468716758_wp, 0.5e-11_wp)
    r = r + assert(azi2, 134.22776532670_wp, 0.5e-11_wp)
    r = r + assert(s12,  19974354.765767_wp, 0.5e-6_wp)

    tstg78 = r

    end function tstg78

    integer function tstg80()
    ! Some tests to add code coverage: computing scale in special cases +
    ! zero length geodesic (includes GeodSolve80 - GeodSolve83).
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 4
    r = 0

    call inverse(a, f,0.0_wp, 0.0_wp, 0.0_wp, 90.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(MM12, -0.00528427534_wp, 0.5e-10_wp)
    r = r + assert(MM21, -0.00528427534_wp, 0.5e-10_wp)

    call inverse(a, f,0.0_wp, 0.0_wp, 1e-6_wp, 1e-6_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(MM12, 1.0_wp, 0.5e-10_wp)
    r = r + assert(MM21, 1.0_wp, 0.5e-10_wp)

    outmask = 15
    call inverse(a, f, 20.001_wp, 0.0_wp, 20.001_wp, 0.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(a12,  0.0_wp,   1.0e-13_wp)
    r = r + assert(s12,  0.0_wp,   1.0e-8_wp )
    r = r + assert(azi1, 180.0_wp, 1.0e-13_wp)
    r = r + assert(azi2, 180.0_wp, 1.0e-13_wp)
    r = r + assert(m12,  0.0_wp,   1.0e-8_wp )
    r = r + assert(MM12, 1.0_wp,   1.0e-15_wp)
    r = r + assert(MM21, 1.0_wp,   1.0e-15_wp)
    r = r + assert(SS12, 0.0_wp,   1.0e-10_wp)
    r = r + assert(sign(1.0_wp, a12), 1.0_wp, 0.0_wp)
    r = r + assert(sign(1.0_wp, s12), 1.0_wp, 0.0_wp)
    r = r + assert(sign(1.0_wp, m12), 1.0_wp, 0.0_wp)

    call inverse(a, f, 90.0_wp, 0.0_wp, 90.0_wp, 180.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(a12,  0.0_wp,   1.0e-13_wp)
    r = r + assert(s12,  0.0_wp,   1.0e-8_wp )
    r = r + assert(azi1, 0.0_wp,   1.0e-13_wp)
    r = r + assert(azi2, 180.0_wp, 1.0e-13_wp)
    r = r + assert(m12,  0.0_wp,   1.0e-8_wp )
    r = r + assert(MM12, 1.0_wp,   1.0e-15_wp)
    r = r + assert(MM21, 1.0_wp,   1.0e-15_wp)
    r = r + assert(SS12, 127516405431022.0_wp, 0.5_wp)

    tstg80 = r

    end function tstg80

    integer function tstg84()
    ! Tests for python implementation to check fix for range errors with
    ! {fmod,sin,cos}(inf) (includes GeodSolve84 - GeodSolve86).
    real(wp) lat2, lon2, azi2, a12, m12, MM12, MM21, SS12
    real(wp) a, f, nan, inf
    integer r, outmask, flags

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    flags = 0
    inf = 1.0_wp/LatFix(0.0_wp)
    nan = LatFix(91.0_wp)
    r = 0
    call direct(a, f,0.0_wp, 0.0_wp, 90.0_wp, inf, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(lat2)
    r = r + chknan(lon2)
    r = r + chknan(azi2)
    call direct(a, f,0.0_wp, 0.0_wp, 90.0_wp, nan, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(lat2)
    r = r + chknan(lon2)
    r = r + chknan(azi2)
    call direct(a, f,0.0_wp, 0.0_wp, inf, 1000.0_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(lat2)
    r = r + chknan(lon2)
    r = r + chknan(azi2)
    call direct(a, f,0.0_wp, 0.0_wp, nan, 1000.0_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(lat2)
    r = r + chknan(lon2)
    r = r + chknan(azi2)
    call direct(a, f, 0.0_wp, inf, 90.0_wp, 1000.0_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(lat2,0.0_wp, 0.0_wp)
    r = r + chknan(lon2)
    r = r + assert(azi2, 90.0_wp, 0.0_wp)
    call direct(a, f, 0.0_wp, nan, 90.0_wp, 1000.0_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(lat2,0.0_wp, 0.0_wp)
    r = r + chknan(lon2)
    r = r + assert(azi2, 90.0_wp, 0.0_wp)
    call direct(a, f, inf, 0.0_wp, 90.0_wp, 1000.0_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(lat2)
    r = r + chknan(lon2)
    r = r + chknan(azi2)
    call direct(a, f, nan, 0.0_wp, 90.0_wp, 1000.0_wp, &
        flags, lat2, lon2, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(lat2)
    r = r + chknan(lon2)
    r = r + chknan(azi2)

    tstg84 = r

    end function tstg84

    integer function tstg92()
    ! Check fix for inaccurate hypot with python 3.[89].  Problem reported
    ! by agdhruv https://github.com/geopy/geopy/issues/466 ; see
    ! https://bugs.python.org/issue43088
    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f, &
        37.757540000000006_wp, -122.47018_wp, &
        37.75754_wp,           -122.470177_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(azi1, 89.99999923_wp, 1.0e-7_wp)
    r = r + assert(azi2, 90.00000106_wp, 1.0e-7_wp)
    r = r + assert(s12,   0.264_wp,      0.5e-3_wp)

    tstg92 = r

    end function tstg92

    integer function tstg94()
    ! Check fix for lat2 = nan being treated as lat2 = 0 (bug found
    ! 2021-07-26)

    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    outmask = 0
    r = 0
    call inverse(a, f,0.0_wp, 0.0_wp, LatFix(91.0_wp), 90.0_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + chknan(azi1)
    r = r + chknan(azi2)
    r = r + chknan(s12)

    tstg94 = r

    end function tstg94

    integer function tstg96()
    ! Failure with long doubles found with test case from Nowak + Nowak Da
    ! Costa (2022).  Problem was using somg12 > 1 as a test that it needed
    ! to be set when roundoff could result in somg12 slightly bigger that 1.
    ! Found + fixed 2022-03-30.

    real(wp) azi1, azi2, s12, a12, m12, MM12, MM21, SS12
    real(wp) a, f
    integer r, outmask

    a = 6378137.0_wp
    f = 1/298.257222101_wp
    outmask = 8
    r = 0
    call inverse(a,f, 0.0_wp, 0.0_wp, 60.0832522871723_wp, 89.8492185074635_wp, &
        s12, azi1, azi2, outmask, a12, m12, MM12, MM21, SS12)
    r = r + assert(SS12, 42426932221845.0_wp, 0.5_wp)

    tstg96 = r

    end function tstg96

    integer function tstp0()
    ! Check fix for pole-encircling bug found 2011-03-16
    real(wp) lata(4), lona(4)
    data lata / 89.0_wp, 89.0_wp, 89.0_wp, 89.0_wp /
    data lona / 0.0_wp, 90.0_wp, 180.0_wp, 270.0_wp /
    real(wp) latb(4), lonb(4)
    data latb / -89.0_wp, -89.0_wp, -89.0_wp, -89.0_wp /
    data lonb / 0.0_wp, 90.0_wp, 180.0_wp, 270.0_wp /
    real(wp) latc(4), lonc(4)
    data latc / 0.0_wp, -1.0_wp, 0.0_wp, 1.0_wp /
    data lonc / -1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp /
    real(wp) latd(3), lond(3)
    data latd / 90.0_wp,0.0_wp, 0.0_wp /
    data lond /0.0_wp, 0.0_wp, 90.0_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lata, lona, 4, AA, PP)
    r = r + assert(PP, 631819.8745_wp, 1.0e-4_wp)
    r = r + assert(AA, 24952305678.0_wp, 1.0_wp)

    call area(a, f, latb, lonb, 4, AA, PP)
    r = r + assert(PP, 631819.8745_wp, 1.0e-4_wp)
    r = r + assert(AA, -24952305678.0_wp, 1.0_wp)

    call area(a, f, latc, lonc, 4, AA, PP)
    r = r + assert(PP, 627598.2731_wp, 1.0e-4_wp)
    r = r + assert(AA, 24619419146.0_wp, 1.0_wp)

    call area(a, f, latd, lond, 3, AA, PP)
    r = r + assert(PP, 30022685.0_wp, 1.0_wp)
    r = r + assert(AA, 63758202715511.0_wp, 1.0_wp)

    tstp0 = r

    end function tstp0

    integer function tstp5()
    ! Check fix for Planimeter pole crossing bug found 2011-06-24
    real(wp) lat(3), lon(3)
    data lat / 89.0_wp, 89.0_wp, 89.0_wp /
    data lon / 0.1_wp, 90.1_wp, -179.9_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lat, lon, 3, AA, PP)
    r = r + assert(PP, 539297.0_wp, 1.0_wp)
    r = r + assert(AA, 12476152838.5_wp, 1.0_wp)

    tstp5 = r
    return
    end

    integer function tstp6()
    ! Check fix for pole-encircling bug found 2011-03-16
    real(wp) lata(3), lona(3)
    data lata / 9.0_wp, 9.0_wp, 9.0_wp /
    data lona / -0.00000000000001_wp, 180.0_wp, 0.0_wp /
    real(wp) latb(3), lonb(3)
    data latb / 9.0_wp, 9.0_wp, 9.0_wp /
    data lonb / 0.00000000000001_wp, 0.0_wp, 180.0_wp /
    real(wp) latc(3), lonc(3)
    data latc / 9.0_wp, 9.0_wp, 9.0_wp /
    data lonc / 0.00000000000001_wp, 180.0_wp, 0.0_wp /
    real(wp) latd(3), lond(3)
    data latd / 9.0_wp, 9.0_wp, 9.0_wp /
    data lond / -0.00000000000001_wp, 0.0_wp, 180.0_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lata, lona, 3, AA, PP)
    r = r + assert(PP, 36026861.0_wp, 1.0_wp)
    r = r + assert(AA, 0.0_wp, 1.0_wp)

    tstp6 = r
    return
    end

    integer function tstp12()
    ! Area of arctic circle (not really -- adjunct to rhumb-area test)
    real(wp) lat(3), lon(3)
    data lat / 66.562222222_wp, 66.562222222_wp, 66.562222222_wp /
    data lon / 0.0_wp, 180.0_wp, 360.0_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lat, lon, 3, AA, PP)
    r = r + assert(PP, 10465729.0_wp, 1.0_wp)
    r = r + assert(AA, 0.0_wp, 1.0_wp)

    tstp12 = r
    return
    end

    integer function tstp12r()
    ! reverse area of arctic circle
    real(wp) lat(3), lon(3)
    data lat / 66.562222222_wp, 66.562222222_wp, 66.562222222_wp /
    data lon / -0.0_wp, -180.0_wp, -360.0_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lat, lon, 3, AA, PP)
    r = r + assert(PP, 10465729.0_wp, 1.0_wp)
    r = r + assert(AA, 0.0_wp, 1.0_wp)

    tstp12r = r
    return
    end

    integer function tstp13()
    ! Check encircling pole twice
    real(wp) lat(6), lon(6)
    data lat / 89.0_wp, 89.0_wp, 89.0_wp, 89.0_wp, 89.0_wp, 89.0_wp /
    data lon / -360.0_wp, -240.0_wp, -120.0_wp, 0.0_wp, 120.0_wp, 240.0_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lat, lon, 6, AA, PP)
    r = r + assert(PP, 1160741.0_wp, 1.0_wp)
    r = r + assert(AA, 32415230256.0_wp, 1.0_wp)

    tstp13 = r
    return
    end

    integer function tstp15()
    ! Coverage tests, includes Planimeter15 - Planimeter18 (combinations of
    ! reverse and sign).  But flags aren't supported in the Fortran
    ! implementation.
    real(wp) lat(3), lon(3)
    data lat / 2.0_wp, 1.0_wp, 3.0_wp /
    data lon / 1.0_wp, 2.0_wp, 3.0_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lat, lon, 3, AA, PP)
    r = r + assert(AA, 18454562325.45119_wp, 1.0_wp)
    ! Interchanging lat and lon is equivalent to traversing the polygon
    ! backwards.
    call area(a, f, lon, lat, 3, AA, PP)
    r = r + assert(AA, -18454562325.45119_wp, 1.0_wp)

    tstp15 = r
    return
    end

    integer function tstp19()
    ! Coverage tests, includes Planimeter19 - Planimeter20 (degenerate
    ! polygons).
    real(wp) lat(1), lon(1)
    data lat / 1.0_wp /
    data lon / 1.0_wp /
    real(wp) a, f, AA, PP
    integer r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0

    call area(a, f, lat, lon, 1, AA, PP)
    r = r + assert(AA,0.0_wp, 0.0_wp)
    r = r + assert(PP,0.0_wp, 0.0_wp)

    tstp19 = r
    return
    end

    integer function tstp21()
    ! Some test to add code coverage: multiple circlings of pole (includes
    ! Planimeter21 - Planimeter28).
    real(wp),parameter :: lat(12) = 45.0_wp
    real(wp),parameter :: lon(12) = [60.0_wp, 180.0_wp, -60.0_wp, &
                                    60.0_wp, 180.0_wp, -60.0_wp, &
                                    60.0_wp, 180.0_wp, -60.0_wp, &
                                    60.0_wp, 180.0_wp, -60.0_wp]
    real(wp),parameter :: lonr(12) = [-60.0_wp, 180.0_wp, 60.0_wp, &
                                    -60.0_wp, 180.0_wp, 60.0_wp, &
                                    -60.0_wp, 180.0_wp, 60.0_wp, &
                                    -60.0_wp, 180.0_wp, 60.0_wp]

    real(wp) :: a, f, AA, PP, AA1
    integer :: i, r

    ! WGS84 values
    a = 6378137.0_wp
    f = 1.0_wp/298.257223563_wp
    r = 0
    ! Area for one circuit
    AA1 = 39433884866571.4277_wp

    do i = 3,4
    call area(a, f, lat, lon, 3*i, AA, PP)
    r = r + assert(AA, AA1*i, 0.5_wp)
    call area(a, f, lat, lonr, 3*i, AA, PP)
    r = r + assert(AA, -AA1*i, 0.5_wp)
    end do

    tstp21 = r

    end function tstp21

end program