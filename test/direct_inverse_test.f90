
!*****************************************************************************************
!>
!  Unit test for the [[direct_vincenty]] and [[inverse_vincenty]] geodetic routines.
!
!@note This was originally from the
!      [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit)

    program direct_inverse_test

    use iso_fortran_env, only: wp => real64
    use geodesic_module

    implicit none

    !Ellipsoid : GRS80 / WGS84 (NAD83)
    real(wp),parameter :: a  = 6378137.0000_wp      !! Equatorial radius
    real(wp),parameter :: rf = 298.25722210088_wp   !! Inverse flattening
    real(wp),parameter :: f  = 1.0_wp / rf          !! flattening

    real(wp) :: glat1,glon1,glat2,glon2,faz,baz,s,sig,lam,glat2_,glon2_,baz_
    integer :: it, kind

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' direct_inverse_test'
    write(*,*) '---------------'
    write(*,*) ''

    !specify two points:
    glat1 = 0.523599_wp
    glon1 = 1.74533_wp
    glat2 = 0.698132_wp
    glon2 = 2.0944_wp

    call inverse_vincenty(a,rf,glat1,glon1,glat2,glon2,faz,baz,s,it,sig,lam,kind)
    call direct_vincenty(a,f,glat1,glon1,faz,s,glat2_,glon2_,baz_)

    write(*,*) ' vincenty:'
    write(*,*) 'lat error: ', glat2_ - glat2
    write(*,*) 'lon error: ', glon2_ - glon2
    write(*,*) 'baz error: ', baz_ - baz

    ! call invers(a,rf,glat1,glon1,glat2,glon2,faz,baz,s,it,sig,lam,kind) ! fix these ...
    ! call direct(a,f,glat1,glon1,faz,s,glat2_,glon2_,baz_)

    ! write(*,*) ' geodeticlib:'
    ! write(*,*) 'lat error: ', glat2_ - glat2
    ! write(*,*) 'lon error: ', glon2_ - glon2
    ! write(*,*) 'baz error: ', baz_ - baz

    end program direct_inverse_test
!*****************************************************************************************
