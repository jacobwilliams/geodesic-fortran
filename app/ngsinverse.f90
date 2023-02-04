!b::inverse
!
program inverse
  use geodesic_module
!
!********1*********2*********3*********4*********5*********6*********7**
!
! name:      inverse
! version:   201211.29
! author:    stephen j. frakes
! last mod:  Charles Karney
! purpose:   to compute a geodetic inverse
!            and then display output information
!
! input parameters:
! -----------------
!
! output parameters:
! ------------------
!
! local variables and constants:
! ------------------------------
! answer           user prompt response
! b                semiminor axis polar (in meters)
! baz              azimuth back (in radians)
! buff             input char buffer array
! dd,dm,ds         temporary values for degrees, minutes, seconds
! dlon             temporary value for difference in longitude (radians)
! dmt,d_mt         char constants for meter units
! edist            ellipsoid distance (in meters)
! elips            ellipsoid choice
! esq              eccentricity squared for reference ellipsoid
! faz              azimuth forward (in radians)
! filout           output file name
! finv             reciprocal flattening
! hem              hemisphere flag for lat & lon entry
! ierror           error condition flag with d,m,s conversion
! lgh              length of buff() array
! option           user prompt response
! r1,r2            temporary variables
! ss               temporary variable
! tol              tolerance for conversion of seconds
!
! name1            name of station one
! ld1,lm1,sl1      latitude  sta one - degrees,minutes,seconds
! ald1,alm1,sl1    latitude  sta one - degrees,minutes,seconds
! lat1sn           latitude  sta one - sign (+/- 1)
! d_ns1            latitude  sta one - char ('N','S')
! lod1,lom1,slo1   longitude sta one - degrees,minutes,seconds
! alod1,alom1,slo1 longitude sta one - degrees,minutes,seconds
! lon1sn           longitude sta one - sign (+/- 1)
! d_ew1            longitude sta one - char ('E','W')
! iaz1,maz1,saz1   forward azimuth   - degrees,minutes,seconds
! isign1           forward azimuth   - flag  (+/- 1)
! glat1,glon1      station one       - (lat & lon in radians )
! p1,e1            standpoint one    - (lat & lon in radians )
!
! name2            name of station two
! ld2,lm2,sl2      latitude  sta two - degrees,minutes,seconds
! ald2,alm2,sl2    latitude  sta two - degrees,minutes,seconds
! lat2sn           latitude  sta two - sign (+/- 1)
! d_ns2            latitude  sta two - char ('N','S')
! lod2,lom2,slo2   longitude sta two - degrees,minutes,seconds
! alod2,alom2,slo2 longitude sta one - degrees,minutes,seconds
! lon2sn           longitude sta two - sign (+/- 1)
! d_ew2            longitude sta two - char ('E','W')
! iaz2,maz2,saz2   back azimuth      - degrees,minutes,seconds
! isign2           back azimuth      - flag  (+/- 1)
! glat2,glon2      station two       - (lat & lon in radians )
! p2,e2            forepoint two     - (lat & lon in radians )
!
! global variables and constants:
! -------------------------------
! a                semimajor axis equatorial (in meters)
! f                flattening
! pi               constant 3.14159....
! rad              constant 180.0/pi
!
!    this module called by:  n/a
!
!    this module calls:      elipss, getdeg, inver1, todmsp
!    gethem, trim,   bufdms, gvalr8, gvali4, fixdms, gpnhri ***********
!    gethem, trim,   bufdms, gvalr8, gvali4, fixdms, invers <----------
!    datan,  write,  read,   dabs,   open,   stop
!
!    include files used:     n/a
!
!    common blocks used:     const, elipsoid
!
!    references:             see comments within subroutines
!
!    comments:
!
!********1*********2*********3*********4*********5*********6*********7**
!::modification history
!::1990mm.dd, sjf, creation of program
!::199412.15, bmt, creation of program on viper
!::200203.08, crs, modified by c.schwarz to correct spelling of Clarke
!::200207.15, rws, modified i/o & standardized program documentation
!::                added subs trim, bufdms, gethem, gvali4, gvalr8
!::200207.23, rws, replaced sub inver1 with gpnarc, gpnloa, gpnhri
!::200208.15, rws, fixed an error in bufdms
!::              - renamed ellips to elipss "common error" with dirct2
!::              - added FAZ & BAZ to printed output
!::200208.19, rws, added more error flags for web interface code
!::              - added logical nowebb
!::200208.xx, sjf, program version number 2.0
!::201105.xx, dgm, program version number 3.0
!::              - replaced sub gpnarc, gpnloa, gpnhri with invers
!::              - tested for valid antipodal solutions (+/- 0.1 mm)
!::              - tested for polar solutions (+/- 0.1 mm)
!::              - needs improvement for long-line/antipodal boundary
!::201211.29, cffk, program version numer 3.1
!::              - drop in replacement routines from
!::                "Algorithms for Geodesics"
!********1*********2*********3*********4*********5*********6*********7**
!e::inverse
!
implicit double precision (a-h, o-z)
implicit integer (i-n)
!
logical  nowebb
!
character*1  answer,option,dmt,buff(50),hem
character*6  d_ns1, d_ew1, d_ns2, d_ew2, d_mt
character*30 filout,name1,name2,elips
!
integer*4    ierror
integer*4    lgh
!
common/const/pi,rad
common/elipsoid/a,f
!
!     ms_unix      0 = web version
!                  1 = ms_dos or unix version
!
parameter   ( ms_unix = 0 )
!
pi   = 4.d0*datan(1.d0)
rad  = 180.d0/pi
dmt  = 'm'
d_mt = 'Meters'
!
if( ms_unix==1 )then
  nowebb = .true.
else
  nowebb = .false.
endif
!
1 do 2 i=1,25
  write(*,*) '  '
2 continue
!
5 write(*,*) '  Program Inverse  -  Version 3.1 '
write(*,*) '  '
write(*,*) '  Ellipsoid options : '
write(*,*) '  '
write(*,*) '  1) GRS80 / WGS84  (NAD83) '
write(*,*) '  2) Clarke 1866    (NAD27) '
write(*,*) '  3) Any other ellipsoid '
write(*,*) '  '
write(*,*) '  Enter choice : '
read(*,10) option
10 format(a1)
!
if(option=='1') then
  a=6378137.d0
  f=1.d0/298.257222100882711243162836600094d0
  elips='GRS80 / WGS84  (NAD83)'
elseif(option=='2') then
  a=6378206.4d0
  f=1.d0/294.9786982138d0
  elips='Clarke 1866    (NAD27)'
elseif(option=='3') then
  call elipss (elips)
else
  write(*,*) '  Enter 1, 2, or 3 !   Try again --'
  goto 5
endif
!
esq = f*(2.0d0-f)
!
write(*,*) '  '
write(*,*) '  '
write(*,*) '  '
write(*,*) '  '
!
15 write(*,*) '  Enter First Station '
write(*,16)
16 format(18x,'(Separate D,M,S by blanks or commas)')
write(*,*) 'hDD MM SS.sssss  Latitude :        (h default = N )'
11 format(50a1)
!
22 hem='N'
read(*,11) buff
call trim (buff,lgh,hem)
call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
!
if( ierror==0 )then
  irlat1 = 0
else
  irlat1 = 1
  write(*,*) ' Invalid Latitude ... Please re-enter it '
  write(*,*) '  '
  if( nowebb )then
    goto 22
  endif
endif
!
ald1 = dd
alm1 = dm
sl1  = ds
!
if( hem=='N' ) then
  lat1sn = +1
else
  lat1sn = -1
endif
!
write(*,*) 'hDDD MM SS.sssss Longitude :       (h default = W )'
!
23 hem='W'
read(*,11) buff
call trim (buff,lgh,hem)
call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
!
if( ierror==0 )then
  irlon1 = 0
else
  irlon1 = 1
  write(*,*) ' Invalid Longitude ... Please re-enter it '
  write(*,*) '  '
  if( nowebb )then
    goto 23
  endif
endif
!
alod1 = dd
alom1 = dm
slo1  = ds
!
if( hem=='E' ) then
  lon1sn = +1
else
  lon1sn = -1
endif
!
call getdeg(ald1, alm1, sl1, lat1sn, glat1)
call getdeg(alod1,alom1,slo1,lon1sn, glon1)
!
write(*,*) '  '
write(*,*) '  '
!
20 write(*,*) '  Enter Second Station '
write(*,16)
write(*,*) 'hDD MM SS.sssss  Latitude :        (h default = N )'
!
24 hem='N'
read(*,11) buff
call trim (buff,lgh,hem)
call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
!
if( ierror==0 )then
  irlat2 = 0
else
  irlat2 = 1
  write(*,*) ' Invalid Latitude ... Please re-enter it '
  write(*,*) '  '
  if( nowebb )then
    goto 24
  endif
endif
!
ald2 = dd
alm2 = dm
sl2  = ds
!
if( hem=='N' ) then
  lat2sn = +1
else
  lat2sn = -1
endif
!
write(*,*) 'hDDD MM SS.sssss Longitude :       (h default = W )'
!
25 hem='W'
read(*,11) buff
call trim (buff,lgh,hem)
call bufdms (buff,lgh,hem,dd,dm,ds,ierror)
!
if( ierror==0 )then
  irlon2 = 0
else
  irlon2 = 1
  write(*,*) ' Invalid Longitude ... Please re-enter it '
  write(*,*) '  '
  if( nowebb )then
    goto 25
  endif
endif
!
alod2 = dd
alom2 = dm
slo2  = ds
!
if( hem=='E' )then
  lon2sn = +1
else
  lon2sn = -1
endif
!
call getdeg(ald2, alm2, sl2, lat2sn, glat2)
call getdeg(alod2,alom2,slo2,lon2sn, glon2)
!
p1 = glat1
e1 = glon1
p2 = glat2
e2 = glon2
!
if( e1<0.0d0 )then
  e1 = e1+2.0d0*pi
endif
if( e2<0.0d0 )then
  e2 = e2+2.0d0*pi
endif
!
!     compute the geodetic inverse
!
call invers(a, f, p1, e1, p2, e2, &
    edist, faz, baz, 0, dummy, dummy, dummy, dummy, dummy)
if (baz >= 0) then
  baz = baz - 180
else
  baz = baz + 180
end if
!
!     set the tolerance (in seconds) for the azimuth conversion
!
tol = 0.00005d0
!
call todmsp(faz,iaz1,maz1,saz1,isign1)
if(isign1<0) then
  iaz1=359-iaz1
  maz1=59-maz1
  saz1=60.d0-saz1
endif
call fixdms ( iaz1, maz1, saz1, tol )
!
call todmsp(baz,iaz2,maz2,saz2,isign2)
if(isign2<0) then
  iaz2=359-iaz2
  maz2=59-maz2
  saz2=60.d0-saz2
endif
call fixdms ( iaz2, maz2, saz2, tol )
!
call todmsp(glat1, ld1,  lm1,  sl1,  lat1sn)
call todmsp(glon1, lod1, lom1, slo1, lon1sn)
call todmsp(glat2, ld2,  lm2,  sl2,  lat2sn)
call todmsp(glon2, lod2, lom2, slo2, lon2sn)
!
call hem_ns ( lat1sn, d_ns1 )
call hem_ew ( lon1sn, d_ew1 )
call hem_ns ( lat2sn, d_ns2 )
call hem_ew ( lon2sn, d_ew2 )
!
write(*,*) '  '
write(*,*) '  '
write(*,*) '  '
write(*,*) '  '
write(*,*) '  '
write(*,49) elips
49 format('  Ellipsoid : ',a30)
finv=1.d0/f
b=a*(1.d0-f)
write(*,50) a,b,finv
50 format('  Equatorial axis,    a   = ',f15.4,/, &
       '  Polar axis,         b   = ',f15.4,/, &
       '  Inverse flattening, 1/f = ',f16.11)
!
18 format('    LAT = ',i3,1x,i2,1x,f8.5,1x,a6)
19 format('    LON = ',i3,1x,i2,1x,f8.5,1x,a6)
!
write(*,*) '  '
write(*,*) '  First  Station : '
write(*,*) '  ---------------- '
write(*,18) ld1, lm1, sl1, d_ns1
write(*,19) lod1,lom1,slo1,d_ew1
!
write(*,*) '  '
write(*,*) '  Second Station : '
write(*,*) '  ---------------- '
write(*,18) ld2, lm2, sl2, d_ns2
write(*,19) lod2,lom2,slo2,d_ew2
!
32 format('  Ellipsoidal distance     S = ',f14.4,1x,a1)
34 format('  Forward azimuth        FAZ = ',i3,1x,i2,1x,f7.4, &
       ' From North')
35 format('  Back azimuth           BAZ = ',i3,1x,i2,1x,f7.4, &
       ' From North')
!
write(*,*) '  '
write(*,34) iaz1,maz1,saz1
write(*,35) iaz2,maz2,saz2
write(*,32) edist,dmt
write(*,*) '  '
write(*,*) '  Do you want to save this output into a file (y/n)?'
read(*,10) answer
!
if( answer=='Y'.or.answer=='y' )then
39   write(*,*) '  Enter the output filename : '
  read(*,40) filout
40   format(a30)
  open(3,file=filout,status='new',err=99)
  goto 98
!
99   write(*,*) '  File already exists, try again.'
  go to 39
!
98   continue
  write(3,*) '  '
  write(3,49) elips
  finv=1.d0/f
  b=a*(1.d0-f)
  write(3,50) a,b,finv
  write(*,*) '  Enter the First Station name : '
  read(*,40) name1
  write(*,*) '  Enter the Second Station name : '
  read(*,40) name2
!
41   format('  First  Station : ',a30)
42   format('  Second Station : ',a30)
84   format('  Error: First  Station ... Invalid Latitude  ')
85   format('  Error: First  Station ... Invalid Longitude ')
86   format('  Error: Second Station ... Invalid Latitude  ')
87   format('  Error: Second Station ... Invalid Longitude ')
88   format(1x,65(1h*))
91   format('         DD(0-89) MM(0-59) SS(0-59.999...)  ')
92   format('         DDD(0-359) MM(0-59) SS(0-59.999...)  ')
!
  write(3,*) '  '
  write(3,41) name1
  write(3,*) '  ---------------- '

  if( irlat1==0 )then
    write(3,18) ld1, lm1, sl1, d_ns1
  else
    write(3,88)
    write(3,84)
    write(3,91)
    write(3,88)
  endif
!
  if( irlon1==0 )then
    write(3,19) lod1,lom1,slo1,d_ew1
  else
    write(3,88)
    write(3,85)
    write(3,92)
    write(3,88)
  endif
!
  write(3,*) '  '
  write(3,42) name2
  write(3,*) '  ---------------- '
!
  if( irlat2==0 )then
    write(3,18) ld2, lm2, sl2, d_ns2
  else
    write(3,88)
    write(3,86)
    write(3,91)
    write(3,88)
  endif
!
  if( irlon2==0 )then
    write(3,19) lod2,lom2,slo2,d_ew2
  else
    write(3,88)
    write(3,87)
    write(3,92)
    write(3,88)
  endif
!
  write(3,*) '  '
  write(3,34) iaz1,maz1,saz1
  write(3,35) iaz2,maz2,saz2
  write(3,32) edist,dmt
  write(3,*) '  '
endif
!
80 write(*,*) '  '
write(*,*) '  '
write(*,*) '  '
write(*,*) '  1) Another inverse, different ellipsoid.'
write(*,*) '  2) Same ellipsoid, two new stations.'
write(*,*) '  3) Same ellipsoid, same First Station.'
write(*,*) '  4) Done.'
write(*,*) '  '
write(*,*) '  Enter choice : '
read(*,10) answer
!
if(     answer=='1' )then
  goto 1
elseif( answer=='2' )then
  goto 15
elseif( answer=='3' )then
  goto 20
else
  write(*,*) '  Thats all, folks!'
endif

!     stop
end
