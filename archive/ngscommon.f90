module ngscommon
  use geodesic_module, wp => geodesic_wp
  private :: wp

  contains

subroutine bufdms (buff,lgh,hem,dd,dm,ds,ierror)
implicit real(wp) (a-h, o-z)
implicit integer (i-n)
!
logical     done,flag
!
character*1 buff(*),abuf(21)
character*1 ch
character*1 hem
integer*4   ll,lgh
integer*4   i4,id,im,is,icond,ierror
real(wp)      x(5)
!
!     set the "error flag"
!
ierror = 0
icond  = 0
!
!     set defaults for dd,dm,ds
!
dd = 0.0_wp
dm = 0.0_wp
ds = 0.0_wp
!
!     set default limits for "hem" flag
!
if(     hem=='N' .or. hem=='S' )then
  ddmax = 90.0_wp
elseif( hem=='E' .or. hem=='W' )then
  ddmax = 360.0_wp
elseif( hem=='A' )then
  ddmax = 360.0_wp
elseif( hem=='Z' )then
  ddmax = 180.0_wp
elseif( hem=='*' )then
  ddmax  = 0.0_wp
  ierror = 1
else
  ddmax = 360.0_wp
endif
!
do 1 i=1,5
  x(i) = 0.0_wp
1 continue
!
icolon = 0
ipoint = 0
icount = 0
flag   = .true.
jlgh   = lgh
!
do 2 i=1,jlgh
  if( buff(i)==':' )then
    icolon = icolon+1
  endif
  if( buff(i)=='.' )then
    ipoint = ipoint+1
    flag   = .false.
  endif
  if( flag )then
    icount = icount+1
  endif
2 continue
!
if( ipoint==1 .and. icolon==0 )then
!
!       load temp buffer
!
  do 3 i=1,jlgh
    abuf(i) = buff(i)
3   continue
  abuf(jlgh+1) = '$'
  ll = jlgh
!
  call gvalr8 (abuf,ll,r8,icond)
!
  if( icount>=5 )then
!
!         value is a packed decimal of ==>  DDMMSS.sssss
!
    ss = r8/10000.0_wp
    id = int( ss )
!
    r8 = r8-10000.0_wp*dble(float(id))
    ss = r8/100.0_wp
    im = int( ss )
!
    r8 = r8-100.0_wp*dble(float(im))
  else
!
!         value is a decimal of ==>  .xx   X.xxx   X.
!
    id = int( r8 )
    r8 = (r8-id)*60.0_wp
    im = int( r8 )
    r8 = (r8-im)*60.0_wp
  endif
!
!       account for rounding error
!
  is = nint( r8*1.0d5 )
  if( is>=6000000 )then
     r8 = 0.0_wp
     im = im+1
  endif
!
  if( im>=60 )then
    im = 0
    id = id+1
  endif
!
  dd = dble( float( id ) )
  dm = dble( float( im ) )
  ds = r8
else
!
!       buff() value is a d,m,s of ==>  NN:NN:XX.xxx
!
  k    = 0
  next = 1
  done = .false.
  ie   = jlgh
!
  do 100 j=1,5
    ib = next
    do 90 i=ib,ie
      ch   = buff(i)
      last = i
      if( i==jlgh .or. ch==':' )then
        if( i==jlgh )then
          done = .true.
        endif
        if( ch==':' )then
          last = i-1
        endif
        goto 91
      endif
90     continue
    goto 98
!
91     ipoint = 0
    ik     = 0
    do 92 i=next,last
      ik = ik+1
      ch = buff(i)
      if( ch=='.' )then
        ipoint = ipoint+1
      endif
      abuf(ik) = buff(i)
92     continue
    abuf(ik+1) = '$'
!
    ll = ik
    if( ipoint==0 )then
      call gvali4 (abuf,ll,i4,icond)
      r8 = dble(float( i4 ))
    else
      call gvalr8 (abuf,ll,r8,icond)
    endif
!
    k    = k+1
    x(k) = r8
!
98     if( done )then
      goto 101
    endif
!
    next = last
99     next = next+1
    if( buff(next)==':' )then
      goto 99
    endif
100   continue
!
!       load dd,dm,ds
!
101   if( k>=1 )then
    dd = x(1)
  endif
!
  if( k>=2 )then
    dm = x(2)
  endif
!
  if( k>=3 )then
    ds = x(3)
  endif
endif
!
if( dd>ddmax  .or. &
    dm>=60.0_wp .or. &
    ds>=60.0_wp )then
  ierror = 1
  dd = 0.0_wp
  dm = 0.0_wp
  ds = 0.0_wp
endif
!
if( icond/=0 )then
  ierror = 1
endif
!
return
end

subroutine elipss (elips)
implicit real(wp)(a-h,o-z)
character*1  answer
character*30 elips
common/elipsoid/a,f
write(*,*) '  Other Ellipsoids.'
write(*,*) '  -----------------'
write(*,*) '  '
write(*,*) '  A) Airy 1858'
write(*,*) '  B) Airy Modified'
write(*,*) '  C) Australian National'
write(*,*) '  D) Bessel 1841'
write(*,*) '  E) Clarke 1880'
write(*,*) '  F) Everest 1830'
write(*,*) '  G) Everest Modified'
write(*,*) '  H) Fisher 1960'
write(*,*) '  I) Fisher 1968'
write(*,*) '  J) Hough 1956'
write(*,*) '  K) International (Hayford)'
write(*,*) '  L) Krassovsky 1938'
write(*,*) '  M) NWL-9D (WGS 66)'
write(*,*) '  N) South American 1969'
write(*,*) '  O) Soviet Geod. System 1985'
write(*,*) '  P) WGS 72'
write(*,*) '  Q-Z) User defined.'
write(*,*) '  '
write(*,*) '  Enter choice : '
read(*,10) answer
10 format(a1)
!
if(answer=='A'.or.answer=='a') then
  a=6377563.396d0
  f=1.d0/299.3249646d0
  elips='Airy 1858'
elseif(answer=='B'.or.answer=='b') then
  a=6377340.189d0
  f=1.d0/299.3249646d0
  elips='Airy Modified'
elseif(answer=='C'.or.answer=='c') then
  a=6378160.d0
  f=1.d0/298.25d0
  elips='Australian National'
elseif(answer=='D'.or.answer=='d') then
  a=6377397.155d0
  f=1.d0/299.1528128d0
  elips='Bessel 1841'
elseif(answer=='E'.or.answer=='e') then
  a=6378249.145d0
  f=1.d0/293.465d0
  elips='Clarke 1880'
elseif(answer=='F'.or.answer=='f') then
  a=6377276.345d0
  f=1.d0/300.8017d0
  elips='Everest 1830'
elseif(answer=='G'.or.answer=='g') then
  a=6377304.063d0
  f=1.d0/300.8017d0
  elips='Everest Modified'
elseif(answer=='H'.or.answer=='h') then
  a=6378166.d0
  f=1.d0/298.3d0
  elips='Fisher 1960'
elseif(answer=='I'.or.answer=='i') then
  a=6378150.d0
  f=1.d0/298.3d0
  elips='Fisher 1968'
elseif(answer=='J'.or.answer=='j') then
  a=6378270.d0
  f=1.d0/297.d0
  elips='Hough 1956'
elseif(answer=='K'.or.answer=='k') then
  a=6378388.d0
  f=1.d0/297.d0
  elips='International (Hayford)'
elseif(answer=='L'.or.answer=='l') then
  a=6378245.d0
  f=1.d0/298.3d0
  elips='Krassovsky 1938'
elseif(answer=='M'.or.answer=='m') then
  a=6378145.d0
  f=1.d0/298.25d0
  elips='NWL-9D  (WGS 66)'
elseif(answer=='N'.or.answer=='n') then
  a=6378160.d0
  f=1.d0/298.25d0
  elips='South American 1969'
elseif(answer=='O'.or.answer=='o') then
  a=6378136.d0
  f=1.d0/298.257d0
  elips='Soviet Geod. System 1985'
elseif(answer=='P'.or.answer=='p') then
  a=6378135.d0
  f=1.d0/298.26d0
  elips='WGS 72'
else
  elips = 'User defined.'
!
  write(*,*) '  Enter Equatorial axis,   a : '
  read(*,*) a
  a = abs(a)
!
  write(*,*) '  Enter either Polar axis, b or '
  write(*,*) '  Reciprocal flattening,   1/f : '
  read(*,*) ss
  ss = abs(ss)
!
  f = 0.0_wp
  if( 200.0_wp<=ss .and. ss<=310.0_wp )then
    f = 1.d0/ss
  elseif( 6000000.0_wp<ss .and. ss<a )then
    f = (a-ss)/a
  else
    elips = 'Error: default GRS80 used.'
    a     = 6378137.0_wp
    f     = 1.0_wp/298.25722210088d0
  endif
endif
!
return
end

subroutine fixdms (ideg, min, sec, tol )
!
implicit real(wp) (a-h, o-z)
implicit integer (i-n)
!
!     test for seconds near 60.0-tol
!
if( sec>=( 60.0_wp-tol ) )then
  sec  = 0.0_wp
  min  = min+1
endif
!
!     test for minutes near 60
!
if( min>=60 )then
  min  = 0
  ideg = ideg+1
endif
!
!     test for degrees near 360
!
if( ideg>=360 )then
  ideg = 0
endif
!
return
end

subroutine hem_ns ( lat_sn, hem )
implicit integer (i-n)
character*6  hem
!
if( lat_sn==1 ) then
  hem = 'North '
else
  hem = 'South '
endif
!
return
end

subroutine hem_ew ( lon_sn, hem )
implicit integer (i-n)
character*6  hem
!
if( lon_sn==1 ) then
  hem = 'East  '
else
  hem = 'West  '
endif
!
return
end

subroutine getdeg(d,m,sec,isign,val)

!** comvert deg, min, sec to degrees

implicit real(wp)(a-h,j-z)

val=(d+m/60.d0+sec/3600.d0)
val=dble(isign)*val

return
end

subroutine gvali4 (buff,ll,vali4,icond)
implicit     integer (i-n)
!
logical      plus,sign,done,error
character*1  buff(*)
character*1  ch
!
!     integer*2    i
!     integer*2    l1
!
integer*4    ich,icond
integer*4    ll
integer*4    vali4
!
l1    = ll
vali4 = 0
icond = 0
plus  = .true.
sign  = .false.
done  = .false.
error = .false.
!
i = 0
10 i = i+1
if( i>l1 .or. done )then
  go to 1000
else
  ch  = buff(i)
  ich = ichar( buff(i) )
endif
!
if(     ch=='+' )then
!
!       enter on plus sign
!
  if( sign )then
    goto 150
  else
    sign = .true.
    goto 10
  endif
elseif( ch=='-' )then
!
!       enter on minus sign
!
  if( sign )then
    goto 150
  else
    sign = .true.
    plus = .false.
    goto 10
  endif
elseif( ch>='0' .and. ch<='9' )then
  goto 100
elseif( ch==' ' )then
!
!       enter on space -- ignore leading spaces
!
  if( .not.sign )then
    goto 10
  else
    buff(i) = '0'
    ich = 48
    goto 100
  endif
elseif( ch==':' )then
!
!       enter on colon -- ignore
!
  if( .not.sign )then
    goto 10
  else
    goto 1000
  endif
elseif( ch=='$' )then
!
!       enter on dollar "$"
!
  done = .true.
  goto 10
else
!
!       something wrong
!
  goto 150
endif
!
!     enter on numeric
!
100 vali4 = 10*vali4+(ich-48)
sign  = .true.
goto 10
!
!     treat illegal character
!
150 buff(i) = '0'
vali4 = 0
icond = 1
!
1000 if( .not.plus )then
  vali4 = -vali4
endif
!
return
end

subroutine gvalr8 (buff,ll,valr8,icond)
implicit     integer (i-n)
!
logical      plus,sign,dpoint,done
!
character*1  buff(*)
character*1  ch
!
!     integer*2    i, ip
!     integer*2    l1
!     integer*2    nn, num, n48
!
integer*4    ich,icond
integer*4    ll
!
real(wp)       ten
real(wp)       valr8
real(wp)       zero
!
data zero,ten/0.0_wp,10.0_wp/
!
n48     =  48
l1      =  ll
icond   =   0
valr8   =  zero
plus    = .true.
sign    = .false.
dpoint  = .false.
done    = .false.
!
!     start loop thru buffer
!
i = 0
10 i = i+1
if( i>l1 .or. done )then
  go to 1000
else
  ch  = buff(i)
  nn  = ichar( ch )
  ich = nn
endif
!
if(     ch=='+' )then
!
!       enter on plus sign
!
  if( sign )then
    goto 150
  else
    sign = .true.
    goto 10
  endif
elseif( ch=='-' )then
!
!       enter on minus sign
!
  if( sign )then
    goto 150
  else
    sign = .true.
    plus = .false.
    goto 10
  endif
elseif( ch=='.' )then
!
!       enter on decimal point
!
  ip     = 0
  sign   = .true.
  dpoint = .true.
  goto 10
elseif( ch>='0' .and. ch<='9' )then
  goto 100
elseif( ch==' ' )then
!
!       enter on space
!
  if( .not.sign )then
    goto 10
  else
    buff(i) = '0'
    ich = 48
    goto 100
  endif
elseif( ch==':' .or. ch=='$' )then
!
!       enter on colon or "$" sign
!
  done = .true.
  goto 10
else
!
!       something wrong
!
  goto 150
endif
!
!     enter on numeric
!
100 sign = .true.
if( dpoint )then
  ip = ip+1
endif
!
num   = ich
valr8 = ten*valr8+dble(float( num-n48 ))
goto 10
!
!     treat illegal character
!
150 buff(i) = '0'
valr8   =  0.0_wp
icond   =  1
!
1000 if( dpoint )then
  valr8 =  valr8/(ten**ip)
endif
!
if( .not.plus )then
  valr8 = -valr8
endif
!
return
end

subroutine todmsp(val,id,im,s,isign)

!** convert position degrees to deg,min,sec
!** range is [-180 to +180]

implicit real(wp)(a-h,o-z)

1 if(val>180) then
  val=val-180-180
  go to 1
endif

2 if(val<-180) then
  val=val+180+180
  go to 2
endif

if(val<0.d0) then
  isign=-1
else
  isign=+1
endif

s=abs(val)
id=int(s)
s=(s-id)*60.d0
im=int(s)
s=(s-im)*60.d0

!** account for rounding error

is=nint(s*1.d5)
if(is>=6000000) then
  s=0.d0
  im=im+1
endif
if(im>=60) then
  im=0
  id=id+1
endif

return
end

subroutine trim (buff,lgh,hem)
!
implicit integer (i-n)
character*1 ch,hem
character*1 buff(*)
integer*4   lgh
!
ibeg = 1
do 10 i=1,50
  if( buff(i)/=' ' )then
    goto 11
  endif
  ibeg = ibeg+1
10 continue
11 continue
if( ibeg>=50 )then
  ibeg = 1
  buff(ibeg) = '0'
endif
!
iend = 50
do 20 i=1,50
  j = 51-i
  if( buff(j)==' ' )then
    iend = iend-1
  else
    goto 21
  endif
20 continue
21 continue
!
ch = buff(ibeg)
if( hem=='N' )then
  if( ch=='N' .or. ch=='n' .or. ch=='+' )then
    hem = 'N'
    ibeg = ibeg+1
  endif
  if( ch=='S' .or. ch=='s' .or. ch=='-' )then
    hem = 'S'
    ibeg = ibeg+1
  endif
!
!       check for wrong hemisphere entry
!
  if( ch=='E' .or. ch=='e' )then
    hem = '*'
    ibeg = ibeg+1
  endif
  if( ch=='W' .or. ch=='w' )then
    hem = '*'
    ibeg = ibeg+1
  endif
elseif( hem=='W' )then
  if( ch=='E' .or. ch=='e' .or. ch=='+' )then
    hem = 'E'
    ibeg = ibeg+1
  endif
  if( ch=='W' .or. ch=='w' .or. ch=='-' )then
    hem = 'W'
    ibeg = ibeg+1
  endif
!
!       check for wrong hemisphere entry
!
  if( ch=='N' .or. ch=='n' )then
    hem = '*'
    ibeg = ibeg+1
  endif
  if( ch=='S' .or. ch=='s' )then
    hem = '*'
    ibeg = ibeg+1
  endif
elseif( hem=='A' )then
  if( .not.('0'<=ch .and. ch<='9') )then
    hem = '*'
    ibeg = ibeg+1
  endif
else
!        do nothing
endif
!
!
do 30 i=ibeg,iend
  ch = buff(i)
!
  if(     ch==':' .or. ch=='.' )then
    goto 30
  elseif( ch==' ' .or. ch==',' )then
    buff(i) = ':'
  elseif( '0'<=ch .and. ch<='9' )then
    goto 30
  else
    buff(i) = ':'
  endif
!
30 continue
!
!     left justify buff() array to its first character position
!     also check for a ":" char in the starting position,
!     if found!!  skip it
!
j  = 0
ib = ibeg
ie = iend
!
do 40 i=ib,ie
  if( i==ibeg .and. buff(i)==':' )then
!
!         move the 1st position pointer to the next char &
!         do not put ":" char in buff(j) array where j=1
!
    ibeg = ibeg+1
    goto 40
  endif
  j = j+1
  buff(j) = buff(i)
40 continue
!
!
lgh = iend-ibeg+1
j   = lgh+1
buff(j) = '$'
!
!     clean-up the rest of the buff() array
!
do 50 i=j+1,50
  buff(i) = ' '
50 continue
!
!     save a maximum of 20 characters
!
if( lgh>20 )then
  lgh = 20
  j   = lgh+1
  buff(j) = '$'
endif
!
return
end

end module ngscommon