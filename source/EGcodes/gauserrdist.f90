!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program creates random Gaussian noise dependent on 
!! distance and adds it to an existing traveltime file.
!! The standard deviation of the added noise is given by 
!! sigma=a*distance+b, where distance is either the length of 
!! the raypath (defined in the input file 'raylength') or the 
!! source-receiver distance, and 'a' and 'b' are input parameters.
!! N.B. The distance is give in DEGREES!!!
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM gauserrdist

IMPLICIT NONE
REAL, EXTERNAL :: GASDEV,ran3
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER seed,nsrc,nrec,ns,nr,nn,nrays,idum,distdep
CHARACTER (LEN=30) gefile,sources,receivers,otimes,raylength
REAL(KIND=i10) er,dist,a,b,sdan,err
REAL(KIND=i10), PARAMETER :: pi = 3.1415926535898
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: src,rec
INTEGER, DIMENSION(:), ALLOCATABLE :: raystat
REAL(KIND=i10), DIMENSION(:), ALLOCATABLE :: ttimes


OPEN(UNIT=10,FILE='gauserrdist.in',status='old')
READ(10,1)gefile
READ(10,1)sources
READ(10,1)receivers
READ(10,1)otimes
READ(10,*)distdep
READ(10,1)raylength
READ(10,*)seed
READ(10,*)er
READ(10,*)a,b
1   FORMAT(a30)
CLOSE(10)

OPEN(UNIT=40,FILE=sources,status='old')
OPEN(UNIT=50,FILE=receivers,status='old')
READ(40,*)nsrc
READ(50,*)nrec
ALLOCATE(src(nsrc,2))
ALLOCATE(rec(nrec,2))
DO ns=1,nsrc
  READ(40,*)src(ns,1),src(ns,2)
END DO
DO nr=1,nrec
  READ(50,*)rec(nr,1),rec(nr,2)
END DO
CLOSE(40)
CLOSE(50)

nrays=nsrc*nrec

OPEN(UNIT=60,FILE=otimes,status='old')
ALLOCATE(raystat(nrays))
ALLOCATE(ttimes(nrays))
DO nn=1,nrays
  READ(60,*)raystat(nn),ttimes(nn)
END DO
CLOSE(60)

src=src*pi/180
rec=rec*pi/180

OPEN(UNIT=30,FILE=gefile,status='replace')
OPEN(UNIT=20,FILE=otimes,status='replace')
IF (distdep.EQ.1) OPEN(UNIT=10,FILE=raylength,status='old')
nn=0
DO ns=1,nsrc
  DO nr=1,nrec
    nn=nn+1
    IF (distdep.EQ.1) THEN
      READ(10,*)idum,dist
    ELSE IF (distdep.EQ.0) THEN
      dist=2*180/pi*asin(sqrt((sin((src(ns,1)-rec(nr,1))/2))**2+(cos(rec(nr,1)))*(cos(src(ns,1)))*(sin((src(ns,2)-rec(nr,2))/2))**2))
    ELSE
      WRITE(*,*)'ERROR: invalid input at line 5!!'
      STOP
    END IF
    sdan=a*dist+b
    err=GASDEV(seed)*sdan
    IF (raystat(nn).EQ.1) THEN
      WRITE(30,*)err,sdan
      WRITE(20,*)raystat(nn),ttimes(nn)+err,sdan
    ELSE
      WRITE(30,*)0,0
      WRITE(20,*)0,100.0,0.1
    END IF
  END DO
END DO
CLOSE(30)
CLOSE(20)
IF (distdep.EQ.1) CLOSE(10)

DEALLOCATE(src,rec)
DEALLOCATE(raystat)
DEALLOCATE(ttimes)

END PROGRAM



!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator for 
!       a Gaussian distribution
!
! ------------------------------------------------------------------

      FUNCTION GASDEV(idum)

!     ..Arguments..
      INTEGER          idum
      REAL GASDEV

!     ..Local..
      REAL v1,v2,r,fac
      REAL ran3

      IF (idum.lt.0) iset=0
10     v1=2*ran3(idum)-1
       v2=2*ran3(idum)-1
       r=v1**2+v2**2
       IF(r.ge.1.or.r.eq.0) GOTO 10
       fac=sqrt(-2*log(r)/r)
       GASDEV=v2*fac

      RETURN
      END
!
!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator
!
! ------------------------------------------------------------------
      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
!      WRITE(*,*)' idum ',idum
      IF(idum.lt.0.or.iff.eq.0)THEN
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        DO 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          IF(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      CONTINUE
        DO 13 k=1,4
          DO 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            IF(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        CONTINUE
13      CONTINUE
        inext=0
        inextp=31
        idum=1
      END IF
      inext=inext+1
      IF(inext.eq.56)inext=1
      inextp=inextp+1
      IF(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      IF(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      RETURN
      END