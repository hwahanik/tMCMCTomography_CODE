!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program creates random Gaussian noise that can be
!! added to a traveltime dataset.
!! Noise can be crated independently for different datasets.
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM gauserr

IMPLICIT NONE
REAL, EXTERNAL :: GASDEV,ran3
INTEGER, PARAMETER :: ii5=SELECTED_REAL_KIND(5,10)
INTEGER seed,totr,rr,nds,nss
CHARACTER (LEN=30) dsfile,gefile
REAL(KIND=ii5), DIMENSION(:), ALLOCATABLE :: sigmas



OPEN(UNIT=10,FILE='gauserr.in',status='old')
READ(10,1)dsfile
READ(10,1)gefile
READ(10,*)seed
READ(10,*)totr
READ(10,*)nds
ALLOCATE(sigmas(nds))
DO nss=1,nds
    READ(10,*)sigmas(nss)
END DO
1   FORMAT(a30)
CLOSE(10)

IF(nds.EQ.1) WRITE(*,*)'Create Gaussian error for ',nds,' dataset'
IF(nds.GT.1) WRITE(*,*)'Create Gaussian error for ',nds,' datasets'

OPEN(UNIT=20,FILE=dsfile,status='old')
OPEN(UNIT=30,FILE=gefile,status='replace')
DO rr=1,totr
      READ(20,*)nss
      IF(nds.GT.1.AND.nss.NE.0) WRITE(30,*)GASDEV(seed)*sigmas(nss),sigmas(nss)
      IF(nds.GT.1.AND.nss.EQ.0) WRITE(30,*)0,0
      IF(nds.EQ.1.AND.nss.NE.0) WRITE(30,*)GASDEV(seed)*sigmas(1),sigmas(1)
      IF(nds.EQ.1.AND.nss.EQ.0) WRITE(30,*)0,0
END DO
CLOSE(20)
CLOSE(30)

DEALLOCATE(sigmas)

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