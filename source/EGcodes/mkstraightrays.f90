!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This code crates simple 2-point raypaths with the receiver
!! as first point and the source as the second point.
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM mkstraightrays

IMPLICIT NONE
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER :: checkstat
CHARACTER (LEN=30) :: sources,receivers,txtout,binout
INTEGER	:: nsrc,nrc,npaths,i,j
REAL(KIND=i5), DIMENSION (:), ALLOCATABLE :: scx,scz,rcx,rcz

OPEN(UNIT=10,FILE='mkstraightrays.in',status='old')
READ(10,*)sources
READ(10,*)receivers
READ(10,*)txtout
READ(10,*)binout
CLOSE(10)
!
! Read in all source coordinates.
!
Open(UNIT=10,FILE=sources,STATUS='old')
READ(10,*)nsrc
ALLOCATE(scx(nsrc),scz(nsrc), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL scx,scz'
ENDIF
DO i=1,nsrc
   READ(10,*)scx(i),scz(i)
ENDDO
CLOSE(10)

!
! Read in all receiver coordinates if required
!
OPEN(UNIT=10,FILE=receivers,status='old')
READ(10,*)nrc
ALLOCATE(rcx(nrc),rcz(nrc), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin2d: REAL rcx,rcz'
ENDIF
DO i=1,nrc
   READ(10,*)rcx(i),rcz(i)
END DO
CLOSE(10)

npaths=nsrc*nrc

!
! Write raypath files
!
OPEN(UNIT=10,FILE=txtout,FORM='formatted',status='replace')
OPEN(UNIT=20,FILE=binout,FORM='unformatted',status='replace')

WRITE(10,*)npaths
WRITE(20)npaths

DO i=1,nsrc
   DO j=1,nrc
   
      WRITE(10,*)2
      WRITE(10,*)rcx(j),rcz(j)
      WRITE(10,*)scx(i),scz(i)
      
      WRITE(20)2
      WRITE(20)rcx(j),rcz(j)
      WRITE(20)scx(i),scz(i)
   
   END DO
END DO 

CLOSE(10)
CLOSE(20)

DEALLOCATE(scx,scz,rcx,rcz)

END PROGRAM mkstraightrays
