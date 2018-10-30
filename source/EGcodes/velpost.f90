!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This code reads the binary file of velocity posteriors
!! produced by the rj_mcmc tomography code and writes text
!! files of posterior distributions for certain points.
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM velpost
IMPLICIT NONE
INTEGER :: i,j,ii,jj
INTEGER :: nvt,nvp,nvd,npt

INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i10) :: lat,long,latmax,longmin,plat,plong
INTEGER :: pi,pj
REAL(KIND=i10), DIMENSION(:), ALLOCATABLE :: vels,lats,longs!,plat,plong
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: points
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: posts

CHARACTER (LEN=30) :: ifile
CHARACTER (LEN=30) :: ofile
CHARACTER (LEN=40) :: ofilept
CHARACTER (LEN=30) :: iic


OPEN(UNIT=10,FILE='velpost.in',STATUS='old')

READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a30)')ifile
READ(10,'(a30)')ofile

READ(10,*)npt
ALLOCATE(points(npt,2))

DO ii=1,npt
	READ(10,*)points(ii,1),points(ii,2)
END DO

CLOSE(10)


OPEN(UNIT=20,FILE=ifile,FORM='unformatted',STATUS='old')
READ(20)nvt,nvp
READ(20)latmax,longmin
READ(20)lat,long
READ(20)nvd

ALLOCATE(vels(nvd))
DO jj=1,nvd
	READ(20)vels(jj)
END DO

ALLOCATE(longs(nvp))
ALLOCATE(lats(nvt))
DO i=1,nvp
	longs(i)=longmin+(i-1)*long
END DO
DO j=1,nvt
	lats(j)=latmax-(nvt-j)*lat
END DO

ALLOCATE(posts(nvp,nvt,nvd))
DO i=1,nvp
	DO j=1,nvt
		DO jj=1,nvd
			READ(20)posts(i,j,jj)
		END DO
	END DO
END DO
CLOSE(20)

ofile=ADJUSTL(ofile)
DO ii=1,npt

	WRITE(iic,*)ii
	iic=ADJUSTL(iic)
	ofilept=trim(ofile) // '_point_' // trim(iic) // '.out'
	
	OPEN(UNIT=54,FILE=ofilept,STATUS='replace')
	
	plat=points(ii,1)
	pj=MINLOC(ABS(lats-plat),DIM=1)
	plong=points(ii,2)
	pi=MINLOC(ABS(longs-plong),DIM=1)
	
	DO jj=1,nvd
		WRITE(54,*)vels(jj),posts(pi,pj,jj)
	END DO
	
	CLOSE(54)
		
END DO

END PROGRAM velpost
