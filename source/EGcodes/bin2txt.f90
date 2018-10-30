!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This code reads the binary raypath file compatible with the
!! rj_mcmc tomography code and produces a text output.
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM bin2txt

IMPLICIT NONE
INTEGER, PARAMETER :: 	i5=SELECTED_REAL_KIND(5,10)
CHARACTER (LEN=30) :: 	txtout,binin
INTEGER	:: 		npaths,i,j,npt
REAL(KIND=i5)	::	lat,long

OPEN(UNIT=10,FILE='bin2txt.in',status='old')
READ(10,*)binin
READ(10,*)txtout
CLOSE(10)

OPEN(UNIT=20,FILE=binin,FORM='unformatted',status='unknown')
OPEN(UNIT=10,FILE=txtout,FORM='formatted',status='replace')

READ(20)npaths
WRITE(10,*)npaths

DO i=1,npaths
	READ(20)npt
	WRITE(10,*)npt
	DO j=1,npt
		READ(20)lat,long
		WRITE(10,*)lat,long
	END DO
END DO

CLOSE(10)
CLOSE(20)


END PROGRAM bin2txt
