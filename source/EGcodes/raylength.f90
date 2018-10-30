!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This code reads the binary raypath file compatible with the
!! fm2dss code and produces a text output with raypath lengths.
!! NOTE: All raypaths must be calculated with fm2dss, so the 
!! number of rays in the fm2dss output must match the number
!! of records in the otimes file.
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM raylength

IMPLICIT NONE
INTEGER, PARAMETER 		::	i5=SELECTED_REAL_KIND(5,10)
INTEGER, PARAMETER 		::	i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10), PARAMETER 	:: 	pi = 3.1415926535898
CHARACTER (LEN=30) 		::	txtout,binin,otimes
INTEGER				::	nrays,nr,npoints,np,raypts,raystat
REAL(KIND=i5)			::	ptlat,ptlong
REAL(KIND=i10)			::	Rlat,Rlong,Slat,Slong,tetaR,phiR,tetaS,phiS
REAL(KIND=i10)			::	alpha,angle



OPEN(UNIT=10,FILE='raylength.in',status='old')
READ(10,*)binin
READ(10,*)otimes
READ(10,*)txtout
CLOSE(10)

OPEN(UNIT=30,FILE=otimes,FORM='formatted',status='unknown')
OPEN(UNIT=20,FILE=binin,FORM='unformatted',status='unknown')
OPEN(UNIT=10,FILE=txtout,FORM='formatted',status='replace')

READ(20)nrays
DO nr=1,nrays

  READ(20)npoints
  READ(30,*)raystat
  
  IF (npoints.NE.0.AND.raystat.EQ.1) THEN
  
    angle=0
    READ(20)ptlat,ptlong 	! read coordinates of the first point of the ray
    Rlat = REAL(ptlat,KIND=i10)*pi/180
    Rlong = REAL(ptlong,KIND=i10)*pi/180
      
    DO np=1,npoints-1
      
      READ(20)ptlat,ptlong
      
      Slat = REAL(ptlat,KIND=i10)*pi/180
      Slong = REAL(ptlong,KIND=i10)*pi/180
      tetaS = Slat		! latitude
      phiS = Slong		! longitude
      tetaR = Rlat		! latitude
      phiR = Rlong		! longitude 

      ! Do some trigonometry on the sphere to get angle alpha between two consecutive points.
      ! Use the Haversine formula (better conditioned for small distances than the spherical law of cosines)
      alpha=2*asin(sqrt((sin((tetaS-tetaR)/2))**2+(cos(tetaR))*(cos(tetaS))*(sin((phiS-phiR)/2))**2))

      ! Add alpha to angle
      angle = angle + alpha
      
      ! S becomes R
      Rlat = Slat
      Rlong = Slong
      
    END DO
    
    angle = angle*(180/pi)
    
    WRITE(10,*)raystat,angle
 
  ELSE IF (npoints.EQ.0.AND.raystat.EQ.0) THEN
    WRITE(10,*)raystat,0
  ELSE
    WRITE(*,*)'ERROR: Number of rays from fm2dss does not match number of records in otimes file!!'
    STOP
  END IF

END DO

CLOSE(30)
CLOSE(20)
CLOSE(10)


END PROGRAM raylength
