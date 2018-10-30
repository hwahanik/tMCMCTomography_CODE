
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!                   MODULE                    !!!!!!!!!!!!!! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%f
MODULE fm2dssvars

!!! This module includes variables that are used when updating raypaths using the subroutines in library updrays
	IMPLICIT NONE
	
	! Define two types of real numbers (for consistency with FMST code)
	INTEGER, PARAMETER :: ii5=SELECTED_REAL_KIND(5,10)	! Added by EG to make variables compatible with the fm2dssRJ code
	INTEGER, PARAMETER :: ii10=SELECTED_REAL_KIND(10,100)	! Added by EG to make variables compatible with the fm2dssRJ code
	
	REAL(KIND=ii10), PARAMETER :: pii = 3.1415926535898
	
	! Latitude and longitude range and dicing
	REAL(KIND=ii10), SAVE	:: longmin,longmax
  	REAL(KIND=ii10), SAVE	:: latmin,latmax
	REAL(KIND=ii10), SAVE	:: latr,longr
	
	! Number of velocity nodes in theta (latitude) and phi (longitude)
	INTEGER, SAVE 	:: nvt,nvp,nvtr,nvpr
	
	! Parameters for raypath modelling
	INTEGER, SAVE 		:: gdt,gdp,sgref,dl,erg,fms,sigdep
	REAL(KIND=ii10), SAVE	:: er,nbs
	
	! Number and location of sources and receivers
	INTEGER, SAVE 	:: nsou,nrec
	REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE, SAVE :: 	sou,rec
	
	! Ray info
	INTEGER, SAVE 	:: nrays
	INTEGER		:: wurf,test
	INTEGER, SAVE 	:: uar
	INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE 	:: 	raystat
	INTEGER, DIMENSION(:), ALLOCATABLE 		::	lengthrayv
	REAL(KIND=ii10), DIMENSION(:), ALLOCATABLE 	:: 	dis,srdist,srdist_prop
	
	! Data
	REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE, SAVE :: 	dat
	
	! Variables for Delaunay triangulation (construction of Voronoi cells)
	INTEGER, SAVE	:: ncell_min,ncell_max
	INTEGER, SAVE	:: nt_max,nv_max,nmax,ptsok
	REAL(KIND=ii10)	:: eps = 0.000001
	INTEGER		:: nt,walk
	INTEGER, DIMENSION(:,:), ALLOCATABLE 	::	vertices,neighbour
	INTEGER, DIMENSION(:), ALLOCATABLE 	::	nnn
	INTEGER, DIMENSION(:), ALLOCATABLE 	::	nnlist,ntwork
	INTEGER, DIMENSION(:), ALLOCATABLE 	::	worki1,worki2,worki3
	LOGICAL, DIMENSION(:), ALLOCATABLE 	::	ldummy
	
	CHARACTER (LEN=30) ::	errout
	CHARACTER (LEN=50) ::	errvtx,errvor
	INTEGER	:: sse

END MODULE fm2dssvars


