!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program calculates raypaths through a Voronoi-tessellated
!! 2D space and outputs the following information:
!! - rayjour.dat: each line in this file corresponds to the path
!!                taken by each raypath in Voronoi cell numbers, 
!!                with the ray traced from the receiver to the 
!!                source
!! - raylengths.dat: this file gives the length of each raypath
!!                   in each Voronoi cell
!! - raycell.dat: this file indicates if a ray passes through a
!!                certain cell with a 1, and if it does not with
!!                a 0
!! - rayneigh.dat: this file contains the number of the raypaths
!!                 passing through the cell indicated at line 30
!!                 of the input file and its neighbours
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE fm2dssvars

!!! This module includes variables that are used when calculating raypaths
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



PROGRAM vororay

USE fm2dssvars
USE fm2dout
IMPLICIT NONE

INTEGER :: i,j,ii,jj,iii,jjj,nr
INTEGER :: checkstat
REAL(KIND=ii10) :: dum1,dum2
INTEGER :: dum3,dum4
CHARACTER (LEN=30) :: ifilev,sources,receivers,otimes
CHARACTER (LEN=30) :: ofilev,rayjourfile,raylengthfile,raycelfile,rpfile,rayneighfile

REAL(KIND=ii10)	p1,p2,vlat,vlong,vvel
INTEGER		ncell,node,ptloc
INTEGER		p,nrr,s,r,v,wvgf,nray_max
INTEGER		istart,iend,k,ind
REAL(KIND=ii10) n,t
INTEGER		rayid,raylength,pti,ptf
REAL(KIND=ii10)	tetaS,tetaR,phiS,phiR,alpha,angle
REAL(KIND=ii10)	point(2)
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE	::	VEL
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE	::	dato
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE	::	RayCel
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE	::	Mi
REAL(KIND=ii5), DIMENSION(:,:), ALLOCATABLE	::	raypoints
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 	::	voro
LOGICAL, DIMENSION(:), ALLOCATABLE 		::	cellogic


OPEN(UNIT=10,FILE='vororay.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a30)')ifilev
READ(10,'(a30)')sources
READ(10,'(a30)')receivers
READ(10,'(a30)')otimes
!
! Area (velocity grid parameters)
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)latmax,latmin
READ(10,*)longmin,longmax
READ(10,*)nvtr,nvpr
READ(10,*)gdt,gdp
READ(10,*)sgref
READ(10,*)dl,erg
READ(10,*)er
READ(10,*)fms
READ(10,*)nbs
!
! Output files
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a30)')rayjourfile
READ(10,'(a30)')raylengthfile
READ(10,'(a30)')raycelfile
READ(10,*)wvgf
READ(10,'(a30)')ofilev
READ(10,*)wurf
READ(10,'(a30)')rpfile
READ(10,*)ind
READ(10,'(a30)')rayneighfile
!
! Now read in parameters for the Voronoi tesselation
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ncell_max
READ(10,*)nt_max
READ(10,*)nv_max
READ(10,*)nray_max
CLOSE(10)

! Allocate arrays
nmax=3*nt_max+ncell_max
ALLOCATE(vertices(3,nt_max))
ALLOCATE(neighbour(3,nt_max))
ALLOCATE(nnn(ncell_max+1))
ALLOCATE(nnlist(nmax))
ALLOCATE(ntwork(nmax))
ALLOCATE(worki1(nv_max))
ALLOCATE(worki2(nv_max))
ALLOCATE(worki3(nv_max))
ALLOCATE(ldummy(ncell_max))

!
! Calculate GMT bounds files. Start off by reading in
! velocity grid.
!
OPEN(UNIT=20,FILE=ifilev,status='old')

READ(20,*)dum3,dum4
READ(20,*)dum1,dum2
READ(20,*)dum1,dum2
READ(20,*)
READ(20,*)ncell

WRITE(*,*)'Model has',ncell,'cells.'

ALLOCATE(voro(3,ncell), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM voroslice: REAL voro'
ENDIF

DO i=1,ncell
	READ(20,*)vlat,vlong,vvel
	voro(1,i)=vlong
	voro(2,i)=vlat
	voro(3,i)=vvel
END DO

CLOSE(20)


! Read file of sources
OPEN(UNIT=2,FILE=sources,STATUS='old')
READ(2,*) nsou
! sou(:,1) is latitude and sou(:,2) is longitude
ALLOCATE(sou(nsou,2))
DO i=1,nsou
	READ(2,*)sou(i,1),sou(i,2)
	IF(sou(i,1).GT.latmax.OR.sou(i,1).LT.latmin.OR.sou(i,2).GT.longmax.OR.sou(i,2).LT.longmin) THEN
		WRITE(*,*)'Source',i,'is located outside model boundaries!!'
		WRITE(*,*)'Terminating program!!'
		STOP
	END IF
END DO
CLOSE(2)

! Read file of receivers
OPEN(UNIT=1,FILE=receivers,STATUS='old')
READ(1,*) nrec
! rec(:,1) is latitude and rec(:,2) is longitude
ALLOCATE(rec(nrec,2))
DO i=1,nrec
	READ(1,*)rec(i,1),rec(i,2)
	IF(rec(i,1).GT.latmax.OR.rec(i,1).LT.latmin.OR.rec(i,2).GT.longmax.OR.rec(i,2).LT.longmin) THEN
		WRITE(*,*)'Receiver',i,'is located outside model boundaries!!'
		WRITE(*,*)'Terminating program!!'
		STOP
	END IF
END DO
CLOSE(1)

! Allocate arrays
ALLOCATE(raystat(nsou*nrec,2)) ! raystat contains the switch in column 1 and the valid ray number in column 2
ALLOCATE(dato(nray_max,6))
raystat=0
dato=0
! Read file of traveltime associations
OPEN(UNIT=3,FILE=otimes,STATUS='old')
p = 0
nrr = 0 ! "loops" over rays
DO s=1,nsou
   DO r=1,nrec
        nrr = nrr + 1 
	READ(3,*)v,t,n   ! read validity of the data (switch), traveltime , noise (error in traveltime)
	IF (v.ne.0) THEN
		p = p+1 
                dato(p,1:2)=sou(s,1:2)
                dato(p,3:4)=rec(r,1:2)
                dato(p,5)=t
		dato(p,6)=n
		raystat(nrr,1)=1
		raystat(nrr,2)=p
	ELSE 
		raystat(nrr,1)=0
		raystat(nrr,2)=0
	END IF
   END DO
END DO
CLOSE(3)! close the file

nrays = p
ALLOCATE(dat(nrays,6))
dat=dato(1:nrays,:)
DEALLOCATE(dato)

!
! Node spacing
!
longr=(longmax-longmin)/(nvpr-1)
latr=(latmax-latmin)/(nvtr-1)
! Allocate array for velocity grid
IF(ALLOCATED(VEL)) DEALLOCATE(VEL)
ALLOCATE(VEL(nvpr+2,nvtr+2))


! Calculate Delaunay triangulation
ptsok=0
CALL delaun(DBLE(voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
                      worki1,worki2,worki3,eps,nv_max,&
                      0,ldummy,0,0,0,ptsok)
IF (ptsok.EQ.1) THEN
	WRITE(*,*)'All',ncell,'points are in a line.'
	WRITE(*,*)'Terminating program!'
	STOP
END IF
CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
              neighbour,nnn,nnlist,ntwork)

IF(wvgf.EQ.1) THEN
	OPEN(UNIT=30,FILE=ofilev,STATUS='unknown')
	WRITE(30,*)nvtr,nvpr
	WRITE(30,*)latmax,longmin
	WRITE(30,*)latr,longr
	WRITE(30,'(1X)')
END IF
iii=0
DO i=0,nvpr+1
	jjj=0
	iii=iii+1
	DO j=0,nvtr+1 !nvtr+1,0,-1
		jjj=jjj+1
        
		ii=i
		jj=j
		IF  (i.EQ.0) ii=1
		IF  (j.EQ.0) jj=1
		IF  (i.EQ.nvpr+1) ii=nvpr
		IF  (j.EQ.nvtr+1) jj=nvtr
		
		p2=latmax-(jj-1)*latr	! theta (latitude)
		p1=longmin+(ii-1)*longr	! phi (longitude)
		
		node=1
				
		CALL find_node2D(DBLE([p1,p2]),node,DBLE(voro(1:2,:)),nnn,nnlist,walk)
			
		IF(wvgf.EQ.1) WRITE(30,*)voro(3,node)
		VEL(iii,jjj)=voro(3,node)
	ENDDO
	IF(wvgf.EQ.1) WRITE(30,'(1X)')
ENDDO
IF(wvgf.EQ.1) CLOSE(30)

!
! Model the raypaths
!
IF(ALLOCATED(npoints)) DEALLOCATE(npoints)
IF(ALLOCATED(raycoords)) DEALLOCATE(raycoords)
IF(ALLOCATED(raypoints)) DEALLOCATE(raypoints)
crazyray=0
nru=sum(raystat(:,1),1)
IF(nru.NE.nrays) THEN
	WRITE(*,*)'There is a problem with the number of active rays!!'
	STOP
END IF
WRITE(*,*)			
WRITE(*,*)'Calling MODRAYS to model raypaths'
WRITE(*,*)
CALL modrays(nsou,sou(:,1),sou(:,2),&
	nrec,rec(:,1),rec(:,2),&
	raystat,-1,&
	nvtr,nvpr,latmax,longmin,latr,longr,VEL,&
	gdt,gdp,&
	sgref,&
	dl,erg,&
	er,&
	fms,&
	nbs,&
	1)
IF(crazyray.GT.0) THEN
	WRITE(*,*)'There is a crazy ray!!'
	STOP
END IF

! Calculate the length of each path in each cell - fill up array RayCel	
IF(ALLOCATED(Mi)) DEALLOCATE(Mi)
IF(ALLOCATED(RayCel)) DEALLOCATE(RayCel)
ALLOCATE(lengthrayv(nrays))
ALLOCATE(dis(nrays))
ALLOCATE(srdist(nrays))
ALLOCATE(RayCel(nrays,ncell_max))
lengthrayv=0
dis=0
RayCel=0


! Write output files
OPEN(UNIT=50,FILE=rayjourfile,STATUS='replace')
IF(wurf.EQ.1) THEN
	OPEN(UNIT=40,FILE=rpfile,FORM='unformatted',STATUS='replace')
	WRITE(40)nru
END IF

pti=0
ptf=0
DO nr=1,nru
	rayid=npoints(nr,2)
	raylength=npoints(nr,1)

	pti=pti+1
	ptf=pti+raylength-1
	ALLOCATE(raypoints(raylength,2))
	raypoints=raycoords(pti:ptf,:)
	pti=ptf
	ALLOCATE(Mi(raylength-1,3))
	lengthrayv(rayid)=raylength
	
	angle=0
	
	point(2) = raypoints(1,1)	! latitude
	point(1) = raypoints(1,2)	! longitude
	point(1) = point(1)*pii/180
	point(2) = point(2)*pii/180
	
	p = 0
	Mi = 0
	
	RayCel(rayid,:) = 0
	 
	IF(wurf.EQ.1) THEN
		WRITE(40)raylength
		WRITE(40)raypoints(1,1),raypoints(1,2)
	END IF
	
	DO i=2,raylength
	
		IF(wurf.EQ.1) THEN
			WRITE(40)raypoints(i,1),raypoints(i,2)
		END IF
	
		p=p+1
		Mi(p,2) = raypoints(i,1)	! latitude
		Mi(p,1) = raypoints(i,2)	! longitude
		Mi(p,1) = Mi(p,1)*pii/180
		Mi(p,2) = Mi(p,2)*pii/180
		tetaS = Mi(p,2)			! latitude
		phiS = Mi(p,1)			! longitude
		IF (i.ne.2) THEN
			tetaR = Mi(p-1,2)	! latitude
			phiR = Mi(p-1,1)	! longitude
		ELSE
			tetaR = point(2)	! latitude
			phiR = point(1)		! longitude
		END IF 
		
		! Do some trigonometry on the sphere to get angle alpha between two consecutive points.
		! Use the Haversine formula (better conditioned for small distances than the spherical law of cosines)
		alpha=2*asin(sqrt((sin((tetaS-tetaR)/2))**2+(cos(tetaR))*(cos(tetaS))*(sin((phiS-phiR)/2))**2))
		
		
		! Store the distance between 2 points in Mi(:,3)
		Mi(p,3)=alpha*er
		angle = angle + alpha
		p1 = (180/pii)*Mi(p,1)
		p2 = (180/pii)*Mi(p,2)
		
		! Find in which Voronoi cell the current point is (i.e. to which node the current point is closest)
		CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro(1:2,:)),nnn,nnlist,walk)
		
		! Compute the distance of ray 'rayid' in cell 'node' and store it in RayCel(rayid,node)
		RayCel(rayid,node)=RayCel(rayid,node)+Mi(p,3)
		
		IF (i.EQ.2) ptloc=node
		IF (ptloc.NE.node) THEN
			IF (i.LT.raylength) WRITE(50,'(1I5)',advance='no')ptloc 
			IF (i.EQ.raylength) WRITE(50,'(1I5)')ptloc
			ptloc=node
		ELSE IF (ptloc.EQ.node.AND.i.EQ.raylength) THEN
			WRITE(50,'(1I5)')ptloc
		END IF
		
	END DO
	
	dis(rayid)=angle*er
	
	DEALLOCATE(Mi)
	
	DEALLOCATE(raypoints)
	
END DO
CLOSE(50)
IF(wurf.EQ.1) THEN
	CLOSE(40)
END IF

DEALLOCATE(raycoords)
DEALLOCATE(npoints)
DEALLOCATE(VEL)

! Find the neighbours of cell 'ind'
ALLOCATE(cellogic(ncell))
 cellogic=.false.
 cellogic(ind)=.true.
istart = nnn(ind)
iend = nnn(ind+1)-1
DO k=istart,iend
	IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
END DO
WRITE(*,*)'The following cells are neighbours of cell',ind,':'
DO i=1,ncell
	IF (cellogic(i).EQ..true..AND.i.NE.ind) WRITE(*,*)i
END DO
WRITE(*,*)

! Write output files
OPEN(UNIT=60,FILE=raylengthfile,STATUS='replace')
OPEN(UNIT=70,FILE=raycelfile,STATUS='replace')
OPEN(UNIT=80,FILE=rayneighfile,STATUS='replace')
DO nr=1,nrays
	DO i=1,ncell
		! Ray lenghts file
		IF (i.LT.ncell) WRITE(60,'(1f10.4)',advance='no')RayCel(nr,i)
		IF (i.EQ.ncell) WRITE(60,'(1f10.4)')RayCel(nr,i)
		! Ray-cell associations file
		IF (RayCel(nr,i).NE.0.AND.i.LT.ncell) WRITE(70,'(1I5)',advance='no') INT(1)
		IF (RayCel(nr,i).NE.0.AND.i.EQ.ncell) WRITE(70,'(1I5)') INT(1)
		IF (RayCel(nr,i).EQ.0.AND.i.LT.ncell) WRITE(70,'(1I5)',advance='no') INT(0)
		IF (RayCel(nr,i).EQ.0.AND.i.EQ.ncell) WRITE(70,'(1I5)') INT(0)
	END DO
END DO
DO i=1,ncell
	IF (cellogic(i).EQ..true..AND.SUM(RayCel(:,i)).NE.0) write(80,'(A,1I5,A)',advance='no')'Cell ',i,'   rays: '
	IF (cellogic(i).EQ..true..AND.SUM(RayCel(:,i)).EQ.0) write(80,'(A,1I5,A)')'Cell ',i,'   empty'
	DO nr=1,nrays
		! Rays-through-cell-'ind' file
		IF (cellogic(i).EQ..true..AND.RayCel(nr,i).NE.0.AND.SUM(RayCel(nr+1:nrays,i)).NE.0) write(80,'(1I5)',advance='no') nr
		IF (cellogic(i).EQ..true..AND.RayCel(nr,i).NE.0.AND.SUM(RayCel(nr+1:nrays,i)).EQ.0) write(80,'(1I5)') nr
	END DO
END DO
CLOSE(60)
CLOSE(70)
CLOSE(80)

DEALLOCATE(sou,rec)
DEALLOCATE(raystat)
DEALLOCATE(lengthrayv)
DEALLOCATE(dis,srdist)
DEALLOCATE(dat)
DEALLOCATE(RayCel)

DEALLOCATE(voro,neighbour,vertices,nnn,nnlist,ntwork, STAT=checkstat)
IF(checkstat > 0)THEN
	WRITE(6,*)'Error with DEALLOCATE: PROGRAM voroslice: voro,neighbour,vertices,nnn,nnlist'
ENDIF
DEALLOCATE(worki1,worki2,worki3,ldummy, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM voroslice: worki1,worki2,worki3,ldummy'
ENDIF

STOP
END PROGRAM vororay
