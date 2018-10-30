!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TYPE: PROGRAM
!! CODE: FORTRAN 90
!! This program is designed to convert velocity from a file 
!! containing coordinates and velocities of Voronoi cells into
!! a form suitable for input to fm2dss.
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM slice
IMPLICIT NONE
INTEGER :: i,j,ii,jj
INTEGER :: checkstat
INTEGER :: nnt,nnp
INTEGER :: nnx,nnz
INTEGER :: ddt,ddp
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i10) :: got,gop,rgst,rgsp
REAL(KIND=i10) :: lft,rgt,btm,top
REAL(KIND=i10) :: rd1,rd2
REAL(KIND=i10) :: dum1,dum2
INTEGER :: dum3,dum4
CHARACTER (LEN=30) :: ifilev
CHARACTER (LEN=30) :: ofileb,ofilev
CHARACTER (LEN=30) :: ifiles,ofiles,ifilerc,ofilerc

! Added by EG
REAL(KIND=i10)	p1,p2,vlat,vlong,vvel
INTEGER		ncell
INTEGER		ptsok
! Variables for Delaunay triangulation (construction of Voronoi cells)
INTEGER, SAVE	:: ncell_min,ncell_max
INTEGER, SAVE	:: nt_max,nv_max,nmax
REAL(KIND=i10)	:: eps = 0.000001
INTEGER		:: nt,walk,node
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE 	::	voro	
INTEGER, DIMENSION(:,:), ALLOCATABLE 		::	vertices,neighbour
INTEGER, DIMENSION(:), ALLOCATABLE 		::	nnn
INTEGER, DIMENSION(:), ALLOCATABLE 		::	nnlist,ntwork
INTEGER, DIMENSION(:), ALLOCATABLE 		::	worki1,worki2,worki3
LOGICAL, DIMENSION(:), ALLOCATABLE 		::	ldummy


OPEN(UNIT=10,FILE='vorogrid.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a30)')ifilev
!
! Area (velocity grid parameters)
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)top,btm
READ(10,*)lft,rgt
READ(10,*)nnt,nnp
got=top
gop=lft
nnx=nnt
nnz=nnp
!
! Output file
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a22)')ofilev
!
! Now read in parameters for the Voronoi tesselation
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ncell_max
READ(10,*)nt_max
READ(10,*)nv_max
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

!
! Node spacing
!
rd1=(rgt-lft)/(nnp-1)
rd2=(top-btm)/(nnt-1)


! Calculate Delaunay triangulation
ptsok=0
CALL delaun(DBLE(voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
                      worki1,worki2,worki3,eps,nv_max,&
                      0,ldummy,0,0,0,ptsok)
IF (ptsok.EQ.1) THEN
	WRITE(*,*)'All',ncell,'points are in a line!'
	WRITE(*,*)'Terminating program!'
END IF
CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
              neighbour,nnn,nnlist,ntwork)

OPEN(UNIT=30,FILE=ofilev,STATUS='unknown')
WRITE(30,*)nnx,nnz
WRITE(30,*)got,gop
WRITE(30,*)rd2,rd1
WRITE(30,'(1X)')
DO i=0,nnz+1
	DO j=0,nnx+1 !nnx+1,0,-1
        
	ii=i
	jj=j
        IF  (i.EQ.0) ii=1
	IF  (j.EQ.0) jj=1
	IF  (i.EQ.nnz+1) ii=nnz
	IF  (j.EQ.nnx+1) jj=nnx
	
	p2=got-(jj-1)*rd2	! theta (latitude)
	p1=gop+(ii-1)*rd1	! phi (longitude)

	node=1
			
	CALL find_node2D(DBLE([p1,p2]),node,DBLE(voro(1:2,:)),nnn,nnlist,walk)
			
	WRITE(30,*)voro(3,node)
 
	ENDDO
	WRITE(30,'(1X)')
ENDDO
CLOSE(30)

DEALLOCATE(voro,neighbour,vertices,nnn,nnlist, STAT=checkstat)
IF(checkstat > 0)THEN
	WRITE(6,*)'Error with DEALLOCATE: PROGRAM voroslice: voro,neighbour,vertices,nnn,nnlist'
ENDIF
DEALLOCATE(worki1,worki2,worki3,ldummy, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM voroslice: worki1,worki2,worki3,ldummy'
ENDIF

STOP
END PROGRAM slice

