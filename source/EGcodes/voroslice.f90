!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TYPE: PROGRAM
!! CODE: FORTRAN 90
!! This program is designed to convert velocity from a file 
!! containing coordinates and velocities of Voronoi cells into
!! a form suitable for input to GMT.
!!
!! Adapted from tslicess.f90 in the FMST code package to plot 
!! Voronoi cells avoiding interpolation.
!!
!! Erica Galetti, April 2014
!! erica.galetti@ed.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM slice
IMPLICIT NONE
INTEGER :: i,j
INTEGER :: checkstat
INTEGER :: nnt,nnp
INTEGER :: nnx,nnz
INTEGER :: ddt,ddp
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i10) :: got,gop,rgst,rgsp
REAL(KIND=i10) :: lft,rgt,btm,top
REAL(KIND=i10) :: rd1,rd2
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
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE 	:: 	voro	
INTEGER, DIMENSION(:,:), ALLOCATABLE 		::	vertices,neighbour
INTEGER, DIMENSION(:), ALLOCATABLE 		::	nnn
INTEGER, DIMENSION(:), ALLOCATABLE 		::	nnlist,ntwork
INTEGER, DIMENSION(:), ALLOCATABLE 		::     	worki1,worki2,worki3
LOGICAL, DIMENSION(:), ALLOCATABLE 		:: 	ldummy


OPEN(UNIT=10,FILE='voroslice.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a30)')ifilev
!
! Bounding box of plot
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a22)')ofileb
!
! Now read in velocity grid parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ddt,ddp
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

READ(20,*)nnt,nnp
READ(20,*)got,gop
READ(20,*)rgst,rgsp
READ(20,*)ncell

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
! Calculate GMT bounds file for depth slice if required.
!
lft=gop
rgt=gop+(nnp-1)*rgsp
btm=got-(nnt-1)*rgst
top=got
rd1=rgsp/ddp
rd2=rgst/ddt
OPEN(UNIT=50,FILE=ofileb,STATUS='unknown')
WRITE(50,'(f16.10)')lft
WRITE(50,'(f16.10)')rgt
WRITE(50,'(f16.10)')btm
WRITE(50,'(f16.10)')top
WRITE(50,'(f16.10)')rd1
WRITE(50,'(f16.10)')rd2
CLOSE(50)

!
! Extract velocity slice 
!
nnx=(nnt-1)*ddt+1
nnz=(nnp-1)*ddp+1

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
DO i=1,nnz
	DO j=nnx,1,-1
       
	p2=got-(j-1)*rd2
	p1=gop+(i-1)*rd1

	node=1
			
	CALL find_node2D(DBLE([p1,p2]),node,DBLE(voro(1:2,:)),nnn,nnlist,walk)
			
	WRITE(30,*)voro(3,node)
 
	ENDDO
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

