!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module. 
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! 	Code modified by Erica Galetti (EG) in order to produce a text raypath file called raypath.out compatible with the RJ_MCMC tomography code.
!!!!!	[In the original version of fm2dss.f90, the raypath.out file was binary and incompatible with the RJ_MCMC code.]

MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER :: checkstat
INTEGER, SAVE :: nvx,nvz,nnx,nnz,nrc,fom,gdx,gdz
INTEGER, SAVE :: vnl,vnr,vnt,vnb,nrnx,nrnz,sgdl,rbint
INTEGER, SAVE :: nnxr,nnzr,asgr
INTEGER, DIMENSION (:,:), ALLOCATABLE :: nsts,nstsr,srs
REAL(KIND=i10), SAVE :: gox,goz,dnx,dnz,dvx,dvz,snb,earth
REAL(KIND=i10), SAVE :: goxd,gozd,dvxd,dvzd,dnxd,dnzd
REAL(KIND=i10), SAVE :: drnx,drnz,gorx,gorz
REAL(KIND=i10), SAVE :: dnxr,dnzr,goxr,gozr
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: velv,veln,velnb
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE, SAVE :: ttn,ttnr
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE, SAVE :: rcx,rcz
REAL(KIND=i10), PARAMETER :: pi=3.1415926535898

INTEGER, DIMENSION (:,:), ALLOCATABLE :: srsv
INTEGER, DIMENSION (:,:), ALLOCATABLE :: npts
REAL(KIND=i5), DIMENSION (:,:), ALLOCATABLE :: raypts
INTEGER	:: crazyrp
!
! nvx,nvz = B-spline vertex values
! dvx,dvz = B-spline vertex separation
! velv(i,j) = velocity values at control points
! nnx,nnz = Number of nodes of grid in x and z
! nnxr,nnzr = Number of nodes of refined grid in x and z
! gox,goz = Origin of grid (theta,phi)
! goxr, gozr = Origin of refined grid (theta,phi)
! dnx,dnz = Node separation of grid in  x and z
! dnxr,dnzr = Node separation of refined grid in x and z
! veln(i,j) = velocity values on a refined grid of nodes
! velnb(i,j) = Backup of veln required for source grid refinement
! ttn(i,j) = traveltime field on the refined grid of nodes
! ttnr(i,j) = ttn for refined grid
! nsts(i,j) = node status (-1=far,0=alive,>0=close)
! nstsr(i,j) = nsts for refined grid
! checkstat = check status of memory allocation
! fom = use first-order(0) or mixed-order(1) scheme
! snb = Maximum size of narrow band as fraction of nnx*nnz
! nrc = number of receivers
! rcx(i),rcz(i) = (x,z) coordinates of receivers
! earth = radius of Earth (in km)
! goxd,gozd = gox,goz in degrees
! dvxd,dvzd = dvx,dvz in degrees
! dnzd,dnzd = dnx,dnz in degrees
! gdx,gdz = grid dicing in x and z
! vnl,vnr,vnb,vnt = Bounds of refined grid
! nrnx,nrnz = Number of nodes in x and z for refined grid
! gorx,gorz = Grid origin of refined grid
! sgdl = Source grid dicing level
! rbint = Ray-boundary intersection (0=no, 1=yes).
! asgr = Apply source grid refinement (0=no,1=yes)
! srs = Source-receiver status (0=no path, 1=path exists)
!
END MODULE globalp
