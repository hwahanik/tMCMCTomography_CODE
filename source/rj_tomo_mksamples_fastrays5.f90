!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This progranm performs transdimensional 2D tomography using the reversible-jump 	!!
!! algorithm by Bodin & Sambridge (2009).						!!
!!											!!
!! All samples accepted along the Markov chain are saved and can be processed using	!!
!! the accompanying code rj_tomo_procsamples.f90 (procsamples).				!!
!!											!!
!! Inputs, outputs and parameters are defined in the input file mksamples.in		!!
!! 											!!
!! This code is based on the rj_tomo.f90 code written by Thomas Bodin, and uses an	!!
!! adapted version of the Eikonal solver fm2dss.f90 by Nick Rawlinson to update raypath	!!
!! geometries within the code.								!!
!!											!!
!! Random numbers are generated using the Mersenne-Twister algorithm that can be found 	!!
!! on the following page: http://jblevins.org/mirror/amiller/mt19937.f90		!!
!! Gaussian distribution sequences are generated using an adapted version of the GASDEV	!!
!! algorithm in Numerical Recipes.							!!
!!											!!
!! Erica Galetti (EG), April 2013							!!
!! Contact: erica.galetti@ed.ac.uk							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI Modification by Helmut Wahanik (HW), Research Scientist at Schlumberger Research !!
!! Version 22 Oct-Nov 2014								!!
!! Contact: HWahanik@slb.com											!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Original code mofified by Helmut Wahanik, Research Scientist, Schlumberger Research.  !!
!											!!
! hybrid embarrasingly-parallel and mpi code						!!
! 											!!
! Original code mofified by Helmut Wahanik (HW).  Tays-routines MPI parallelization	!!
! Here we build 1 Markov chain instead of an mpi parallelization of nbproc processors	!!
! The mpi parallelization will be held at the ray algorithms stage			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!                   MODULES                   !!!!!!!!!!!!!! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        INTEGER :: wrgf ! Write rays to file?
	
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
	INTEGER, DIMENSION(:,:), ALLOCATABLE	 	:: 	raystat
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Parallelization variables

MODULE mksamples_prllvars

  IMPLICIT NONE

  INTEGER, PARAMETER :: iiii5=SELECTED_REAL_KIND(5,10)	        ! Added by EG to make variables compatible with the fm2dssRJ code
  INTEGER, PARAMETER :: iiii10=SELECTED_REAL_KIND(10,100)	! Added by EG to make variables compatible with the fm2dssRJ code
!! Variables used in embarrisingly parallel outer loop for creating chains.

  INTEGER ra, rank, nbproc, ranseed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! For MPI inner loop (rays routines)
 
  INTEGER	:: ierror,rr
  REAL		:: clock
  INTEGER	:: date_time(8)
  CHARACTER*10	:: datevals(3)
  INTEGER	:: wcount
  INTEGER, SAVE	:: run_worker 
  INTEGER	:: run_w
  INTEGER	:: TEST_MODRAYS


!!!! Setting message types for MPI.  Variables introduced by HW 
!!!! for parallelization of raytracer.  

  INTEGER	:: TAG
  INTEGER,SAVE 	:: MASTER        ! taskid of first task 
  INTEGER,SAVE	:: FROM_MASTER   ! setting a message type for works coming from the master
  INTEGER,SAVE	:: FROM_WORKER   ! setting a message type for works coming from the workers

  INTEGER	:: numtasks		! total number of tasks (processes) available from mpi partition 
  INTEGER	:: taskid		! task identifier
  INTEGER	:: numworkers		! number of worker tasks available for the ray-routines processing
  INTEGER	:: from_work		! task id of message source 
  INTEGER	:: dest			! task id of message destination
  INTEGER	:: origin
 
  INTEGER	:: counter_1, counter_2, counter_3, counter_4, counter_5    ! counters for MPI send.
  !INTEGER, PARAMETER :: iii5=SELECTED_REAL_KIND(5,10)	! data type  

  INTEGER	:: offset, offsetpts, offsetcoords
  INTEGER	:: offset_end, offsetpts_end, offsetcoords_end  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: offset_master
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: offset_worker

  INTEGER :: nsou_sub
  INTEGER :: avesou, extrasou		! used to determine the number of sources to send to each worker

  REAL(KIND=iiii10), DIMENSION(8)	::	params1
  INTEGER, DIMENSION(11)		::	params2  

  REAL(KIND=iiii10), DIMENSION(:,:), ALLOCATABLE  :: 	sou_worker    ! sources to send to workers for calculating a "fraction" of the modrays algorithm.
  REAL(KIND=iiii10), DIMENSION(:,:), ALLOCATABLE  :: 	sou_subset    ! used by workers for calculating a "fraction" of the modrays algorithm.
  REAL(KIND=iiii10), DIMENSION(:,:), ALLOCATABLE  :: 	rec_worker    ! (all) receivers used in workers for calculating a "fraction" of the modrays algorithm.
  REAL(KIND=iiii10), DIMENSION(:,:), ALLOCATABLE  :: 	rec_subset    ! used by workers for calculating a "fraction" of the modrays algorithm. "Subset version".

  INTEGER, DIMENSION(:), ALLOCATABLE 		::	nsou_master  
  INTEGER, DIMENSION(:), ALLOCATABLE 		::	nsou_worker  
  INTEGER			 		::	nsou_subset  

  INTEGER, DIMENSION(:,:), ALLOCATABLE 		::	raystat_worker  
  INTEGER, DIMENSION(:,:), ALLOCATABLE 		::	raystat_subset  
  INTEGER, DIMENSION(:,:), ALLOCATABLE		::	switch
  
  REAL(KIND=iiii10), DIMENSION(:,:), ALLOCATABLE  :: 	VELI_worker    ! VELI matrix for workers
              
  INTEGER 					:: 	nru_master
  INTEGER, DIMENSION(:,:), ALLOCATABLE		::	npoints_master

  INTEGER					::	numraycoords_master
  REAL(KIND=iiii10), DIMENSION(:,:), ALLOCATABLE	::	raycoords_master 
  INTEGER					::	crazyray_master

  !new variables for base case test
  REAL(KIND=iiii10)				:: 	t_rays_test_base, t1_rays_test_base, t2_rays_test_base
  INTEGER 					:: 	npoints_diff1, npoints_diff2, npoints_base_sum1, npoints_base_sum2, npoints_sum1, npoints_sum2
  REAL(KIND=iiii10)				:: 	raycoords_diffx, raycoords_diffy, raycoords_base_sumx, raycoords_base_sumy
  REAL(KIND=iiii10)				:: 	raycoords_sumx, raycoords_sumy

END MODULE mksamples_prllvars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!                MAIN PROGRAM                 !!!!!!!!!!!!!! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM mksamples_fastrays

  USE mksamples_prllvars
  USE MPI
  USE fm2dssvars
  USE fm2dout
  USE mt19937
  IMPLICIT NONE  

  INTEGER :: status(MPI_STATUS_SIZE)
!*****************************************************************************************
! DECLARATION OF VARIABLES
!*****************************************************************************************
  CHARACTER (LEN=30)	:: sources,receivers,otimes,rayini,cdum
  CHARACTER (LEN=30)	:: sampleout,ountc,updout,rankc,runnumc
  CHARACTER (LEN=30)	:: sampfile,ncelfile,sigfile,misffile,aratfile
  CHARACTER (LEN=40)	:: updvtx,updvor,samplevtx,samplevor,string,rayupdout,last_tmp,del_tmp,sigma_tmp,arat_tmp
  LOGICAL		:: file_exists,file_open
!-----------------------------------------------------------------------------------------
! Markov chain parameters
! The program can be run in one go or as separate runs (in case there are restrictions in
! the job runtime). In the latter case, the state of the chain is saved to disk and the
! chain can be restarted from the last accepted model after changing the run number at
! line 13 of the input file mksamples.in.
!-----------------------------------------------------------------------------------------
  INTEGER 	:: nsample  			! Post burn-in
  INTEGER 	:: runnum			! Run number
! Each chain is run for 'nsample' steps in total. 
! The convergence of the algorithm is monitored with a number of indicators such as 
! acceptance ratios, and sampling efficiency is optimized by adjusting the variance of 
! the Gaussian proposal functions. 
!-----------------------------------------------------------------------------------------
! Some parameters for the PRIOR distribution
!-----------------------------------------------------------------------------------------
  ! Bounds in the number of cells (already declared in the fm2dssvars module)
  !INTEGER 		:: ncell_min
  !INTEGER 		:: ncell_max
  REAL(KIND=ii10)	:: mean			! Mean velocity of the uniform prior
  REAL(KIND=ii10) 	:: theta		! Lower and upper bounds of the velocity prior are [mean-theta , mean+theta]
  REAL(KIND=ii10) 	:: sigma_min,sigma_max 	! Lower and upper bounds for noise parameter (if noise is constant along all raypaths)
  REAL(KIND=ii10) 	:: aa_min,aa_max,bb_min,bb_max	! Parameters for gradient and y-intercept of noise if noise proportional to source-receiver distance or raypath length
  REAL(KIND=ii10) 	:: lambda_min,lambda_max	! Parameters for scaling factor of recorded noise
  ! Bounds of the 2D region (already declared in the fm2dssvars module)
  !REAL(KIND=ii10)	:: longmin, longmax
  !REAL(KIND=ii10)	:: latmin, latmax
!-----------------------------------------------------------------------------------------
! Standard deviations for PROPOSAL distributions
!-----------------------------------------------------------------------------------------
! These values have to be "tuned" so that the acceptance rates written in mpi.out (or to 
! the screen) are as close as possible to 44%. 
! This determines the efficiency of posterior sampling.
! If AR larger than 44%, increase the StDev for less acceptance.
! If AR smaller than 44%, decrease the StDev for more acceptance.
  REAL(KIND=ii10) 	:: pv    		! Proposal on velocity
  REAL(KIND=ii10) 	:: pv2  		! Proposal on velocity for DR proposal
  REAL(KIND=ii10) 	:: pd    		! Proposal on change in position
  REAL(KIND=ii10) 	:: pd2   		! Proposal on change in position for DR proposal
  REAL(KIND=ii10) 	:: ps,pa,pb,pl		! Proposal on change in data noise parameters
  REAL(KIND=ii10) 	:: sigmav 		! Proposal on velocity when Birth move
!-----------------------------------------------------------------------------------------
! Parameters for displaying results 
!-----------------------------------------------------------------------------------------
  INTEGER 	:: display    			! Display results at this interval (in mpi.out if parallelized) 
!-----------------------------------------------------------------------------------------
! Parameters for the Voronoi tesselation
!-----------------------------------------------------------------------------------------
  INTEGER 	:: npoint_max			! Maximum number of ray points
  INTEGER 	:: nray_max 			! Maximum number rays
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  INTEGER		thin
  INTEGER		node
  REAL(KIND=ii10)	tetaS,tetaR,phiS,phiR,alpha,angle
  REAL(KIND=ii10) 	lat,long
  DOUBLE PRECISION, EXTERNAL	:: GASDEV
  REAL(KIND=ii10) 	log,sqrt
  REAL(KIND=ii10) 	t1,t2,t3,t4, t1_rays, t2_rays, t_rays, t_uprays, t_uprays_test, t_rays_test, t1_rays_test, t2_rays_test
  REAL(KIND=ii10) 	lgsigma
  REAL(KIND=ii10) 	u,TT1,TT2,TT2a,TT2b,TT3,TT3t,TT3b,alpha1,alpha2,one,rannum,rannum2
  INTEGER 		w,as
  REAL(KIND=ii10) 	AR(3),ARB(3),ARM(3),ARV(3),ARD(3),AR1p(3),AR2p(3),AR1v(3),AR2v(3),ARS(3),ARSa(3),ARSb(3)
  REAL(KIND=ii10) 	like,like_prop,like_pprop,misfit
  INTEGER 		p,pp,cha,chb
  INTEGER 		nraysr,nrr,nr,sampletotal
  INTEGER 		ount,ounti,totacc,totacci,randcount
  INTEGER 		ncell,ncell_prop,i,j,s,r,sample,ind,a,a2,c
  REAL(KIND=ii10) 	p1,p2,prob
  INTEGER		v
  REAL(KIND=ii10) 	n,t
  REAL(KIND=ii10) 	point(2)
  REAL(KIND=ii5) :: 	ptlat,ptlong	! Added by EG to read input raypath file (binary output file from fm2dssRJ)
  INTEGER 		np,ii,jj,iii,jjj
  INTEGER		istart,iend,jstart,jend,k
  INTEGER		rayid,raylength,pti,ptf
  INTEGER		urps,upthin,ops,oilm,opsi,slmur,uglpd
  INTEGER 		birth,death,move,velch,sigmach,nodata
  LOGICAL           	accept,DRv,DRp
  REAL(KIND=ii10)	latp,longp,latminp,latmaxp,longminp,longmaxp
  ! Allocatable Arrays
  LOGICAL, DIMENSION(:), ALLOCATABLE 			:: 	cellogic,RV
  INTEGER , DIMENSION(:), ALLOCATABLE 			::	nnM
  INTEGER, DIMENSION(:), ALLOCATABLE 			::	lengthray
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE		:: 	dato
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		:: 	RayCel,RayCel_prop,RayCel_pprop
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		::	Voro,Voro_prop,Voro_pprop
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		::	firstVoro
  INTEGER, DIMENSION(:,:), ALLOCATABLE 			::	samplesStep
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		::	samplesVal
  INTEGER, DIMENSION(:), ALLOCATABLE 			::	samplesInd
  INTEGER, DIMENSION(:), ALLOCATABLE 			::	samplesNcel
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		::	samplesSigmas
  REAL(KIND=ii10), DIMENSION(:), ALLOCATABLE 		::	samplesMisfit
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		::	samplesAratios
  REAL(KIND=ii10), DIMENSION(:), ALLOCATABLE 		::	ttime,ttime_prop,ttime_pprop
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		::	sigma,sigma_prop
  REAL(KIND=ii10), DIMENSION(:), ALLOCATABLE		::	aa,bb,aa_prop,bb_prop,lambda,lambda_prop
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		:: 	VELI,VELU
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 		:: 	M,Mi
  REAL(KIND=ii5), DIMENSION(:,:), ALLOCATABLE 		::	raypoints
  REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE		::	lastacc
  INTEGER 						::	done,nlastacc,sw
  INTEGER 				:: nds,nss,dssval
  CHARACTER (LEN=30) 			:: dsfile
  INTEGER, DIMENSION(:), ALLOCATABLE 	:: dsstat,dsstato



INTERFACE
SUBROUTINE updaterays(RayCel_prop,RayCel,Voro_prop,ncell,ncell_prop,cellogic,ind,pn,sn, t_uprays, t_uprays_test)

		USE mksamples_prllvars
		USE MPI
		USE fm2dssvars
		USE fm2dout
		
		IMPLICIT NONE
		
		INTEGER			::	i,ii,iii,j,jj,jjj,node,pp,pu,nr,nc,nrr,ind,rayid,raylength,pti,ptf,pn,sn
		INTEGER			::	nsin,nrin,ncell,ncell_prop
		REAL(KIND=ii10)		::	p1,p2,point(2), time1_upr, time2_upr, time1_upr_test, time2_upr_test
		LOGICAL 		::	cellogic(ncell_max)
		REAL(KIND=ii10)	  	::	tetaS,tetaR,phiS,phiR,alpha,angle,RayCel(nrays,ncell)

	        INTEGER :: status(MPI_STATUS_SIZE)
		
		REAL(KIND=ii5), DIMENSION(:,:), ALLOCATABLE 	::	raypoints
		REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 	::	Mu
		REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 	::	RAYUP
		REAL(KIND=ii10), DIMENSION(1:,1:) 		:: 	Voro_prop
		REAL(KIND=ii10), DIMENSION(1:,1:)		::	RayCel_prop
		CHARACTER (LEN=30) 				::	pnc,snc

		REAL(KIND=ii10), INTENT(out) :: t_uprays, t_uprays_test
		
END SUBROUTINE updaterays
END INTERFACE

!****************************************************************
!////////////////////////////////////////////////////////////////
!****************************************************************

! This is just to test the code using testrays.exe
! The WRITE and READ lines can be commented out in the final version of the code
test=0
!WRITE(*,*)'Test? 1=yes,0=no'
!READ(*,*)test

WRITE(*,*)'Hello from reverse-jump MCMC tomography!'

!-------------------------------------------------------------
! Read input file rj_tomo.in and get values of variables
!-------------------------------------------------------------

OPEN(UNIT=20,FILE='mksamples.in',STATUS='old')
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)sources
READ(20,1)receivers
READ(20,1)otimes
READ(20,1)rayini
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)runnum
READ(20,*)nsample
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)ncell_min,ncell_max
READ(20,*)mean
READ(20,*)theta
READ(20,*)sigdep
READ(20,*)uglpd
READ(20,*)sigma_min,sigma_max
READ(20,*)aa_min,aa_max,bb_min,bb_max
READ(20,*)lambda_min,lambda_max
READ(20,*)latmax,latmin
READ(20,*)longmin,longmax
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)nds
READ(20,1)dsfile
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)pv
READ(20,*)pv2
READ(20,*)pd
READ(20,*)pd2
READ(20,*)ps
READ(20,*)pa,pb
READ(20,*)pl
READ(20,*)sigmav
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)display
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)npoint_max
READ(20,*)nray_max
READ(20,*)nt_max
READ(20,*)nv_max
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)sampfile
READ(20,1)ncelfile
READ(20,1)sigfile
READ(20,1)misffile
READ(20,1)aratfile
READ(20,*)ops
READ(20,*)opsi
READ(20,1)sampleout
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)urps
READ(20,*)upthin
READ(20,*)uar
READ(20,*)slmur
READ(20,1)updout
READ(20,*)wurf
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)nvtr,nvpr
READ(20,*)gdt,gdp
READ(20,*)sgref
READ(20,*)dl,erg
READ(20,*)er
READ(20,*)fms
READ(20,*)nbs
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)sse
READ(20,1)errout
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)nodata
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)rank
READ(20,*)nbproc
READ(20,*)TEST_MODRAYS
1   FORMAT(a30)
CLOSE(20)

sse=0	! Do not save samples that cause errors

nmax=3*nt_max+ncell_max

! Check input values are ok
IF(display.LE.0) THEN
	IF(rank.EQ.0) WRITE(*,*)'Display interval cannot be less than or equal to zero!!'
	STOP
ELSE IF(opsi.LT.0) THEN
	IF(rank.EQ.0) WRITE(*,*)'Output sample interval cannot be less than zero!!'
	STOP	
ELSE IF(upthin.LE.0.AND.urps.EQ.1) THEN
	IF(rank.EQ.0) WRITE(*,*)'Ray update interval cannot be less than or equal to zero!!'
	STOP
ELSE IF(runnum.LE.0) THEN
	IF(rank.EQ.0) WRITE(*,*)'Run number cannot be less than or equal to zero!!'
	STOP
ELSE IF(nds.LE.0) THEN
	IF(rank.EQ.0) WRITE(*,*)'Dataset number can only be equal to or greater than one!!'
	STOP
END IF
	
	
!-------------------------------------------------------------
! Allocate allocatable arrays using the input values
!-------------------------------------------------------------
ALLOCATE(vertices(3,nt_max))
ALLOCATE(neighbour(3,nt_max))
ALLOCATE(nnn(ncell_max+1))
ALLOCATE(nnlist(nmax))
ALLOCATE(ntwork(nmax))
ALLOCATE(worki1(nv_max))
ALLOCATE(worki2(nv_max))
ALLOCATE(worki3(nv_max))
ALLOCATE(lengthray(nray_max))
ALLOCATE(Voro(3,ncell_max))
ALLOCATE(Voro_prop(3,ncell_max))
ALLOCATE(Voro_pprop(3,ncell_max))
ALLOCATE(firstVoro(3,ncell_max))
ALLOCATE(dato(nray_max,6))
ALLOCATE(dsstato(nray_max))
ALLOCATE(cellogic(ncell_max))
ALLOCATE(ldummy(ncell_max))
ALLOCATE(RV(nray_max))
vertices=0
neighbour=0
nnn=0
nnlist=0
ntwork=0
worki1=0
worki2=0
worki3=0
lengthray=0
Voro=0
Voro_prop=0
Voro_pprop=0
firstVoro=0
dato=0
dsstato=0
ldummy=0
RV=0

IF (runnum.EQ.1) THEN
	ALLOCATE(samplesStep(nsample+1,2))
	ALLOCATE(samplesVal(nsample+1,3))
	ALLOCATE(samplesInd(nsample+1))
	ALLOCATE(samplesNcel(nsample+1))
	IF (sigdep.EQ.1.OR.sigdep.EQ.2) ALLOCATE(samplesSigmas(nsample+1,nds*2))
	IF (sigdep.EQ.0.OR.sigdep.EQ.3) ALLOCATE(samplesSigmas(nsample+1,nds))
	ALLOCATE(samplesMisfit(nsample+1))
	ALLOCATE(samplesAratios(nsample+1,6))
ELSE IF (runnum.GT.1) THEN
	ALLOCATE(samplesStep(nsample,2))
	ALLOCATE(samplesVal(nsample,3))
	ALLOCATE(samplesInd(nsample))
	ALLOCATE(samplesNcel(nsample))
	IF (sigdep.EQ.1.OR.sigdep.EQ.2) ALLOCATE(samplesSigmas(nsample,nds*2))
	IF (sigdep.EQ.0.OR.sigdep.EQ.3) ALLOCATE(samplesSigmas(nsample,nds))
	ALLOCATE(samplesMisfit(nsample))
	ALLOCATE(samplesAratios(nsample,6))
END IF
samplesStep = 0
samplesVal = 0
samplesInd = 0
samplesNcel = 0
samplesSigmas = 0
samplesMisfit = 0
samplesAratios = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Original code mofified by Helmut Wahanik, Research Scientist, Schlumberger Research.
!
! hybrid embarrasingly-parallel and mpi code
! 
! Original code mofified by Helmut Wahanik (HW).  Tays-routines MPI parallelization
! Here we build 1 Markov chain instead of an mpi parallelization of nbproc processors
! The mpi parallelization will be held at the ray algorithms stage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1. Outer embarrasingly-parallel loop

! In this version the rank would not come from an MPI processor.  Rather, it will come from the submission 
! number to the cluster, such as Gondor (SGR) for example.  This might look as a MPI process distribution, 
! but it is simply an embarrisingly parallel process

! Rank number (character)

rank = 4     ! used for testing seeral single jobs submission
! nbproc = 1   ! number of chains submitted in the embarrasingly parallel loop.									     		

ra=rank ! comment/uncomment for ra assign.
WRITE(rankc,*)rank+1    
rankc = ADJUSTL(rankc)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2. Internal mpi parallelization

! Start parallelization of the code. 

 CALL MPI_INIT(ierror)				   
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierror) ! Originally we had that nbproc = numtasks.  rank = taskid.  
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierror)   


numworkers = numtasks-1       ! The first processor is completely dedicated to the Markov chain creation.  
	       		      ! The remaining processors are available for running subsets of the ray-routines processing.
wcount = 1		      ! wcount is the number of samples and is used as flag for the mpi workers.	
run_worker = 0		      ! run_worker the flag for the while loop of the mpi workers	

!Nbproc corresponds now to embarrasingly parallel processes, and rank is the index of the independent job submitted to the cluster
      

MASTER = 0 
FROM_MASTER = 1 
FROM_WORKER = 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 CALL cpu_time(t1)	! TIC - Start counting time 
 t_rays = 0 		! Counter for modray routine and to calculate improvement from parallelization
 t_uprays = 0 		! Counter for modrays routine inside of update rays

 t_rays_test = 0 
 t_uprays_test = 0 

!! ---------------------------------------------- master section ----------------------------------------------------------!!

IF (taskid.EQ.MASTER) THEN  ! this conditional goes up to the end of mksamples loop, and after the outputs have been exported

WRITE(*,*) 'I am the master of', numtasks, 'processes'

WRITE(*,*) 'Number of exclusively worker tasks is', numworkers
  
IF (runnum.EQ.1) THEN
	last_tmp = 'last.tmp' // trim(rankc)
	sigma_tmp = 'sigma.tmp' // trim(rankc)
	arat_tmp = 'arats.tmp' // trim(rankc)
  
	INQUIRE(FILE=last_tmp,EXIST=file_exists)
	IF (file_exists) THEN
		del_tmp='rm ' // last_tmp
		CALL SYSTEM( del_tmp )
	END IF

	INQUIRE(FILE=sigma_tmp,EXIST=file_exists)
	IF (file_exists) THEN
		del_tmp='rm ' // sigma_tmp
		CALL SYSTEM( del_tmp )
	END IF

	INQUIRE(FILE=arat_tmp,EXIST=file_exists)
	IF (file_exists) THEN
		del_tmp='rm ' // arat_tmp
		CALL SYSTEM( del_tmp )
	END IF
ELSE IF (runnum.GT.1) THEN
	last_tmp = 'last.tmp' // trim(rankc)
	sigma_tmp = 'sigma.tmp' // trim(rankc)
	arat_tmp = 'arats.tmp' // trim(rankc)
END IF
  
! Allocate arrays for last accepted model   
ALLOCATE(lastacc(3,ncell_max))
lastacc=0
nlastacc=0
  
  
CALL date_and_time(datevals(1),datevals(2),datevals(3),date_time)
  clock=real(date_time(7))+real(date_time(8))/1000
!ranseed=NINT(clock*(rank+1))
ranseed=rank+1
CALL init_genrand(ranseed)
WRITE(*,*)'Ranseed=',ranseed,'on rank',rank+1

DRp=.true.
DRv=.true.
RV=.false.

IF (runnum.EQ.1) THEN
	totacc=0
	totacci=0
	ounti=0
	AR=0
	ARB=0
	ARD=0
	ARM=0
	AR1p=0
	AR2p=0
	ARV=0
	AR1v=0
	AR2v=0
	ARS=0
	ARSa=0
	ARSb=0
	randcount=0
ELSE IF (runnum.GT.1) THEN
	OPEN(UNIT=84,FILE=arat_tmp,FORM='unformatted',STATUS='old')
	READ(84)totacc
	totacci=totacc
	READ(84)ounti
	READ(84)AR(1),AR(2),AR(3)
	READ(84)ARB(1),ARB(2),ARB(3)
	READ(84)ARD(1),ARD(2),ARD(3)
	READ(84)ARM(1),ARM(2),ARM(3)
	READ(84)AR1p(1),AR1p(2),AR1p(3)
	READ(84)AR2p(1),AR2p(2),AR2p(3)
	READ(84)ARV(1),ARV(2),ARV(3)
	READ(84)AR1v(1),AR1v(2),AR1v(3)
	READ(84)AR2v(1),AR2v(2),AR2v(3)
	READ(84)ARS(1),ARS(2),ARS(3)
	IF (sigdep.EQ.1.OR.sigdep.EQ.2) READ(84)ARSa(1),ARSa(2),ARSa(3)
	IF (sigdep.EQ.1.OR.sigdep.EQ.2) READ(84)ARSb(1),ARSb(2),ARSb(3)
	READ(84)randcount
	CLOSE(84)
	! Call grnd() as many times as it was called in the previous runs of the code
	! (this ensures that the random sequence restarts from the point it was paused)
	DO i=1,randcount
		u=grnd()
	END DO
END IF

sampletotal = nsample+ounti
prob=0.0
one=1.0

lat = (latmax - latmin)/REAL(nvt-1,KIND=ii10)
long = (longmax - longmin)/REAL(nvp-1,KIND=ii10)
latr = (latmax - latmin)/REAL(nvtr-1,KIND=ii10)
longr = (longmax - longmin)/REAL(nvpr-1,KIND=ii10)


!**************************************************************
!                        READ INPUT FILES
!**************************************************************

! Read the stations coordinates file
OPEN(UNIT=1,FILE=receivers,STATUS='old')
READ(1,*) nrec
! rec(:,1) is latitude and rec(:,2) is longitude
ALLOCATE(rec(nrec,2))
DO i=1,nrec
	READ(1,*)rec(i,1),rec(i,2)
	IF(rec(i,1).GT.latmax.OR.rec(i,1).LT.latmin.OR.rec(i,2).GT.longmax.OR.rec(i,2).LT.longmin) THEN
		WRITE(*,*)'Receiver',i,'is located outside model boundaries!!'
		WRITE(*,*)'Terminating program!!'
		CALL MPI_FINALIZE(ierror) !! uncomment for former MPI parallelization
		STOP
	END IF
END DO
CLOSE(1)
  
! Read the sources coordinates file(in ambient noise, source=station)
OPEN(UNIT=2,FILE=sources,STATUS='old')
READ(2,*) nsou
! sou(:,1) is latitude and sou(:,2) is longitude
ALLOCATE(sou(nsou,2))
DO i=1,nsou
	READ(2,*)sou(i,1),sou(i,2)
	IF(sou(i,1).GT.latmax.OR.sou(i,1).LT.latmin.OR.sou(i,2).GT.longmax.OR.sou(i,2).LT.longmin) THEN
		WRITE(*,*)'Source',i,'is located outside model boundaries!!'
		WRITE(*,*)'Terminating program!!'
		CALL MPI_FINALIZE(ierror) !! uncomment for former MPI parallelization
		STOP
	END IF
END DO
CLOSE(2)
  
! Check that the total number of source-station pairs does not exceed the limit given in the input file
IF(nsou*nrec.GT.nray_max) THEN
	IF(rank.EQ.0) WRITE(*,*)'The number of source-station pairs exceeds the maximum number of raypaths nray_max set in the input file!!'
	CALL MPI_FINALIZE(ierror) !! uncomment for former MPI parallelization
	STOP
END IF

!Get the DATA
OPEN(UNIT=3,FILE=otimes,STATUS='old')			! File containing source-receiver associations and traveltimes
IF (nds.GT.1) OPEN(UNIT=6,FILE=dsfile,STATUS='old')	! File containing source-receiver associations and dataset info

! We test all possible couples of stations 

ALLOCATE(raystat(nsou*nrec,2)) ! raystat contains the switch in column 1 and the valid ray number in column 2

p = 0
nrr = 0 ! "loops" over rays
DO s=1,nsou
   DO r=1,nrec
        nrr = nrr + 1 
	READ(3,*)v,t,n   ! read validity of the data (switch), traveltime , noise (error in traveltime)
	IF (nds.GT.1) READ(6,*)dssval
	IF (nds.EQ.1) dssval=1
 	
	IF (v.ne.0) THEN
		RV(nrr)=.true.
		p = p+1 
                dato(p,1:2)=sou(s,1:2)
                dato(p,3:4)=rec(r,1:2)
                dato(p,5)=t
		dato(p,6)=n
		raystat(nrr,1)=1
		raystat(nrr,2)=p
		dsstato(p)=dssval
	ELSE 
		RV(nrr)=.false.
		raystat(nrr,1)=0
		raystat(nrr,2)=0
	END IF
   END DO
END DO

CLOSE(3)! close the file
IF (nds.GT.1) CLOSE(6)! close the file

nrays = p

ALLOCATE(dat(nrays,6))
dat=dato(1:nrays,:)
DEALLOCATE(dato)
ALLOCATE(dsstat(nrays))
dsstat=dsstato(1:nrays)
DEALLOCATE(dsstato)

! Allocate arrays
ALLOCATE(RayCel(nrays,ncell_max))
ALLOCATE(RayCel_prop(nrays,ncell_max))
ALLOCATE(RayCel_pprop(nrays,ncell_max))
ALLOCATE(ttime(nrays))
ALLOCATE(ttime_prop(nrays))
ALLOCATE(ttime_pprop(nrays))
ALLOCATE(lengthrayv(nrays))
ALLOCATE(dis(nrays))
RayCel=0
RayCel_prop=0
RayCel_pprop=0
ttime=0
ttime_prop=0
ttime_pprop=0
lengthrayv=0
dis=0

! Allocate arrays for source-receiver distances
IF (sigdep.EQ.1.OR.sigdep.EQ.2) ALLOCATE(srdist(nrays))
IF (sigdep.EQ.2) ALLOCATE(srdist_prop(nrays))
IF (sigdep.EQ.1) THEN
	! Calculate source-receiver distances as angle in degrees
        srdist=180/pii*2*asin(sqrt((sin(((dat(:,1)*pii/180)-(dat(:,3)*pii/180))/2))**2+(cos((dat(:,3)*pii/180)))*(cos((dat(:,1)*pii/180)))*(sin(((dat(:,2)*pii/180)-(dat(:,4)*pii/180))/2))**2))
ELSE IF (sigdep.EQ.2) THEN
	srdist=0
	srdist_prop=0
END IF

! Allocate arrays for data noise (sigma)
IF(sigdep.EQ.0) THEN
	ALLOCATE(sigma(1,nds))
	ALLOCATE(sigma_prop(1,nds))
ELSE IF(sigdep.EQ.1.OR.sigdep.EQ.2.OR.sigdep.EQ.3.OR.sigdep.EQ.4) THEN
	IF (sigdep.EQ.1.OR.sigdep.EQ.2) ALLOCATE(aa(nds),aa_prop(nds),bb(nds),bb_prop(nds))
	IF (sigdep.EQ.3) ALLOCATE(lambda(nds),lambda_prop(nds))
	ALLOCATE(sigma(nrays,1))
	ALLOCATE(sigma_prop(nrays,1))
ELSE
	WRITE(*,*)'Wrong input at line 22 of mksamples.in - terminating program!!'
	CALL MPI_FINALIZE(ierror) !! uncomment for former MPI parallelization
	STOP
END IF
sigma=0
sigma_prop=0


! Display information
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)'------------------------------------------------------------------------'
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)'Number of valid rays:',nrays

 
! Compute average observed travel time 
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)'Average obs time:',SUM(dat(:,5))/nrays,'seconds'
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)'------------------------------------------------------------------------'

!****************************************************************************
!    Draw the first model randomly from the prior distribution
!****************************************************************************

! If this is not the first run of the chain...
IF (runnum.NE.1) THEN
	! ...read last accepted model of previous run from file
	OPEN(UNIT=84,FILE=last_tmp,FORM='unformatted',STATUS='old')
	READ(84)ncell
	DO c=1,ncell
		READ(84)Voro(2,c),Voro(1,c),Voro(3,c)
	END DO
	CLOSE(84)
END IF

ptsok=0
! If raypaths are not updated at each iteration...
IF(urps.NE.2) THEN
	! If this is the first run, initialize a random velocity model
111	IF (runnum.EQ.1) THEN
	!* Initial number of cells
		ncell=int(ncell_min+grnd()*(ncell_max-ncell_min))
		randcount=randcount+1
		!* We place the cells randomly and give them random velocity values
		DO i=1,ncell
		 	! Initial location of cells (Voro(1:2,i))
		   	Voro(1:2,i) = [longmin+grnd()*(longmax-longmin),latmin+grnd()*(latmax-latmin)]
			! Initial velocity value for Voronoi cells (Voro(3,i))
		   	Voro(3,i) = mean - theta +2*grnd()*theta
		   	randcount=randcount+3
		END DO
	END IF
	!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
	WRITE(*,*)'Initial number of cells:',ncell,'on rank',rank+1
	!CALL cpu_time(t3)
	!* Compute the Voronoi grid
	CALL delaun(DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,& 	! outputs are neighbour, vertices, nt (number of triangles)
	                      worki1,worki2,worki3,eps,nv_max,&
	                      0,ldummy,0,0,0,ptsok)
	IF (ptsok.EQ.1) THEN
		WRITE(*,*)'All',ncell,'points in a line for initial model on processor',rank+1
		WRITE(*,*)'Proposing new initial model...'
		GOTO 111
	END IF
	CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
	              neighbour,nnn,nnlist,ntwork)
	!CALL cpu_time(t4)
	!write(*,*)'Triangulation time:',t4-t3,'s'

! ...else if raypaths are updated at each iteration...
ELSE IF(urps.EQ.2) THEN
112	crazyray=1
	DO WHILE(crazyray.GT.0)
		IF (runnum.EQ.1) THEN
			! If this is the first run, initialize a random velocity model that does not cause any "crazy" raypaths
			!* Initial number of cells
			ncell=int(ncell_min+grnd()*(ncell_max-ncell_min))
			randcount=randcount+1
			!* We place the cells randomly and give them random velocity values
			DO i=1,ncell
				! Initial location of cells (Voro(1:2,i))
				Voro(1:2,i) = [longmin+grnd()*(longmax-longmin),latmin+grnd()*(latmax-latmin)]
				! Initial velocity value for Voronoi cells (Voro(3,i))
				Voro(3,i) = mean - theta +2*grnd()*theta 
				randcount=randcount+3
			END DO
		END IF
		!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
		WRITE(*,*)'Initial number of cells:',ncell,'on rank',rank+1
		!CALL cpu_time(t3)
		!* Compute the Voronoi grid
		CALL delaun(DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,& 	! outputs are neighbour, vertices, nt (number of triangles)
                      			worki1,worki2,worki3,eps,nv_max,&
                      			0,ldummy,0,0,0,ptsok)
		IF (ptsok.EQ.1) THEN
			WRITE(*,*)'All',ncell,'points in a line for initial model on processor',rank+1
			WRITE(*,*)'Proposing new initial model...'
			GOTO 112
		END IF
		CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
              			neighbour,nnn,nnlist,ntwork)
		!CALL cpu_time(t4)
		!write(*,*)'Triangulation time:',t4-t3,'s'
		
		IF(ALLOCATED(VELI)) DEALLOCATE(VELI)
		ALLOCATE(VELI(nvpr+2,nvtr+2))	! This is the matrix containing the initial model velocities that will be used to calculate the initial raypaths	
	
		node=1
					
		iii=0
		DO i=0,nvpr+1
			jjj=0
			iii=iii+1
			DO j=0,nvtr+1
				jjj=jjj+1					
				ii=i
				jj=j
				IF (i.EQ.0)   	ii=1
				IF (j.EQ.0)   	jj=1
				IF (i.EQ.nvpr+1) 	ii=nvpr
				IF (j.EQ.nvtr+1) 	jj=nvtr
				
				p1=longmin+(ii-1)*longr
				p2=latmax-(jj-1)*latr
			
				node=1
			
				CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro(1:2,:)),nnn,nnlist,walk)
			
				VELI(iii,jjj) = (Voro(3,node))
			
			END DO
		END DO
		
		nru=sum(raystat(:,1),1)
		
		!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                WRITE(*,*) 'Number of valid rays', nru
	
		IF(nru.NE.nrays) THEN
			IF(rank.EQ.0) WRITE(*,*)'There is a problem with the number of active rays!!'
			CALL MPI_FINALIZE(ierror)
			STOP
		END IF
		
		IF(rank.EQ.0) WRITE(*,*)			
		IF(rank.EQ.0) WRITE(*,*)'Calling MODRAYS to model raypaths for the initial model'
		IF(rank.EQ.0) WRITE(*,*)



		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	        ! The section below corresponds to the traditional case for the modrays algorithm. Use for algorithm output comparison.

		IF(TEST_MODRAYS.EQ.1) THEN 

			! Time counter

			t1_rays_test_base = 0
	  	        t2_rays_test_base = 0
			CALL cpu_time(t1_rays_test_base) ! TOC - Start counting time

			IF(ALLOCATED(npoints)) DEALLOCATE(npoints)
			IF(ALLOCATED(raycoords)) DEALLOCATE(raycoords)
			IF(ALLOCATED(raypoints)) DEALLOCATE(raypoints)
			IF(ALLOCATED(Mi)) DEALLOCATE(Mi)

		        !!! We want to send arrays independently

			CALL modrays(nsou,sou(:,1),sou(:,2),&
				nrec,rec(:,1),rec(:,2),&
				raystat,-1,&
				nvtr,nvpr,latmax,longmin,latr,longr,VELI,&
				gdt,gdp,&
				sgref,&
				dl,erg,&
				er,&
				fms,&
				nbs,&
				rank+1)
		
		        IF(ALLOCATED(npoints_base)) DEALLOCATE(npoints_base)
			IF(ALLOCATED(raycoords_base)) DEALLOCATE(raycoords_base)
		        nru_base = nru
			crazyray_base = crazyray

			!WRITE(*,*) 'Done with test!', SIZE(npoints), SIZE(raycoords)

		        ALLOCATE( npoints_base( (SIZE(npoints))/2,2 ) )
		        ALLOCATE( raycoords_base( (SIZE(raycoords)/2),2 ) )
		        npoints_base = npoints
		        raycoords_base = raycoords

			CALL cpu_time(t2_rays_test_base) ! TOC - Stop counting time 
			! Calculate cumulative time spent in modrays routine when parallelized
	 		t_rays_test_base =  t_rays_test_base  +  (t2_rays_test_base-t1_rays_test_base)
		END IF

	        ! The section above corresponds to the traditional case for the modrays algorithm.   Use for algorithm output benchmark.
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Inner parallelization starts here.  
		! We start by sending the source sets to the workers for the modrays routine for calculation.  We will send the partial  
		! information required by the modrays routines to every worker.      	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF(ALLOCATED(npoints)) DEALLOCATE(npoints)
		IF(ALLOCATED(raycoords)) DEALLOCATE(raycoords)
		IF(ALLOCATED(raypoints)) DEALLOCATE(raypoints)
		IF(ALLOCATED(Mi)) DEALLOCATE(Mi)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Time counter

		t1_rays_test = 0
 	        t2_rays_test = 0
		CALL cpu_time(t1_rays_test) ! TOC - Start counting time

		!t1_rays = 0
  	        !t2_rays = 0
		!CALL cpu_time(t1_rays) ! TOC - Start counting time
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!WRITE(*,*) 'Start parallelization 1'

		crazyray = 0 
		IF(ALLOCATED(offset_master)) DEALLOCATE(offset_master)
		IF(ALLOCATED(nsou_master)) DEALLOCATE(nsou_master)

		ALLOCATE(offset_master(numworkers,2))
		ALLOCATE(nsou_master(numworkers))

		! (We store data in several matrices and send to the worker tasks).  

		!WRITE(*,*) 'Number of sources', nsou
		avesou = nsou/numworkers  ! Number of sources in each subset, i.e. per worker excluding remainder.   e.g. 11/3 = 3. remainder 2.  
		extrasou = MOD(nsou,numworkers)		

                ! The offsets are the indices used for transiting along the arrays
		! In the fisrt loop we calculate the offsets that will be used by the workers.   
		! In the fisrt loop we calculate the nsources that will be send to each worker.   
  
                offset = 0        

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               DO dest=1,numworkers  
		    ! Here we need to count the number of rays corresponding to the sources sent to the workers.  Number of sources sent is nsou_sent.             
                    IF(dest.LE.extrasou) THEN     														
		        nsou_sub = avesou+1  	! nsou_subset are the number of sources sent to the worker for calculating the modrays algorithm.       		
                    ELSE 
                        nsou_sub = avesou     			
	            END IF
                    offset_end = offset + nsou_sub   ! This defines the interval we are sending.  Ex:  offset = 0, nsou_subset = 3, offset_end = 3.                      
		    nsou_master(dest) = nsou_sub
		    offset_master(dest, 1) = offset
 		    offset_master(dest, 2) = offset_end
		    ! ATTENTION offsetpts and offsetcoords are calculated and looped in the second master processor loop, but after receiving the npoints and raycoords results from the worker.
                    offset = offset_end	
  	        END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! TAG FROM_MASTER = 1              
                ! TO: DEST 
		! Send to worker ===> (dest,1)

                !Brodcast

		!nvx,nvz,goxd,gozd,dvxd,dvzd,velv
		params1(1) = nvtr
		params1(2) = nvpr
		params1(3) = latmax
		params1(4) = longmin
		params1(5) = latr
		params1(6) = longr
		params1(7) = er ! Earth radius
		params1(8) = nbs

                wrgf = -1 	                      
		params2(1) = nsou
		params2(2) = nrec
		params2(3) = wrgf          	! wrgf = write geometries to file?
		params2(4) = gdt		! Grid dicing
		params2(5) = gdp		! Grid dicing
		params2(6) = sgref		! Apply source grid refinement? (0=no,1=yes)
		params2(7) = dl			! Dicing level and extent of refined grid
		params2(8) = erg
		params2(9) = fms		! use first-order(0) or mixed-order(1) scheme
		params2(10)= rank+1		! Processor rank (in this case for embarrasingly parallel processing).
		params2(11)= wcount		! Sample number! 


  	        counter_1 = 2*nsou*nrec               
                counter_2 = (nvpr+2)*(nvtr+2)
		!CALL MPI_BCAST(params1, 8, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 		! Brodcast parameters 1 to all processes 

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               DO dest=1,numworkers  
		    CALL MPI_SEND(0, 1, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)			 	! Send run_w flag for running worker
		    CALL MPI_SEND(params1, 8, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror) 		! Brodcast parameters 1 to all processes 
    	            CALL MPI_SEND(params2, 11, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror) 			! Brodcast parameters 2 to all processes 
	            CALL MPI_SEND(offset_master, 2*numworkers, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)	
                    CALL MPI_SEND(nsou_master, numworkers, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)
		    CALL MPI_SEND(rec, nrec*2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)
		    CALL MPI_SEND(sou, nsou*2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)		   		   
		    CALL MPI_SEND(raystat, counter_1, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)
		    CALL MPI_SEND(VELI, counter_2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)
               END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
                !WRITE(*,*) 'All variables and arrays have been broadcasted to the workers'
		!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

	        ! TAG = FROM_WORKER
                ! FROM: Worker
                ! Receive in master from worker ===> (from_worker,2)		

                offsetpts 	= 0  
                offsetcoords 	= 0  

		! Now we loop over all workers receiving the information for building the complete arrays.            

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	        DO  i=1,numworkers
		    ! ATTENTION offsetpts and offsetcoords are calculated and looped here, but are updated receiving the nru, npoints and raycoords results from the worker.
                    IF(ALLOCATED(offset_master)) DEALLOCATE(offset_master)
                    IF(ALLOCATED(npoints_master)) DEALLOCATE(npoints_master)
  		    ALLOCATE( offset_master(numworkers,2) )  		    
		    from_worker = i		
 
		    !TAG pair:  (i,2)		
                    CALL MPI_RECV(offset_master, 2*numworkers, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)			    		   
                    CALL MPI_RECV(nru_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)	                    
		    ALLOCATE( npoints_master( nru_master , 2 ) )			! array of length nru_subset, filled w. number of points that each ray contains and ray ID corresponding to the subset chosen.
		    CALL MPI_RECV(npoints_master, 2*nru_master, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)  

		    !numraycoords_master = sum( npoints_master( 1:nru_master , 1 ) )	
    	            CALL MPI_RECV(numraycoords_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)  
                  
                    IF(ALLOCATED(raycoords_master)) DEALLOCATE(raycoords_master)
                    ALLOCATE( raycoords_master( numraycoords_master, 2  )  )		! allocate array of ray coordinates for the sources subset.                   
		    CALL MPI_RECV(raycoords_master, numraycoords_master*2, MPI_DOUBLE_PRECISION, from_worker, 2, MPI_COMM_WORLD, status, ierror)
		    CALL MPI_RECV(crazyray_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)

		    offsetpts_end = offsetpts + nru_master
                    offsetcoords_end = offsetcoords + numraycoords_master

		    ! Next we fill in the raycoords array

		    IF (i.EQ.1) THEN		                             		    
                       ALLOCATE( npoints(nru_master, 2 ) )			! npoints is allocated here for the first time
		       npoints = npoints_master					! npoints is filled here for the first time
                       ALLOCATE( raycoords( numraycoords_master, 2  )  )      	! allocate array of ray coordinates for the sources subset for the first time
                       raycoords = raycoords_master				! raycoords is filled here for the first time				
                       IF(ALLOCATED(npoints_prev))DEALLOCATE( npoints_prev )  
                       IF(ALLOCATED(raycoords_prev))DEALLOCATE( raycoords_prev )  
		       ALLOCATE( npoints_prev( nru_master , 2 )  )
		       ALLOCATE( raycoords_prev( numraycoords_master , 2 ) )
		       npoints_prev   = npoints					  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.
		       raycoords_prev = raycoords				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.

		    ELSE IF (i.GT.1) THEN
                       DEALLOCATE( raycoords ) 
                       DEALLOCATE( npoints ) 
                       ALLOCATE( npoints( offsetpts_end, 2  )  )  	 	   ! allocate array of ray coordinates for the sources subset.		
                       ALLOCATE( raycoords( offsetcoords_end, 2  )  )      	   ! allocate array of ray coordinates for the sources subset.		
                       npoints(1:offsetpts, 1:2) = npoints_prev			   ! npoints is partially filled here with the new npoints subset that has been received recently.	                          
		       npoints(offsetpts+1: offsetpts_end, 1:2) = npoints_master   ! npoints is partially filled here with the new npoints subset that has been received recently.	                          
 		       raycoords(1:offsetcoords, 1:2 ) = raycoords_prev
                       raycoords(offsetcoords+1: offsetcoords_end, 1:2) = raycoords_master ! raycoords is partially filled here with the new raycoords subset that has been received recently.	                          

                       IF(ALLOCATED(npoints_prev))DEALLOCATE( npoints_prev )  
                       IF(ALLOCATED(raycoords_prev))DEALLOCATE( raycoords_prev )
	    	       ALLOCATE( npoints_prev( offsetpts_end , 2 )  )
		       ALLOCATE( raycoords_prev( offsetcoords_end , 2 ) )
		       raycoords_prev = raycoords				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.
		       npoints_prev   = npoints				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.

		    END IF 

                    offsetpts    = offsetpts_end 
                    offsetcoords = offsetcoords_end  !OK ALL REVISED UP TO THIS POINT.
     		    numraycoords = sum( npoints( : , 1 ) )                    ! size of raycoords_subset found by adding the first column of the corresponding subset of npoints = number of total ray coordinates in the subset.			
                    crazyray     = crazyray + crazyray_master
  		    !WRITE(*,*) 'crazyray is', crazyray

		    ! Printing
		    !WRITE(*,*)'offset after 26 equals', offset_master(1,1), 'in processor = ',i
		    !WRITE(*,*)'offset after 26 equals', offset_master(3,1), 'in processor = ',i
		    !WRITE(*,*)'offset after 26 equals', offset_master(5,1), 'in processor = ',i
		    !WRITE(*,*)'offset after 26 equals', offset_master(7,1), 'in processor = ',i	
		    !WRITE(*,*)'nru_master is ', nru_master, 'in processor = ',i
		    !WRITE(*,*)'npoints_master is ', npoints_master(1,1), 'in processor = ',i
		    !WRITE(*,*)'numraycoords_master is ', numraycoords_master, 'in processor = ',i
		    !WRITE(*,*)'raycoords_master is ', raycoords_master(1,1), 'in processor = ',i
		    !WRITE(*,*)'crazyray_master is ', crazyray_master, 'in processor = ',i
   		    !WRITE(*,*)'raycoords is ', crazyray_master, 'in processor = ',i
		    !WRITE(*,*)'crazyray_master is ', crazyray_master, 'in processor = ',i
		    !WRITE(*,*)'Partial numraycoords at master is', numraycoords, 'processor', taskid
		    !WRITE(*,*)'Partial size of npoints at master is ', SHAPE(npoints)
		    !WRITE(*,*)'Partial size of raycoords', SHAPE(raycoords)
	        END DO

                !WRITE(*,*) 'Finished parallelization'

		CALL cpu_time(t2_rays_test) ! TOC - Stop counting time 
		! Calculate cumulative time spent in modrays routine when parallelized
  		t_rays_test =  t_rays_test  +  (t2_rays_test-t1_rays_test)
	           
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Inner parallelization ends here.  
		! We obtain as a result:
		! npoints, raycoords, crazyray.     	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!Here finally we compare the results of both algorithms 

		IF(TEST_MODRAYS.EQ.1) THEN 

		    WRITE(*,*)'crazyray_master is ', crazyray_master
  		    WRITE(*,*) 'crazyray is', crazyray
		    WRITE(*,*)'nru_master is ', offsetpts_end
		    !WRITE(*,*)'Size of (parallel) npoints at master is ', SHAPE(npoints)
		    !WRITE(*,*)'Size of (parallel) raycoords at master is', SHAPE(raycoords)
		    WRITE(*,*)'original crazyray is ', crazyray_base
		    WRITE(*,*)'original nru is ', nru_base
		    !WRITE(*,*)'Size of (original) npoints is ', SHAPE(npoints_base)
		    !WRITE(*,*)'Size of (original) raycoords is', SHAPE(raycoords_base)
		    !WRITE(*,*)

		    npoints_base = npoints_base - npoints
		    raycoords_base = raycoords_base - raycoords 
		    npoints_diff1 = abs( sum( npoints_base(:, 1) ) ) 
                   npoints_diff2 = abs( sum( npoints_base(:, 2) ) )
		    raycoords_diffx = abs( sum( raycoords_base(:, 1) ) )
		    raycoords_diffy = abs( sum( raycoords_base(:, 2) ) )

		    WRITE(*,*) 'Results parallelization sample', 0
		    WRITE(*,*) 'The sum of the differences is (', npoints_diff1,',', npoints_diff2, ')' 
		    WRITE(*,*) 'The sum of the differences is (', raycoords_diffx,',', raycoords_diffy, ')' 

		END IF 

	        !WRITE(*,*) 'Finished comparison'
  	        !WRITE(*,*) 'crazyray is', crazyray

		! This should not happen, but just in case terminate program if there are problems with the last accepted model of the previous run (which is the first model of this run).
		IF (runnum.GT.1.AND.crazyray.GT.0) THEN
			WRITE(*,*)'Problem tracing rays for the first velocity model (i.e. last model of previous iteration) on rank',rank+1,'!! Terminating program!!'
			CALL MPI_FINALIZE(ierror) !! uncomment for former MPI parallelization
			STOP
		END IF
	END DO
END IF


!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)'------------------------------------------------------------------------'

!WRITE(*,*) 'We now have an INITIAL velocity model made of randomly-located Voronoi cells with random velocity values'


!*************************************************************************************
! If "update raypaths at each proposed model" is selected, then calculate new raypaths
!-------------------------------------------------------------------------------------
IF(urps.EQ.2) THEN
	
	IF(wurf.EQ.-1.OR.wurf.EQ.rank+1) THEN
		rayupdout = 'rayupdate_p' // trim(rankc) // '.dat'
		OPEN(UNIT=40,FILE=rayupdout,FORM='unformatted',STATUS='replace')
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
		 
		IF(wurf.EQ.-1.OR.wurf.EQ.rank+1) THEN
			WRITE(40)raylength
			WRITE(40)raypoints(1,1),raypoints(1,2)
		END IF
			
		DO i=2,raylength
		
			IF(wurf.EQ.-1.OR.wurf.EQ.rank+1) THEN
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
			
		END DO
			
		dis(rayid)=angle*er

		IF (sigdep.EQ.2) THEN
			srdist(rayid)=angle*180/pii !*er
		END IF

		DEALLOCATE(Mi)
				
		DEALLOCATE(raypoints)
				
	END DO
	
	IF(wurf.EQ.-1.OR.wurf.EQ.rank+1) THEN
		CLOSE(40)
	END IF
	
	np=sum(npoints(:,1),1)-nru
		
	DEALLOCATE(raycoords)
	DEALLOCATE(npoints)
	DEALLOCATE(VELI)
	

!-------------------------------------------------------------------------------------
! ...otherwise read the ray geometries from the input file rayini
!-------------------------------------------------------------------------------------
ELSE IF(urps.EQ.0.OR.urps.EQ.1) THEN

	IF(urps.EQ.1) THEN
		IF(ALLOCATED(VELU)) DEALLOCATE(VELU)
		ALLOCATE(VELU(nvpr+2,nvtr+2))
		VELU=0	! This is the matrix containing the model velocities that will be used to update the raypaths
	END IF
	

	IF(ALLOCATED(M)) DEALLOCATE(M)
	IF(ALLOCATED(nnM)) DEALLOCATE(nnM)
	ALLOCATE(M(npoint_max,3))
	ALLOCATE(nnM(nray_max))
	
	! The input file rayini contains the initial geometry of rays between all station pairs.
	! If the raypaths have been updated, then file rayupdout is used instead.
	! The final goal is to store the length of ray i in cell j. 
	! This is given by matrix RayCel(i,j). 
	rayupdout = 'rayupdate_p' // trim(rankc) // '.dat'
	
	IF (urps.EQ.0.OR.(urps.EQ.1.AND.runnum.EQ.1)) OPEN(UNIT=4,FILE=rayini,FORM='unformatted',STATUS='old')
	IF (urps.EQ.1.AND.runnum.GT.1) OPEN(UNIT=4,FILE=rayupdout,FORM='unformatted',STATUS='old')
	
	!IF(rank.EQ.0) WRITE(*,*)'Reading number of source-station pairs'
	READ(4)nraysr ! number of source-station pairs
	
	! Check if nraysr is consistent with the number of source-station pairs.
	IF (nraysr.NE.nsou*nrec) THEN
		IF(rank.EQ.0) WRITE(*,*)'HO HO! Problem 1 here: the number of source-station pairs from the raypath file and from the source and receiver files are not consistent!'
		STOP
	END IF
	
	! Store the ray geometries in matrix M
	
	! column 1: phi (longitude, z)
	! column 2: theta (latitude, x)
	! column 3: distance between the point and the next point. This distance can vary, and is needed for the computation of matrix RayCel
	
	! We start from ray 1 and trace each ray from the receiver to the source
	
	! lengthray(nr) = number of points in ray nr
	
	p=0
	nr=0
	node=1
	RayCel=0
	M=0
	DO nrr=1,nraysr
			
		!IF(rank.EQ.0) WRITE(*,*)'Reading length of ray'
		READ(4)lengthray(nrr) 				! number of points in the ray
		angle=0
		!IF(rank.EQ.0) WRITE(*,*)'Reading ray points'
		READ(4)ptlat,ptlong 				! read coordinates of the first point of the ray
		point(2) = ptlat
		point(1) = ptlong
		!IF(rank.EQ.0) WRITE(*,*)'Done reading ray points'
		point(1) = point(1)*pii/180
		point(2) = point(2)*pii/180

 		IF (RV(nrr).eq..true.) THEN 			! If the ray is valid 

			nr = nr+1 				! nr is the number of valid rays
			lengthrayv(nr)=lengthray(nrr)
			
			! nnM stores the number of points in each valid ray
			nnM(nr) = p+1 
			
			DO i=1,lengthray(nrr)-1
				!IF(rank.EQ.0) WRITE(*,*)'Looping over ray points'
				p=p+1
				!IF(rank.EQ.0) WRITE(*,*)'Reading ray points'
				READ(4)ptlat,ptlong
				
				M(p,2) = ptlat
				M(p,1) = ptlong
				M(p,1) = M(p,1)*pii/180
				M(p,2) = M(p,2)*pii/180
				tetaS = M(p,2)				! latitude
				phiS = M(p,1)				! longitude
				IF (i.ne.1) THEN
					tetaR = M(p-1,2)		! latitude
					phiR = M(p-1,1)			! longitude
				ELSE
					tetaR = point(2)		! latitude
					phiR = point(1)			! longitude
				END IF 
				
				! Do some trigonometry on the sphere to get angle alpha between two consecutive points.
				! Use the Haversine formula (better conditioned for small distances than the spherical law of cosines)
                        	alpha=2*asin(sqrt((sin((tetaS-tetaR)/2))**2+(cos(tetaR))*(cos(tetaS))*(sin((phiS-phiR)/2))**2))

				! Store the distance between 2 points in M(:,3)
				M(p,3)=alpha*er
				angle = angle + alpha
				p1 = (180/pii)*M(p,1)
				p2 = (180/pii)*M(p,2)
	
				! Find in which Voronoi cell the current point is (i.e. to which node the current point is closest)
				CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro(1:2,:)),nnn,nnlist,walk)

				! Compute the distance of ray 'nr' in cell 'node' and store it in Raycel(nr,node) 
				RayCel(nr,node)=RayCel(nr,node)+M(p,3)
			END DO

			dis(nr)=angle*er

			IF (sigdep.EQ.2) THEN
				srdist(nr)=angle*180/pii !*er
			END IF

		ELSE !If the ray is not valid, read points but forget them
			DO i=1,lengthray(nrr)-1
				READ(4)ptlat,ptlong
			END DO
		END IF
 	END DO

	CLOSE(4)

	M(:,1)=(180/pii)*M(:,1)
	M(:,2)=(180/pii)*M(:,2)
	
	np = p !total number of points
	
	! Check that the number of valid rays nr corresponds to the one given by 'otimes.dat'.
	IF (nr.NE.nrays) THEN
		IF(rank.EQ.0) WRITE(*,*)'HO HO! Problem 2 here: the number of valid rays from the raypath file and from the otimes file are not consistent!'
		CALL MPI_FINALIZE(ierror) !! uncomment for former MPI parallelization
		STOP
	END IF
	
	! Check that the total number of points nb does not exceed the limit given in the input file as npoint_max
	IF(np.GT.npoint_max) THEN
  		IF(rank.EQ.0) WRITE(*,*)'The number of ray points exceeds the maximum number of points npoint_max set in the input file!!'
		CALL MPI_FINALIZE(ierror) !! uncomment for former MPI parallelization
		STOP
	END IF
	
	DEALLOCATE(lengthray)
	DEALLOCATE(RV)

!-------------------------------------------------------------------------------------
ELSE

	IF(rank.EQ.0) WRITE(*,*)'Invalid entry at line 69 of the input file!'
	CALL MPI_FINALIZE(ierror) !! un comment for former MPI parallelization
	STOP

END IF
!*************************************************************************************

!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
WRITE(*,*)'Max ray length',MAXVAL(SUM(RayCel,2)),'on rank',rank+1
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
WRITE(*,*)'Min ray length',MINVAL(SUM(RayCel,2)),'on rank',rank+1

!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
WRITE(*,*)'Total number of ray points:',np,'on rank',rank+1
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)'------------------------------------------------------------------------'


!****************************************************************************
!                             Initialize sigma
!****************************************************************************
! If this is not the first run of the chain...
IF (runnum.NE.1) THEN
	! ...read sigma parameters of last accepted model of previous run from file...
	OPEN(UNIT=84,FILE=sigma_tmp,FORM='unformatted',STATUS='old')
	IF (sigdep.EQ.0) THEN
		DO nss=1,nds
			READ(84)sigma(1,nss)
		END DO
	ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
		DO nss=1,nds
			READ(84)aa(nss),bb(nss)
		END DO
		DO nr=1,nrays
			sigma(nr,1) = aa(dsstat(nr))*srdist(nr)+bb(dsstat(nr))
		END DO
	ELSE IF (sigdep.EQ.3) THEN
		DO nss=1,nds
			READ(84)lambda(nss)
		END DO
		DO nr=1,nrays
			sigma(nr,1) = lambda(dsstat(nr))*dat(nr,6)
		END DO
	ELSE IF(sigdep.EQ.4) THEN
		sigma(:,1) = dat(:,6)
		sigma_min = MINVAL(sigma(:,1))
		sigma_max = MAXVAL(sigma(:,1))
	END IF
	CLOSE(84)
!...otherwise, if this is the first run of the chain...
ELSE IF (runnum.EQ.1) THEN
	! ...initialise a random sigma
	IF(sigdep.EQ.0) THEN
		DO nss=1,nds
			sigma(1,nss) = sigma_min+grnd()*(sigma_max-sigma_min)
			randcount=randcount+1
		END DO
	ELSE IF(sigdep.EQ.1.OR.sigdep.EQ.2) THEN
		DO nss=1,nds
			aa(nss) = aa_min+grnd()*(aa_max-aa_min)
			bb(nss) = bb_min+grnd()*(bb_max-bb_min)
			randcount=randcount+2
		END DO
		DO nr=1,nrays
			sigma(nr,1) = aa(dsstat(nr))*srdist(nr)+bb(dsstat(nr))
		END DO
	ELSE IF(sigdep.EQ.3) THEN
		DO nss=1,nds
			lambda(nss) = lambda_min+grnd()*(lambda_max-lambda_min)
			randcount=randcount+1
		END DO
		DO nr=1,nrays
			sigma(nr,1) = lambda(dsstat(nr))*dat(nr,6)
		END DO
	ELSE IF(sigdep.EQ.4) THEN
		sigma(:,1) = dat(:,6)
		sigma_min = MINVAL(sigma(:,1))
		sigma_max = MAXVAL(sigma(:,1))
	END IF
END IF


!****************************************************************************
!                  Get likelihood for the initial model
!****************************************************************************

! ttime(nr) is the estimated travel time for ray nr
! RayCel(nr,i) is the length of ray nr in Voronoi Cell i
! Voro(3,i) is the velocity in Voronoi cell i
! like is the Gaussian or Laplacian likelihood
! The misfit is normalised by the data noise sigma
 
like=0
misfit=0
ttime=0
DO nr=1,nrays
	DO i=1,ncell
		ttime(nr) = ttime(nr) + RayCel(nr,i)/Voro(3,i)
	END DO
	IF(sigdep.EQ.0) THEN
		IF (uglpd.EQ.0) THEN
			like = like + (ttime(nr)-dat(nr,5))**2/(2*(sigma(1,(dsstat(nr)))**2))
		ELSE IF (uglpd.EQ.1) THEN
			like = like + ABS(dat(nr,5)-ttime(nr))/sigma(1,(dsstat(nr)))
		END IF
		misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(1,(dsstat(nr))))**2
	ELSE IF(sigdep.NE.0) THEN
		IF (uglpd.EQ.0) THEN
			like = like + (ttime(nr)-dat(nr,5))**2/(2*(sigma(nr,1)**2))
		ELSE IF (uglpd.EQ.1) THEN
			like = like + ABS(dat(nr,5)-ttime(nr))/sigma(nr,1)
		END IF
		misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(nr,1))**2
	END IF
END DO
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
WRITE(*,*)'Initial average misfit for processor',rank+1,':',misfit/nrays
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)'------------------------------------------------------------------------'
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
IF(rank.EQ.0) WRITE(*,*)
!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization

! If this is the first run of the chain, save initial model
IF (runnum.EQ.1) THEN
	firstVoro=Voro
	samplesNcel(1)=ncell
	IF (sigdep.EQ.0) THEN
		DO nss=1,nds
			samplesSigmas(1,nss)=sigma(1,nss)
		END DO
	ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
		DO nss=1,nds
			samplesSigmas(1,nss)=aa(nss)
			samplesSigmas(1,nss+nds)=bb(nss)
		END DO
	ELSE IF (sigdep.EQ.3) THEN
		DO nss=1,nds
			samplesSigmas(1,nss)=lambda(nss)
		END DO
	END IF
	samplesMisfit(1) = misfit
END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***************************************************************************
!                                                                          !
!    Start the sampling of the posterior distribution with the rj_MCMC     !
!                                                                          !
!***************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ount=0 
sample=0

IF(test.EQ.1) THEN
	CALL SYSTEM('mv First.vor initial.vor')
	CALL SYSTEM('mv First.vtx initial.vtx')
	IF(urps.EQ.2) THEN
		CALL SYSTEM('mv rayupdate.dat raysubi.out')
	END IF
END IF


DO WHILE (ount.LT.nsample)	! <<<----------------------------------- LOOP OVER SAMPLES STARTS HERE !!! -----------------------------------
!write(*,*)'Sample',sample,'on processor',ra
ount=ount+1
Voro_prop = Voro
Voro_pprop = Voro
ncell_prop = ncell
sigma_prop = sigma
RayCel_prop = RayCel
RayCel_pprop = RayCel
ttime_prop = ttime
ttime_pprop = ttime
like_prop = 0
like_pprop = 0
IF (sigdep.EQ.2) srdist_prop = srdist
IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
	aa_prop = aa
	bb_prop = bb
ELSE IF (sigdep.EQ.3) THEN
	lambda_prop = lambda
END IF
a=0
birth=0
death=0
move=0
velch=0
sigmach=0

IF(crazyray.GT.0) THEN
	WRITE(*,*)'Loop cycled! Proposing new sample',ount
	WRITE(*,*)
	!STOP
END IF

!*************************************************************
!                   Propose a new model
!*************************************************************

 cellogic=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 			Choose randomly between 4 or 5 types of moves
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
u=grnd()
randcount=randcount+1

IF (sigdep.NE.4) THEN
	IF (mod(ount,2).EQ.0) THEN
		IF (u<0.333) THEN
			! BIRTH - add a cell
			birth=1
			! check for boundary (do not allow more than ncell_max cells)
			IF (ncell.NE.ncell_max) THEN
				a=1
			END IF
		ELSE IF (u<0.666) THEN
			! DEATH - remove a cell
			death=1
			! check for boundary (do not allow less than ncell_min cells)
			IF (ncell.NE.ncell_min) THEN
				a=1
			END IF
		ELSE
			! MOVE - move a cell
			move=1 
		END IF
	ELSE
		IF (u<0.500) THEN
			! VELOCITY CHANGE - change the velocity of a cell
			velch=1
		ELSE
			! SIGMA - change sigma by either changing its gradient (aa) or its y-intercept (bb)
			sigmach=1
		END IF  
	END IF
ELSE IF (sigdep.EQ.4) THEN
	IF (mod(ount,2).EQ.0) THEN
		IF (u<0.333) THEN
			! BIRTH - add a cell
			birth=1
			! check for boundary (do not allow more than ncell_max cells)
			IF (ncell.NE.ncell_max) THEN
				a=1
			END IF
		ELSE IF (u<0.666) THEN
			! DEATH - remove a cell
			death=1
			! check for boundary (do not allow less than ncell_min cells)
			IF (ncell.NE.ncell_min) THEN
				a=1
			END IF
		ELSE
			! MOVE - move a cell
			move=1 
		END IF
	ELSE
		! VELOCITY CHANGE - change the velocity of a cell
		velch=1
	END IF  
END IF


! Now, depending on move type, update the model

!***************************************************************************************	 
!                  Change the cell geometry (cell birth, death, move)	
!***************************************************************************************	
IF (birth.EQ.1.OR.death.EQ.1.OR.move.EQ.1) THEN

	IF (birth.EQ.1) THEN ! Add one cell
		
		IF (a.EQ.1) THEN ! if inside the prior...
		
			! We're adding a cell
			ncell_prop = ncell+1;
			
			! Generate a new random location
			Voro_prop(1:2,ncell_prop) = [longmin+grnd()*(longmax-longmin),latmin+grnd()*(latmax-latmin)]
			randcount=randcount+2
			
			! Calculate triangulation of the previous (accepted) model
			CALL delaun (DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			
			CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
				
			! Find what node (in the previous model) the new cell is closest to, so that we can get a starting velocity value for the new cell
			node=1
			CALL find_node2D(DBLE([Voro_prop(1,ncell_prop),Voro_prop(2,ncell_prop)]),node,DBLE(Voro(1:2,:)),nnn,nnlist,walk)
		
			! Propose a velocity for the new cell
			Voro_prop(3,ncell_prop) = Voro(3,node) + GASDEV(randcount)*sigmav
			
			! If the proposed velocity is outside the bounds given by theta, then the model is automatically rejected (a=0)
			IF (abs(mean-Voro_prop(3,ncell_prop)).GT.theta) a=0
			
		END IF

		IF (a.EQ.1) THEN ! if inside the prior...
			
			! Set values of Voro_prop beyond ncell_prop to zero
			Voro_prop(:,ncell_prop+1:ncell_max)=0
			
			! Calculate "prob" (one of the terms needed to evaluate alpha)
			prob=log(sigmav*sqrt(2*pii))+((Voro_prop(3,ncell_prop)-Voro(3,node))**2/(2*sigmav**2))
			
			! Compute the proposed Delaunay triangulation: nnn_prop, nnlist_prop  
			CALL delaun (DBLE(Voro_prop(1:2,:)),ncell_prop,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			IF (ptsok.EQ.1) THEN
				WRITE(*,*)'All',ncell_prop,'points in a line for sample',ount,'on processor',rank+1
				WRITE(*,*)'Sample number',ount,'discarded...',birth,death,move,velch,sigmach,ptsok
				!STOP
				ount=ount-1
				CYCLE
			END IF
			
			CALL build_nv(ncell_prop,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
			
			! Find the neighbours of the added cell - these cells are the ones which are affected by the geometry change
			istart = nnn(ncell_prop)
			iend = nnn(ncell_prop+1)-1
			
			! Update cellogic 
			DO k=istart,iend
				IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
			END DO
			
			! Set the length of the rays in the added cell to zero (temporary measure!)
			RayCel_prop(:,ncell_prop) = 0
			RayCel_prop(:,ncell_prop+1:ncell_max) = 0
			
			!!! ******************************************************************************************************
			!!! Original code modified by EG - option to update raypaths added
			
			IF(urps.EQ.2) THEN
			
				! Update raypaths by calling subroutine updaterays
				IF (uar.EQ.1) cellogic=.true.
				ind=1
				!write(*,*)'BIRTH'


		                !!! Parallelization starts here 

				!t1_rays = 0
		  	        !t2_rays = 0
				!CALL cpu_time(t1_rays) ! TOC - Start counting time

				CALL updaterays(RayCel_prop,RayCel,Voro_prop,ncell,ncell_prop,cellogic,ind,rank+1,ount, t_uprays, t_uprays_test)
		
				!CALL cpu_time(t2_rays) ! TOC - Stop counting time 
				! Calculate cumulative time spent in modrays routine 

		  		t_rays_test_base =  t_rays_test_base  +  t_uprays
				t_rays_test =  t_rays_test  +  t_uprays_test


				
				! If there are crazy raypaths for this sample, cycle the loop and propose a new sample
				IF(crazyray.GT.0) THEN
					WRITE(*,*)'Sample number',ount,'discarded...'
					!STOP
					ount=ount-1
					CYCLE
				END IF
		
			ELSE 
				
				! Do not update the raypaths
				
				node=1
				DO nr=1,nrays
					DO i=1,ncell
						IF ((cellogic(i)).and.(RayCel(nr,i).ne.0)) THEN
							
							RayCel_prop(nr,:) = 0
							
							jstart = nnM(nr)
							jend =  nnM(nr)+lengthrayv(nr)-2
							
							DO j=jstart,jend
								CALL find_node2D(DBLE(M(j,1:2)),node,DBLE(Voro_prop(1:2,:)),nnn,nnlist,walk)
								RayCel_prop(nr,node)=RayCel_prop(nr,node)+M(j,3)
							END DO
							
							EXIT
							
						END IF
					END DO
				END DO
			END IF
			!!! ******************************************************************************************************
		
		ELSE	! if outside the prior... (a=0)
		
			prob=0.0
		
		END IF ! a=1 (if values we selected are inside the prior)
	
	ELSE IF (death.EQ.1) THEN ! Remove one cell
	
		IF (a.EQ.1) THEN	! if inside the prior...
			
			! We've lost a cell
			ncell_prop = ncell-1
			
			! Choose a cell to delete
			ind = ceiling(grnd()*ncell)
			randcount=randcount+1
			IF (ind.EQ.0) ind=1
			!WRITE(*,*)ind
			
			! Replace the deleted cell by the last one
			Voro_prop(1:3,ind) = Voro_prop(1:3,ncell)
			Voro_prop(1:3,ncell_prop+1:ncell_max) = 0
			
			! Update cellogic for the deleted cell
			cellogic(ind)=.true. 
		
			! Calculate triangulation of the previous (accepted) model
			CALL delaun (DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			
			CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
			
			! Find the neighbours of the deleted cell - these cells are the ones which are affected by the geometry change
			istart = nnn(ind)
			iend = nnn(ind+1)-1
			
			! Update cellogic for the neighbouring cells
			DO k=istart,iend
				IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
			END DO
			
			! Compute the proposed Delaunay triangulation: nnn_prop, nnlist_prop
			CALL delaun (DBLE(Voro_prop(1:2,:)),ncell_prop,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			IF (ptsok.EQ.1) THEN
				WRITE(*,*)'All',ncell_prop,'points in a line for sample',ount,'on processor',rank+1
				WRITE(*,*)'Sample number',ount,'discarded...',birth,death,move,velch,sigmach,ptsok
				!STOP
				ount=ount-1
				CYCLE
			END IF
			
			CALL build_nv(ncell_prop,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
			
			! Find what node (in the previous model) the deleted cell is closest to, so that we can get a new velocity value
			node=1
			CALL find_node2D(DBLE([Voro(1,ind),Voro(2,ind)]),node,DBLE(Voro_prop(1:2,:)),nnn,nnlist,walk)
			
			! Calculate "prob" (one of the terms needed to evaluate alpha)
			prob=log(1/(sigmav*sqrt(2*pii)))+(-(Voro_prop(3,node)-Voro(3,ind))**2/(2*sigmav**2))
			
			! Replace the length of the rays in the deleted cell by the length of the rays in the last one
			RayCel_prop(:,ind) = RayCel_prop(:,ncell)
			RayCel_prop(:,ncell_prop+1:ncell_max) = 0
			
			!!! ******************************************************************************************************
			!!! Original code modified by EG - option to update raypaths added
			
			IF(urps.EQ.2) THEN
			
				! Update raypaths by calling subroutine updaterays
				IF (uar.EQ.1) cellogic=.true.
				!write(*,*)'DEATH'

		                !!! Parallelization starts here 

				!t1_rays = 0
		  	        !t2_rays = 0
				!CALL cpu_time(t1_rays) ! TOC - Start counting time

				CALL updaterays(RayCel_prop,RayCel,Voro_prop,ncell,ncell_prop,cellogic,ind,rank+1,ount, t_uprays, t_uprays_test)		

				!CALL cpu_time(t2_rays) ! TOC - Stop counting time 
				! Calculate cumulative time spent in modrays routine 
		  				  		
		  		t_rays_test_base =  t_rays_test_base  +  t_uprays
				t_rays_test =  t_rays_test  +  t_uprays_test




			
				! If there are crazy raypaths for this sample, cycle the loop and propose a new sample
				IF(crazyray.GT.0) THEN
					WRITE(*,*)'Sample number',ount,'discarded...'
					!STOP
					ount=ount-1
					CYCLE
				END IF
				
			ELSE 
				
				! Do not update the raypaths
				
				node=1			
				DO nr=1,nrays
					DO i=1,ncell
						IF ((cellogic(i)).and.(RayCel(nr,i).ne.0)) THEN
						
							RayCel_prop(nr,:) = 0
							
							jstart = nnM(nr)
							jend = nnM(nr)+lengthrayv(nr)-2
							
							DO j=jstart,jend
								CALL find_node2D(DBLE(M(j,1:2)),node,DBLE(Voro_prop(1:2,:)),nnn,nnlist,walk)
								RayCel_prop(nr,node)=RayCel_prop(nr,node)+ M(j,3)
							END DO
							
							EXIT
							
						END IF
					END DO
				END DO
			END IF
			!!! ******************************************************************************************************
			
		ELSE	! if outside the prior... (a=0)
		
			prob=0.0
		
		END IF ! a=1 (if values we selected are inside the prior)
	
	ELSE IF (move.EQ.1) THEN ! Move the center of one cell
		
		! The number of cells has not changed
		ncell_prop=ncell
		
		! Choose a cell to move
		ind = ceiling(grnd()*ncell)
		randcount=randcount+1
		IF (ind.EQ.0) ind=1
		
		! Choose new coordinates for the center of the cell
		p1 = Voro(1,ind)+(GASDEV(randcount)*pd*(longmax-longmin)/100)
		p2 = Voro(2,ind)+(GASDEV(randcount)*pd*(latmax-latmin)/100)
		
		! If the values fall within the accepted range...
		IF ((p1.GE.longmin).AND.(p2.GE.latmin).AND.(p1.LE.longmax).AND.(p2.LE.latmax)) THEN
			
			a=1
			Voro_prop(1,ind) = p1
			Voro_prop(2,ind) = p2
			
		END IF
		
		IF(a.EQ.1) THEN
		
			! Calculate triangulation of the previous (accepted) model (next 6 lines added by EG)
			CALL delaun (DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			
			CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
		
			! Update cellogic for the selected cell
			cellogic(ind)=.true.
			
			! Find the neighbours of the selected cell before moving - these cells are affected by the geometry change 
			istart = nnn(ind)
			iend = nnn(ind+1)-1
			
			! Update cellogic for the neighbouring cells
			DO k=istart,iend
				IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
			END DO
			
			! Compute the proposed Delaunay triangulation: nnn_prop, nnlist_prop
			CALL delaun (DBLE(Voro_prop(1:2,:)),ncell_prop,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			IF (ptsok.EQ.1) THEN
				WRITE(*,*)'All',ncell_prop,'points in a line for sample',ount,'on processor',rank+1
				WRITE(*,*)'Sample number',ount,'discarded...',birth,death,move,velch,sigmach,ptsok
				!STOP
				ount=ount-1
				CYCLE
			END IF
			
			CALL build_nv(ncell_prop,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
			
			! Find the neighbours of the selected cell after moving - these cells are affected by the geometry change 
			istart = nnn(ind)
			iend = nnn(ind+1)-1
			
			! Update cellogic for the neighbouring cells
			DO k=istart,iend
				IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
			END DO
			
			!!! ******************************************************************************************************
			!!! Original code modified by EG - option to update raypaths added
			
			IF(urps.EQ.2) THEN
			
				! Update raypaths by calling subroutine updaterays
				IF (uar.EQ.1) cellogic=.true.
				!write(*,*)'MOVE'

				!t1_rays = 0
		  	        !t2_rays = 0
				!CALL cpu_time(t1_rays) ! TOC - Start counting time

				CALL updaterays(RayCel_prop,RayCel,Voro_prop,ncell,ncell_prop,cellogic,ind,rank+1,ount, t_uprays, t_uprays_test)

				!CALL cpu_time(t2_rays) ! TOC - Stop counting time 
				! Calculate cumulative time spent in modrays routine 

		  		t_rays_test_base =  t_rays_test_base  +  t_uprays
				t_rays_test =  t_rays_test  +  t_uprays_test


		
				! If there are crazy raypaths for this sample, cycle the loop and propose a new sample
				IF(crazyray.GT.0) THEN
					WRITE(*,*)'Sample number',ount,'discarded...'
					!STOP
					ount=ount-1
					CYCLE
				END IF
				
			ELSE 
				
				! Do not update the raypaths
				
				node=1
				pp=0
				DO nr=1,nrays
					DO i=1,ncell
						IF ((cellogic(i)).and.(RayCel(nr,i).ne.0)) THEN
							
							pp=pp+1
							
							RayCel_prop(nr,:) = 0
							
							jstart = nnM(nr)
							jend =  nnM(nr)+lengthrayv(nr)-2
							
							DO j=jstart,jend
								CALL find_node2D(DBLE(M(j,1:2)),node,DBLE(Voro_prop(1:2,:)),nnn,nnlist,walk)
								RayCel_prop(nr,node)=RayCel_prop(nr,node)+ M(j,3)
							END DO
							
							EXIT
							
						END IF
					END DO
				END DO
			END IF
			!!! ******************************************************************************************************
			
		END IF ! a=1 (if values we selected are inside the prior)
	END IF
	
	!***************************************************************************************
	!         After changing the cell geometry, get the misfit of the proposed model
	!***************************************************************************************

	! ttime_prop(nr) is the estimated travel time for ray nr
	! RayCel_prop(nr,i) is the length of ray nr in Voronoi Cell i
	! Voro_prop(3,i) is the velocity in Voronoi cell i
	! like_prop is the Gaussian or Laplacian likelihood

	! The parameter a indicates whether we are inside or outside the prior.
	! If we are outside, no need the compute the new misfit since the model is
	! going to be rejected anyway.
	
		
	IF (a.EQ.1) THEN ! compute the likelihood
		
		like_prop=0

		DO nr=1,nrays
			
			ttime_prop(nr)=0
			
			DO i=1,ncell_prop
				ttime_prop(nr) = ttime_prop(nr) + (RayCel_prop(nr,i))/(Voro_prop(3,i))
			END DO
			
			IF(sigdep.EQ.0) THEN
				IF (uglpd.EQ.0) THEN
					like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(1,dsstat(nr))**2))
				ELSE IF (uglpd.EQ.1) THEN
					like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(1,dsstat(nr))
				END IF
			ELSE IF(sigdep.NE.0) THEN
				IF (uglpd.EQ.0) THEN
					like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(nr,1)**2))
				ELSE IF (uglpd.EQ.1) THEN
					like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(nr,1)
				END IF
			END IF
		
		END DO
		
	ELSE ! (a=0) if we are outside the prior
		
		like_prop = like
		 
	END IF
	

!***************************************************************************************	 
!              Change a velocity value (the cell geometry does not change)
!***************************************************************************************
ELSE IF (velch.EQ.1) THEN

	! The number of cells has not changed
	ncell_prop=ncell 
	
	! Pick a cell
	ind = ceiling(grnd()*ncell)
	randcount=randcount+1
	IF (ind.EQ.0) ind=1
	
	! Propose a new velocity for cell ind
	Voro_prop(3,ind) = Voro(3,ind)+GASDEV(randcount)*pv
	
	a=1
	IF (abs(mean-Voro_prop(3,ind)).GT.theta) a=0	! if we are outside the acceptable velocity range
	
	!*************************************************************************************
	! After changing the velocity in one cell, get the misfit of the proposed model
	!*************************************************************************************
	IF (a.EQ.1) THEN	! if we are inside the prior...
		
		!!! ******************************************************************************************************
		!!! Original code modified by EG - option to update raypaths added
			
		IF(urps.EQ.2) THEN
			
			! Update raypaths by calling subroutine updaterays
			!write(*,*)'CHANGE VELOCITY VALUE'
			
			! Update cellogic for the selected cell
			cellogic(ind)=.true.
			IF (uar.EQ.1) cellogic=.true.
			
			! Compute the proposed Delaunay triangulation : nnn_prop,nnlist_prop
			CALL delaun (DBLE(Voro_prop(1:2,:)),ncell_prop,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			IF (ptsok.EQ.1) THEN
				WRITE(*,*)'All',ncell_prop,'points in a line for sample',ount,'on processor',rank+1
				WRITE(*,*)'Sample number',ount,'discarded...',birth,death,move,velch,sigmach,ptsok
				!STOP
				ount=ount-1
				CYCLE
			END IF
			
			CALL build_nv(ncell_prop,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
			
			! Find the neighbours of the selected cell after moving - these cells are affected by the geometry change 
			istart = nnn(ind)
			iend = nnn(ind+1)-1
			
			! Update cellogic for the neighbouring cells
			DO k=istart,iend
				IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
			END DO


			!t1_rays = 0
	  	        !t2_rays = 0
			!CALL cpu_time(t1_rays) ! TOC - Start counting time
			
			CALL updaterays(RayCel_prop,RayCel,Voro_prop,ncell,ncell_prop,cellogic,ind,rank+1,ount, t_uprays, t_uprays_test)

			!CALL cpu_time(t2_rays) ! TOC - Stop counting time 
			! Calculate cumulative time spent in modrays routine 

	  		t_rays_test_base =  t_rays_test_base  +  t_uprays
			t_rays_test =  t_rays_test  +  t_uprays_test

			
			! If there are crazy raypaths for this sample, cycle the loop and propose a new sample
			IF(crazyray.GT.0) THEN
					WRITE(*,*)'Sample number',ount,'discarded...'
					!STOP
				ount=ount-1
				CYCLE
			END IF
				
			like_prop=0
		
			DO nr=1,nrays
			
				ttime_prop(nr)=0
				
				DO i=1,ncell_prop
					ttime_prop(nr) = ttime_prop(nr) + (RayCel_prop(nr,i))/(Voro_prop(3,i))
				END DO
				
				IF(sigdep.EQ.0) THEN
					IF (uglpd.EQ.0) THEN
						like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(1,dsstat(nr))**2))
					ELSE IF (uglpd.EQ.1) THEN
						like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(1,dsstat(nr))
					END IF
				ELSE IF(sigdep.NE.0) THEN
					IF (uglpd.EQ.0) THEN
						like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(nr,1)**2))
					ELSE IF (uglpd.EQ.1) THEN
						like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(nr,1)
					END IF
				END IF
			
			END DO
		
		ELSE
		
			like_prop=0
			
			DO nr=1,nrays 	! Just change the rays going through the cell
				
				IF (RayCel(nr,ind).NE.0) THEN
				
					! If the ray goes through the cell changed 'ind' 
					
					ttime_prop(nr)= 0  
			
					DO i=1,ncell 
						ttime_prop(nr) = ttime_prop(nr) + (RayCel(nr,i))/(Voro_prop(3,i))
					END DO
			
				END IF
			
				IF(sigdep.EQ.0) THEN
					IF (uglpd.EQ.0) THEN
						like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(1,dsstat(nr))**2))
					ELSE IF (uglpd.EQ.1) THEN
						like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(1,dsstat(nr))
					END IF
				ELSE IF(sigdep.NE.0) THEN
					IF (uglpd.EQ.0) THEN
						like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(nr,1)**2))
					ELSE IF (uglpd.EQ.1) THEN
						like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(nr,1)
					END IF
				END IF
				
			END DO
		END IF
		!!! ******************************************************************************************************
		
	ELSE ! (a=0) if we are outside the prior
		
		like_prop = like 
	
	END IF

!***************************************************************************************	 
!                  Change the value of sigma (the cell geometry does not change)	
!***************************************************************************************
ELSE IF (sigmach.EQ.1) THEN

	IF (nds.GT.1) THEN
		nss=ceiling(grnd()*nds)
		randcount=randcount+1
		IF (nss.EQ.0) nss=1
		IF (nss.GT.nds) nss=nds
	ELSE
		nss=1
	END IF
	!* Propose a new sigma
	IF (sigdep.EQ.0) THEN
		sigma_prop(1,nss) = sigma(1,nss) + GASDEV(randcount)*ps
	ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
		u=grnd()
		randcount=randcount+1
		IF (u<0.5) THEN
			cha=1
			chb=0
			aa_prop(nss) = aa(nss) + GASDEV(randcount)*pa
			bb_prop = bb
		ELSE
			cha=0
			chb=1
			aa_prop = aa
			bb_prop(nss) = bb(nss) + GASDEV(randcount)*pb
		END IF
		DO nr=1,nrays
			IF (dsstat(nr).EQ.nss) sigma_prop(nr,1) = aa_prop(nss)*srdist(nr) + bb_prop(nss)
		END DO
	ELSE IF (sigdep.EQ.3) THEN
		lambda_prop(nss) = lambda(nss) + GASDEV(randcount)*pl
		DO nr=1,nrays
			IF (dsstat(nr).EQ.nss) sigma_prop(nr,1) = lambda_prop(nss)*dat(nr,6)
		END DO
	END IF

	!* Check its bounds
	a=1
	IF(sigdep.EQ.0) THEN
		IF ((sigma_prop(1,nss).LT.sigma_min).OR.(sigma_prop(1,nss).GT.sigma_max)) a=0
	ELSE IF(sigdep.EQ.1.OR.sigdep.EQ.2) THEN
		IF (aa_prop(nss).LT.aa_min.OR.aa_prop(nss).GT.aa_max.OR.bb_prop(nss).LT.bb_min.OR.bb_prop(nss).GT.bb_max) a=0
	ELSE IF(sigdep.EQ.3) THEN
		IF (lambda_prop(nss).LT.lambda_min.OR.lambda_prop(nss).GT.lambda_max) a=0
	END IF

	IF (a.EQ.1) THEN	! if we are inside the prior

		like_prop=0
			
		IF(sigdep.EQ.0) THEN
			DO nr=1,nrays
				IF (uglpd.EQ.0) THEN
					like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(1,dsstat(nr))**2))
				ELSE IF (uglpd.EQ.1) THEN
					like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(1,dsstat(nr))
				END IF
			END DO
		ELSE IF(sigdep.NE.0) THEN
			DO nr=1,nrays
				IF (uglpd.EQ.0) THEN
					like_prop = like_prop + (ttime_prop(nr)-dat(nr,5))**2/(2*(sigma_prop(nr,1)**2))
				ELSE IF (uglpd.EQ.1) THEN
					like_prop = like_prop + ABS(dat(nr,5)-ttime_prop(nr))/sigma_prop(nr,1)
				END IF	
			END DO
		END IF

	ELSE ! (a=0) if we are outside the prior
		
		like_prop = like 
	
	END IF

END IF ! Geometry or velocity or sigma change



!*****************************************************************************************************
!  Now, depending on the type of move, compute the acceptance term, and see if we accept or reject
!  the model. If we reject the model for a cell move or velocity change, propose a second move 
!  (this is called delayed rejection).
!*****************************************************************************************************

!***************************************************************************************
! The likelihood needs to be normalised by the determinant of the matrix of data errors
!***************************************************************************************
! lgsigma is the log of the ratio of the determinants of the matrix of data errors
IF(sigdep.EQ.0) THEN
	! For a fixed sigma we have
	lgsigma=0
	DO nr=1,nrays
		lgsigma = lgsigma+log(sigma(1,dsstat(nr))/sigma_prop(1,dsstat(nr)))
	END DO
ELSE IF(sigdep.NE.0) THEN
	! For a variable sigma, we have
	lgsigma=0
	DO nr=1,nrays
		lgsigma = lgsigma+log(sigma(nr,1)/sigma_prop(nr,1))
	END DO
END IF

accept=.false.
alpha1=0
rannum=grnd()
randcount=randcount+1
!DRp=.false.
!DRv=.false.

IF(birth.EQ.1) THEN	
	
	IF (a.EQ.1) THEN	! if we are inside the prior evaluate alpha1, otherwise don't even bother
		IF (nodata.EQ.0) THEN
			alpha1=MINVAL([log(one),prob-log(2*theta)-like_prop+like+lgsigma])
		ELSE IF (nodata.EQ.1) THEN
			alpha1=MINVAL([log(one),prob-log(2*theta)])	! here the likelihood is set to unity (we remove the data)
		END IF
	END IF
	! If we are inside the prior and rand<=alpha1 ACCEPT the model
	IF ((a.EQ.1).AND.(log(rannum).LE.alpha1)) accept=.true.

	!write(*,*)'BIRTH ount=',ount+ounti,'a=',a,'prob=',prob,'log(2*theta)=',log(2*theta),'like_prop=',like_prop,'like=',like,'lgsigma=',lgsigma,'sum=',prob-log(2*theta)-like_prop+like+lgsigma,'log(rannum)=',log(rannum),'alpha1=',alpha1,'accept=',accept

ELSE IF(death.EQ.1) THEN

	IF (a.EQ.1) THEN	! if we are inside the prior evaluate alpha1, otherwise don't even bother
		IF (nodata.EQ.0) THEN
			alpha1=MINVAL([log(one),prob+log(2*theta)-like_prop+like+lgsigma])
		ELSE IF (nodata.EQ.1) THEN
			alpha1=MINVAL([log(one),prob+log(2*theta)])	! here the likelihood is set to unity (we remove the data)
		END IF
	END IF
	
	! If we are inside the prior and rand<=alpha1 ACCEPT the model
	IF ((a.EQ.1).AND.(log(rannum).LE.alpha1)) accept=.true.
	
	!write(*,*)'DEATH ount=',ount+ounti,'a=',a,'prob=',prob,'log(2*theta)=',log(2*theta),'like_prop=',like_prop,'like=',like,'lgsigma=',lgsigma,'sum=',prob+log(2*theta)-like_prop+like+lgsigma,'log(rannum)=',log(rannum),'alpha1=',alpha1,'accept=',accept

ELSE IF (move.EQ.1.OR.velch.EQ.1) THEN
	
	IF (a.EQ.1) THEN	! if we are inside the prior evaluate alpha1, otherwise don't even bother
		IF (nodata.EQ.0) THEN
			alpha1=MINVAL([log(one),-like_prop+like+lgsigma])
		ELSE IF (nodata.EQ.1) THEN
			alpha1=log(one)	! here the likelihood is set to unity (we remove the data)
		END IF
	END IF
	
	! If we are inside the prior and rand<alpha1 ACCEPT the model
	IF ((a.EQ.1).AND.(log(rannum).LE.alpha1)) THEN
	
		accept=.true.
		IF (move.EQ.1) AR1p(1)=AR1p(1)+1		! if we are moving a cell
		IF (velch.EQ.1) AR1v(1)=AR1v(1)+1		! if we are changing a cell velocity

		!IF (move.EQ.1) write(*,*)'MOVE ount=',ount+ounti,'a=',a,'like_prop=',like_prop,'like=',like,'lgsigma=',lgsigma,'sum=',-like_prop+like+lgsigma,'log(rannum)=',log(rannum),'alpha1=',alpha1,'accept=',accept

		!IF (velch.EQ.1) write(*,*)'VELCH ount=',ount+ounti,'a=',a,'like_prop=',like_prop,'like=',like,'lgsigma=',lgsigma,'sum=',-like_prop+like+lgsigma,'log(rannum)=',log(rannum),'alpha1=',alpha1,'accept=',accept

	ELSE IF ((log(rannum).GT.alpha1.OR.a.EQ.0).AND.(velch.EQ.1).AND.(DRv.EQ..true.)) THEN ! if we have rejected a velocity value change 

		AR1v(2)=AR1v(2)+1
	
		!write(*,*)'VELCH ount=',ount+ounti,'a=',a,'like_prop=',like_prop,'like=',like,'lgsigma=',lgsigma,'sum=',-like_prop+like+lgsigma,'log(rannum)=',log(rannum),'alpha1=',alpha1,'accept=',accept

!***********************************************************************************************
!         If we reject a proposed velocity move, start delayed rejection on the value
!***********************************************************************************************
	
		a2=0
		
		Voro_pprop(3,ind) = Voro(3,ind)+GASDEV(randcount)*pv2
		IF (sigdep.EQ.2) srdist_prop = srdist
		
		IF (abs(mean-Voro_pprop(3,ind)).LE.theta) a2=1	! if we are inside the acceptable velocity range
			IF (a2.EQ.1) THEN
		
			!!! ******************************************************************************************************
			!!! Original code modified by EG - option to update raypaths added
				
			IF(urps.EQ.2) THEN
				
				! Update raypaths by calling subroutine updaterays
				!write(*,*)'CHANGE VELOCITY VALUE - DR'
				IF (uar.EQ.1) cellogic=.true.

				CALL updaterays(RayCel_pprop,RayCel,Voro_pprop,ncell,ncell_prop,cellogic,ind,rank+1,ount, t_uprays, t_uprays_test)

		  		t_rays_test_base =  t_rays_test_base  +  t_uprays
				t_rays_test =  t_rays_test  +  t_uprays_test


				
				! If there are crazy raypaths for this sample, cycle the loop and propose a new sample
				IF(crazyray.GT.0) THEN
					WRITE(*,*)'Sample number',ount,'discarded...'
					!STOP
					ount=ount-1
					CYCLE
				END IF
				
				like_pprop=0
			
				DO nr=1,nrays
				
					ttime_pprop(nr)=0
					
					DO i=1,ncell
						ttime_pprop(nr) = ttime_pprop(nr) + (RayCel_pprop(nr,i))/(Voro_pprop(3,i))
					END DO
					
					IF(sigdep.EQ.0) THEN
						IF (uglpd.EQ.0) THEN
							like_pprop = like_pprop + (ttime_pprop(nr)-dat(nr,5))**2/(2*(sigma_prop(1,dsstat(nr))**2))
						ELSE IF (uglpd.EQ.1) THEN
							like_pprop = like_pprop + ABS(dat(nr,5)-ttime_pprop(nr))/sigma_prop(1,dsstat(nr))
						END IF
					ELSE IF(sigdep.NE.0) THEN
						IF (uglpd.EQ.0) THEN
							like_pprop = like_pprop + (ttime_pprop(nr)-dat(nr,5))**2/(2*(sigma_prop(nr,1))**2)
						ELSE IF (uglpd.EQ.1) THEN
							like_pprop = like_pprop + ABS(dat(nr,5)-ttime_pprop(nr))/sigma_prop(nr,1)
						END IF
					END IF
				
				END DO
			
			ELSE
				
				! Compute the misfit of the proposed model
	
				like_pprop=0
				
				DO nr=1,nrays ! Just change the rays going through cell 'ind'
					IF (RayCel(nr,ind).NE.0) THEN
						! If the ray goes through the cell changed 'ind' 
						ttime_pprop(nr)= 0  
						DO i=1,ncell 
							ttime_pprop(nr) = ttime_pprop(nr) + (RayCel(nr,i))/(Voro_pprop(3,i))  
						END DO
					END IF
					
					IF(sigdep.EQ.0) THEN
						IF (uglpd.EQ.0) THEN
							like_pprop = like_pprop + (ttime_pprop(nr)-dat(nr,5))**2/(2*(sigma_prop(1,dsstat(nr))**2))
						ELSE IF (uglpd.EQ.1) THEN
							like_pprop = like_pprop + ABS(dat(nr,5)-ttime_pprop(nr))/sigma_prop(1,dsstat(nr))
						END IF
					ELSE IF(sigdep.NE.0) THEN
						IF (uglpd.EQ.0) THEN
							like_pprop = like_pprop + (ttime_pprop(nr)-dat(nr,5))**2/(2*(sigma_prop(nr,1))**2)
						ELSE IF (uglpd.EQ.1) THEN
							like_pprop = like_pprop + ABS(dat(nr,5)-ttime_pprop(nr))/sigma_prop(nr,1)
						END IF
					END IF
					
				END DO
				
			END IF
				
		ELSE ! (a2=0) if we are outside the prior
		
			like_pprop = like_prop 
		
		END IF

		!*********************************!
		! NOW NEED TO COMPILE THE ALPHA 2 !
		!*********************************!
		
		IF (nodata.EQ.0) THEN
			TT1=-like_pprop+like+lgsigma
		ELSE IF (nodata.EQ.1) THEN
			TT1=log(one)	! ignore the data
		END IF
		TT2=(-(Voro_prop(3,ind)-Voro_pprop(3,ind))**2+(Voro_prop(3,ind)-Voro(3,ind))**2)/(2*pv**2)
		
		alpha2=0
		rannum2=grnd()
		randcount=randcount+1
		
		IF (a2.EQ.1.AND.a.EQ.0) THEN		! if the first proposal was rejected because out of bounds, then TT3=0
			TT3t=0
			TT3b=0
			TT3=0
			alpha2=MINVAL([log(one),TT1+TT2])
		ELSE IF (a2.EQ.1.AND.a.EQ.1) THEN	! if both the first and the delayed proposal are within the bounds
			IF (nodata.EQ.0) THEN
				TT3t=log(1-exp(MINVAL([log(one),-like_prop+like_pprop+lgsigma])))
			ELSE IF (nodata.EQ.1) THEN	! shouldn't even get to this point because for no data all m' are accepted
				TT3t=log(1-exp(log(one)))	! ignore the data
			END IF
			TT3b=log(1-exp(alpha1))
			TT3=TT3t-TT3b
			alpha2=MINVAL([log(one),TT1+TT2+TT3])
		END IF

		! Now compare
		IF (a2.EQ.1.AND.log(rannum2).LE.alpha2) THEN
		
			! Accept the second attempt
			accept=.true.
			AR2v(1)=AR2v(1)+1
			
			like_prop = like_pprop
			Voro_prop = Voro_pprop
			ttime_prop = ttime_pprop
			RayCel_prop = RayCel_pprop
			
		ELSE 
			AR2v(2)=AR2v(2)+1
		END IF
			
		!write(*,*)'VELCH2 ount=',ount+ounti,'a=',a,'a2=',a2,'like_pprop=',like_pprop,'TT1=',TT1,'TT2=',TT2,'TT3t=',TT3t,'TT3b=',TT3b,'TT3=',TT3,'log(rannum2)=',log(rannum2),'alpha2=',alpha2,'accept=',accept
	
!*********************END DELAYED REJECTION FOR THE VALUE *************************
!**********************************************************************************
	ELSE IF ((log(rannum).GT.alpha1.OR.a.EQ.0).AND.(move.EQ.1).AND.(DRp.EQ..true.)) THEN ! if we have rejected a nucleus move 
	
		AR1p(2)=AR1p(2)+1
		
		!write(*,*)'MOVE ount=',ount+ounti,'a=',a,'like_prop=',like_prop,'like=',like,'lgsigma=',lgsigma,'sum=',-like_prop+like+lgsigma,'log(rannum)=',log(rannum),'alpha1=',alpha1,'accept=',accept
	
!**********************************************************************************
!  If we reject a proposed position move, start delayed rejection on the position
!**********************************************************************************
	
		a2=0
		
		Voro_pprop=Voro
		IF (sigdep.EQ.2) srdist_prop = srdist

		cellogic=.false.
		
		p1 = Voro(1,ind)+(GASDEV(randcount)*pd2*(longmax-longmin)/100)
		p2 = Voro(2,ind)+(GASDEV(randcount)*pd2*(latmax-latmin)/100)
		
		IF ((p1.GT.longmin).and.(p2.GT.latmin).and.&		! if we are inside the acceptable range
			(p1.LT.longmax).and.(p2.LT.latmax)) THEN
		
			a2=1
			Voro_pprop(1,ind) = p1
			Voro_pprop(2,ind) = p2
			
		END IF
	
		IF(a2.EQ.1) THEN
		
			! Calculate triangulation of the previous (accepted) model (next 6 lines added by EG)
			CALL delaun (DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			
			CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
		
			! Update cellogic for the selected cell
			cellogic(ind)=.true.
			
			! Find the neighbours of the selected cell before moving - these cells are affected by the geometry change
			istart = nnn(ind)
			iend = nnn(ind+1)-1
			
			! Update cellogic for the neighbouring cells			
			DO k=istart,iend
				IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
			END DO
			
			! Compute the proposed Delaunay triangulation : nnn_prop,nnlist_prop
			CALL delaun (DBLE(Voro_pprop(1:2,:)),ncell_prop,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			IF (ptsok.EQ.1) THEN
				WRITE(*,*)'All',ncell_prop,'points in a line for sample',ount,'on processor',rank+1
				WRITE(*,*)'Sample number',ount,'discarded...',birth,death,move,velch,sigmach,ptsok
				!STOP
				ount=ount-1
				CYCLE
			END IF
			
			CALL build_nv(ncell_prop,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
			
			! Find the neighbours of the selected cell after moving - these cells are affected by the geometry change
			istart = nnn(ind)
			iend = nnn(ind+1)-1
			
			! Update cellogic for the neighbouring cells
			DO k=istart,iend
				IF (nnlist(k).ne.0) cellogic(nnlist(k))=.true. 
			END DO
			
			! Update RayCel_pprop just for rays concerned in the move
			RayCel_pprop = RayCel
			
			!!! ******************************************************************************************************
			!!! Original code modified by EG - option to update raypaths added
			
			IF(urps.EQ.2) THEN
			
				! Update raypaths by calling subroutine updaterays
				!write(*,*)'MOVE - DR'
				IF (uar.EQ.1) cellogic=.true.

				!t1_rays = 0
		  	        !t2_rays = 0
				!CALL cpu_time(t1_rays) ! TOC - Start counting time

				CALL updaterays(RayCel_pprop,RayCel,Voro_pprop,ncell,ncell_prop,cellogic,ind,rank+1,ount, t_uprays, t_uprays_test)

				!CALL cpu_time(t2_rays) ! TOC - Stop counting time 
				! Calculate cumulative time spent in modrays routine 
		  		t_rays_test_base =  t_rays_test_base  +  t_uprays
				t_rays_test =  t_rays_test  +  t_uprays_test


			
				! If there are crazy raypaths for this sample, cycle the loop and propose a new sample
				IF(crazyray.GT.0) THEN
					WRITE(*,*)'Sample number',ount,'discarded...'
					!STOP
					ount=ount-1
					CYCLE
				END IF
				
			ELSE 
			
				! Do not update the raypaths
				
				node=1
				
				DO nr=1,nrays
					DO i=1,ncell
						IF ((cellogic(i)).and.(RayCel(nr,i).ne.0)) THEN
							
							RayCel_pprop(nr,:) = 0
							
							jstart = nnM(nr)
							jend =  nnM(nr)+lengthrayv(nr)-2
							
							DO j=jstart,jend
								CALL find_node2D(DBLE(M(j,1:2)),node,DBLE(Voro_pprop(1:2,:)),nnn,nnlist,walk)
								RayCel_pprop(nr,node)=RayCel_pprop(nr,node)+ M(j,3)
							END DO
							
							EXIT
							
						END IF
					END DO
				END DO
			END IF
			
		END IF ! (a2=0) if we are outside the prior
	
		
		! Get like_pprop for the second try
		IF (a2.EQ.1) THEN ! compute the likelihood
			
			like_pprop=0
			
			DO nr=1,nrays
			
				ttime_pprop(nr)=0
				
				DO i=1,ncell_prop
					ttime_pprop(nr) = ttime_pprop(nr) + (RayCel_pprop(nr,i))/(Voro_pprop(3,i))
				END DO
				
				IF(sigdep.EQ.0) THEN
					IF (uglpd.EQ.0) THEN
						like_pprop = like_pprop + (ttime_pprop(nr)-dat(nr,5))**2/(2*(sigma_prop(1,dsstat(nr))**2))
					ELSE IF (uglpd.EQ.1) THEN
						like_pprop = like_pprop + ABS(dat(nr,5)-ttime_pprop(nr))/sigma_prop(1,dsstat(nr))
					END IF
				ELSE IF(sigdep.NE.0) THEN
					IF (uglpd.EQ.0) THEN
						like_pprop = like_pprop + (ttime_pprop(nr)-dat(nr,5))**2/(2*(sigma_prop(nr,1)**2))
					ELSE IF (uglpd.EQ.1) THEN
						like_pprop = like_pprop + ABS(dat(nr,5)-ttime_pprop(nr))/sigma_prop(nr,1)
					END IF
				END IF
				
			END DO
			
		ELSE ! (a2=0) if we are outside the prior
		
			like_pprop = like_prop 
			
		END IF

		!*********************************!
		! NOW NEED TO COMPILE THE ALPHA 2 !
		!*********************************!
		
		IF (nodata.EQ.0) THEN
			TT1=-like_pprop+like+lgsigma
		ELSE IF (nodata.EQ.1) THEN
			TT1=log(one)	! ignore the data
		END IF
		TT2a=((Voro_prop(1,ind)-Voro(1,ind))**2-(Voro_prop(1,ind)-Voro_pprop(1,ind))**2)/(2*pd**2)
		TT2b=((Voro_prop(2,ind)-Voro(2,ind))**2-(Voro_prop(2,ind)-Voro_pprop(2,ind))**2)/(2*pd**2)
		
		alpha2=0
		rannum2=grnd()
		randcount=randcount+1
		
		IF (a2.EQ.1.AND.a.EQ.0) THEN		! if the first proposal was rejected because out of bounds, then TT3=0
			TT3t=0
			TT3b=0
			TT3=0
			alpha2=MINVAL([log(one),TT1+TT2a+TT2b])
		ELSE IF (a2.EQ.1.AND.a.EQ.1) THEN	! if both the first and the delayed proposal are within the bounds
			IF (nodata.EQ.0) THEN
				TT3t=log(1-exp(MINVAL([log(one),-like_prop+like_pprop+lgsigma])))
			ELSE IF (nodata.EQ.1) THEN	! shouldn't even get to this point because for no data all m' are accepted
				TT3t=log(1-exp(log(one)))	! ignore the data
			END IF
			TT3b=log(1-exp(alpha1))
			TT3=TT3t-TT3b
			alpha2=MINVAL([log(one),TT1+TT2a+TT2b+TT3])
		END IF
		
		! Now compare
		IF (a2.EQ.1.AND.log(rannum2).LE.alpha2) THEN
			
			! Accept the second attempt
			accept=.true.
			
			AR2p(1)=AR2p(1)+1
			
			like_prop = like_pprop
			Voro_prop = Voro_pprop
			ttime_prop = ttime_pprop
			RayCel_prop = RayCel_pprop
			
		ELSE
			AR2p(2)=AR2p(2)+1
		END IF
		
		!write(*,*)'MOVE2 ount=',ount+ounti,'a=',a,'a2=',a2,'like_pprop=',like_pprop,'TT1=',TT1,'TT2a=',TT2a,'TT2b=',TT2b,'TT3t=',TT3t,'TT3b=',TT3b,'TT3=',TT3,'log(rannum2)=',log(rannum2),'alpha2=',alpha2,'accept=',accept
		
	END IF
!*********************END DELAYED REJECTION FOR THE POSITION *************************
!*************************************************************************************

ELSE IF (sigmach.EQ.1) THEN

	IF (a.EQ.1) THEN	! if we are inside the prior evaluate alpha1, otherwise don't even bother
		IF (nodata.EQ.0) THEN
			alpha1=MINVAL([log(one),-like_prop+like+lgsigma])
		ELSE IF (nodata.EQ.1) THEN
			alpha1=log(one)	! here the likelihood is set to unity (we remove the data)
		END IF
	END IF

	! If we are inside the prior and rand<alpha1 ACCEPT the model
	IF ((a.EQ.1).AND.(log(rannum).LE.alpha1)) accept=.true.
	
	!write(*,*)'SIGMACH ount=',ount+ounti,'a=',a,'like_prop=',like_prop,'like=',like,'lgsigma=',lgsigma,'sum=',-like_prop+like+lgsigma,'log(rannum)=',log(rannum),'alpha1=',alpha1,'accept=',accept

END IF


!*************************************************************************************
!   Output proposed model depending on choice made in input file
!*************************************************************************************
IF(ops.EQ.-1.OR.rank+1.EQ.ops) THEN
	IF(mod(ount,opsi).EQ.0) THEN
	
		WRITE(*,*)'Saving proposed model on processor',rank+1
		
		! Plot the sample map
		WRITE(ountc,*)ount+ounti
		ountc = ADJUSTL(ountc)
		sampleout = ADJUSTL(sampleout)
		IF(accept.EQ..true.) THEN
			samplevtx = trim(sampleout) // '_' // trim(ountc) // '_p' // trim(rankc) // '_acc.vtx'
			samplevor = trim(sampleout) // '_' // trim(ountc) // '_p' // trim(rankc) // '_acc.vor'
		ELSE IF(accept.EQ..false.) THEN
			samplevtx = trim(sampleout) // '_' // trim(ountc) // '_p' // trim(rankc) // '_rej.vtx'
			samplevor = trim(sampleout) // '_' // trim(ountc) // '_p' // trim(rankc) // '_rej.vor'
		END IF
		
		IF(test.EQ.1) THEN
			samplevtx = trim(sampleout) // trim(ountc) // '.vtx'
			samplevor = trim(sampleout) // trim(ountc) // '.vor'
			IF(urps.eq.2) THEN
				string = 'cp rayupdate.dat raysub' // trim(ountc) // '.out'
				CALL SYSTEM( string )
			END IF
		END IF
		
		! Calculate Delaunay triangulation
		CALL delaun(DBLE(Voro_prop(1:2,:)),ncell_prop,neighbour,vertices,nt,nt_max,& 	! outputs are neighbour, vertices, nt (number of triangles)
			worki1,worki2,worki3,eps,nv_max,&
			0,ldummy,0,0,0,ptsok)
		CALL build_nv(ncell_prop,vertices,nt,ncell_max,nmax,&
              			neighbour,nnn,nnlist,ntwork)
				
		OPEN(UNIT=50,FILE=samplevtx,STATUS='replace')
		WRITE(50,*)nvtr,nvpr
		WRITE(50,'(3f14.8)')latmax,longmin
		WRITE(50,'(3f14.8)')latr,longr
		WRITE(50,'(1X)')
	
		OPEN(UNIT=51,FILE=samplevor,STATUS='replace')
		WRITE(51,*)nvtr,nvpr
		WRITE(51,'(3f14.8)')latmax,longmin
		WRITE(51,'(3f14.8)')latr,longr
		WRITE(51,'(1X)')
		WRITE(51,*)ncell_prop
	
		DO i=0,nvpr+1
			DO j=0,nvtr+1
			
				ii=i
				jj=j
				IF (i.EQ.0)   ii=1
				IF (j.EQ.0)   jj=1
				IF (i.EQ.nvpr+1) ii=nvpr
				IF (j.EQ.nvtr+1) jj=nvtr
			
				p1=longmin+(ii-1)*longr
				p2=latmax-(jj-1)*latr
				
				node=1
				
				CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro_prop(1:2,:)),nnn,nnlist,walk)
			
				WRITE(50,*)Voro_prop(3,node)
		
			END DO
			WRITE(50,'(1X)')		
		END DO
		
		CLOSE(50) ! close the file
	
		DO c=1,ncell_prop
			WRITE(51,*)Voro_prop(2,c),Voro_prop(1,c),Voro_prop(3,c)
		END DO
		
		CLOSE(51) ! close the file
		
	END IF
END IF


!*************************************************************************************
!   If we accept the proposed model, update the status of the Markov Chain
!   and update the acceptance ratios for each type of move. 
!   Depending on the choice made in the input file, update (or not) the raypaths.
!*************************************************************************************
IF (accept.EQ..true.) THEN

	! Number of accepted models
	totacc=totacc+1

	! Acceptance ratio
	AR(1)=AR(1)+1
	
	! Update the chain
	like = like_prop
	Voro = Voro_prop
	ttime = ttime_prop
	sigma = sigma_prop
	RayCel = RayCel_prop
	IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
		aa = aa_prop
		bb = bb_prop
	ELSE IF (sigdep.EQ.3) THEN
		lambda = lambda_prop
	END IF
   
	IF (birth.EQ.1.OR.death.EQ.1.OR.move.EQ.1) THEN
		ncell = ncell_prop
		IF (birth.EQ.1) ARB(1)=ARB(1)+1 
		IF (death.EQ.1) ARD(1)=ARD(1)+1 
		IF (move.EQ.1)  ARM(1)=ARM(1)+1
		IF (sigdep.EQ.2) srdist = srdist_prop
	ELSE IF (velch.EQ.1) THEN
		ARV(1)=ARV(1)+1
		IF (sigdep.EQ.2) srdist = srdist_prop
	ELSE IF (sigmach.EQ.1) THEN
		ARS(1)=ARS(1)+1
		IF ((sigdep.EQ.1.OR.sigdep.EQ.2).AND.cha.EQ.1) ARSa(1)=ARSa(1)+1
		IF ((sigdep.EQ.1.OR.sigdep.EQ.2).AND.chb.EQ.1) ARSb(1)=ARSb(1)+1
	END IF
	
	! Calculate misfit of accepted model
	misfit = 0
	DO nr=1,nrays
		IF(sigdep.EQ.0) THEN
			misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(1,dsstat(nr)))**2
		ELSE IF(sigdep.NE.0) THEN
			misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(nr,1))**2
		END IF
	END DO 

	! Last accepted model
	nlastacc=ncell
	lastacc(:,1:nlastacc)=Voro(:,1:ncell)   
   
ELSE
	AR(2)=AR(2)+1
        IF (birth.EQ.1.OR.death.EQ.1.OR.move.EQ.1) THEN
		IF (birth.EQ.1) ARB(2)=ARB(2)+1 
		IF (move.EQ.1)  ARM(2)=ARM(2)+1
		IF (death.EQ.1)  ARD(2)=ARD(2)+1
	ELSE IF (velch.EQ.1) THEN
		ARV(2)=ARV(2)+1
	ELSE IF (sigmach.EQ.1) THEN
		ARS(2)=ARS(2)+1
		IF ((sigdep.EQ.1.OR.sigdep.EQ.2).AND.cha.EQ.1) ARSa(2)=ARSa(2)+1
		IF ((sigdep.EQ.1.OR.sigdep.EQ.2).AND.chb.EQ.1) ARSb(2)=ARSb(2)+1
	END IF
END IF

IF ((AR(2)+AR(1)).GT.0) AR(3)=100*AR(1)/(AR(2)+AR(1))
IF ((ARB(2)+ARB(1)).GT.0) ARB(3)=100*ARB(1)/(ARB(2)+ARB(1))
IF ((ARD(2)+ARD(1)).GT.0) ARD(3)=100*ARD(1)/(ARD(2)+ARD(1))
IF ((ARM(2)+ARM(1)).GT.0) ARM(3)=100*ARM(1)/(ARM(2)+ARM(1))
IF ((ARV(2)+ARV(1)).GT.0) ARV(3)=100*ARV(1)/(ARV(2)+ARV(1))
IF ((ARS(2)+ARS(1)).GT.0) ARS(3)=100*ARS(1)/(ARS(2)+ARS(1))
IF ((ARSa(2)+ARSa(1)).GT.0) ARSa(3)=100*ARSa(1)/(ARSa(2)+ARSa(1))
IF ((ARSb(2)+ARSb(1)).GT.0) ARSb(3)=100*ARSb(1)/(ARSb(2)+ARSb(1))
IF ((AR1p(2)+AR1p(1)).GT.0) AR1p(3)=100*AR1p(1)/(AR1p(2)+AR1p(1))
IF ((AR2p(2)+AR2p(1)).GT.0) AR2p(3)=100*AR2p(1)/(AR2p(2)+AR2p(1))
IF ((AR1v(2)+AR1v(1)).GT.0) AR1v(3)=100*AR1v(1)/(AR1v(2)+AR1v(1))
IF ((AR2v(2)+AR2v(1)).GT.0) AR2v(3)=100*AR2v(1)/(AR2v(2)+AR2v(1))

! Save accepted sample to matrix
! The various models are saved in matrices samplesStep, samplesInd and samplesVal:
! - samplesStep contains the step type in column 1 (birth=1, death=2, move=3, velocity=4, sigma=5)
!               and the acceptance switch in column 2 (accept=1, reject=0)
! - samplesInd  contains the index of the modified cell for a death, move or velocity step
! - samplesVal  contains new longitude, latitude and velocity values for the affected cell
! The values are saved to these 3 matrices only if the proposed model is accepted, otherwise only step type
! and acceptance switch are saved (the other values remain =0)
IF (runnum.EQ.1) THEN
! If this is the first run of the chain, then the initial model is also saved and all folowing models are at position ount+1
	IF (birth.EQ.1) THEN
		samplesStep(ount+1,1)=1
		IF (accept.EQ..true.) THEN
			samplesStep(ount+1,2)=1
			samplesInd(ount+1)=0
			samplesVal(ount+1,1)=Voro(1,ncell)
			samplesVal(ount+1,2)=Voro(2,ncell)
			samplesVal(ount+1,3)=Voro(3,ncell)
		END IF
	ELSE IF (death.EQ.1) THEN
		samplesStep(ount+1,1)=2
		IF (accept.EQ..true.) THEN
			samplesStep(ount+1,2)=1
			samplesInd(ount+1)=ind
			samplesVal(ount+1,1)=0
			samplesVal(ount+1,2)=0
			samplesVal(ount+1,3)=0
		END IF
	ELSE IF (move.EQ.1) THEN
		samplesStep(ount+1,1)=3
		IF (accept.EQ..true.) THEN
			samplesStep(ount+1,2)=1
			samplesInd(ount+1)=ind
			samplesVal(ount+1,1)=Voro(1,ind)
			samplesVal(ount+1,2)=Voro(2,ind)
			samplesVal(ount+1,3)=0
		END IF
	ELSE IF (velch.EQ.1) THEN
		samplesStep(ount+1,1)=4
		IF (accept.EQ..true.) THEN
			samplesStep(ount+1,2)=1
			samplesInd(ount+1)=ind
			samplesVal(ount+1,1)=0
			samplesVal(ount+1,2)=0
			samplesVal(ount+1,3)=Voro(3,ind)
		END IF
	ELSE IF (sigmach.EQ.1) THEN
		samplesStep(ount+1,1)=5
		IF (accept.EQ..true.) THEN
			samplesStep(ount+1,2)=1
			samplesInd(ount+1)=0
			samplesVal(ount+1,1)=0
			samplesVal(ount+1,2)=0
			samplesVal(ount+1,3)=0
		END IF
	END IF
	samplesNcel(ount+1)=ncell
	DO nss=1,nds
		IF (sigdep.EQ.0) THEN
			samplesSigmas(ount+1,nss)=sigma(1,nss)
		ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
			samplesSigmas(ount+1,nss)=aa(nss)
			samplesSigmas(ount+1,nss+nds)=bb(nss)
		ELSE IF (sigdep.EQ.3) THEN
			samplesSigmas(ount+1,nss)=lambda(nss)
		END IF
	END DO
	samplesMisfit(ount+1) = misfit
	samplesAratios(ount+1,1) = AR(3)
	samplesAratios(ount+1,2) = ARB(3)
	samplesAratios(ount+1,3) = ARD(3)
	samplesAratios(ount+1,4) = ARM(3)
	samplesAratios(ount+1,5) = ARV(3)
	samplesAratios(ount+1,6) = ARS(3)
ELSE IF (runnum.GT.1) THEN
! If this is not the first run of the chain, each model is saved at its position ount
	IF (birth.EQ.1) THEN
		samplesStep(ount,1)=1
		IF (accept.EQ..true.) THEN
			samplesStep(ount,2)=1
			samplesInd(ount)=0
			samplesVal(ount,1)=Voro(1,ncell)
			samplesVal(ount,2)=Voro(2,ncell)
			samplesVal(ount,3)=Voro(3,ncell)
		END IF
	ELSE IF (death.EQ.1) THEN
		samplesStep(ount,1)=2
		IF (accept.EQ..true.) THEN
			samplesStep(ount,2)=1
			samplesInd(ount)=ind
			samplesVal(ount,1)=0
			samplesVal(ount,2)=0
			samplesVal(ount,3)=0
		END IF
	ELSE IF (move.EQ.1) THEN
		samplesStep(ount,1)=3
		IF (accept.EQ..true.) THEN
			samplesStep(ount,2)=1
			samplesInd(ount)=ind
			samplesVal(ount,1)=Voro(1,ind)
			samplesVal(ount,2)=Voro(2,ind)
			samplesVal(ount,3)=0
		END IF
	ELSE IF (velch.EQ.1) THEN
		samplesStep(ount,1)=4
		IF (accept.EQ..true.) THEN
			samplesStep(ount,2)=1
			samplesInd(ount)=ind
			samplesVal(ount,1)=0
			samplesVal(ount,2)=0
			samplesVal(ount,3)=Voro(3,ind)
		END IF
	ELSE IF (sigmach.EQ.1) THEN
		samplesStep(ount,1)=5
		IF (accept.EQ..true.) THEN
			samplesStep(ount,2)=1
			samplesInd(ount)=0
			samplesVal(ount,1)=0
			samplesVal(ount,2)=0
			samplesVal(ount,3)=0
		END IF
	END IF
	samplesNcel(ount)=ncell
	DO nss=1,nds
		IF (sigdep.EQ.0) THEN
			samplesSigmas(ount,nss)=sigma(1,nss)
		ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
			samplesSigmas(ount,nss)=aa(nss)
			samplesSigmas(ount,nss+nds)=bb(nss)
		ELSE IF (sigdep.EQ.3) THEN
			samplesSigmas(ount,nss)=lambda(nss)
		END IF
	END DO
	samplesMisfit(ount) = misfit
	samplesAratios(ount,1) = AR(3)
	samplesAratios(ount,2) = ARB(3)
	samplesAratios(ount,3) = ARD(3)
	samplesAratios(ount,4) = ARM(3)
	samplesAratios(ount,5) = ARV(3)
	samplesAratios(ount,6) = ARS(3)
END IF

!**********************************************************************
!       Display what is going on every "Display" samples
!**********************************************************************
IF (mod(ount,display).EQ.0) THEN

	misfit=0
	DO nr=1,nrays
		IF(sigdep.EQ.0) THEN
			misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(1,dsstat(nr)))**2
		ELSE IF(sigdep.NE.0) THEN
			misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(nr,1))**2
		END IF
	END DO

        WRITE(*,*)
	WRITE(*,*)'Processor number',rank+1,'/',nbproc
	WRITE(*,*)'Sample:',ount+ounti,'/',sampletotal
	WRITE(*,*)'Number of cells:',ncell
	WRITE(*,*)'Average misfit:',misfit/nrays
	IF (sigdep.EQ.0) THEN
		DO nss=1,nds
			WRITE(*,*)'Data noise for dataset',nss,':',sigma(1,nss) 
		END DO
	ELSE IF (sigdep.NE.0) THEN
		WRITE(*,*)'Average data noise:',SUM(sigma(:,1))/nrays
		IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
			DO nss=1,nds
				WRITE(*,*)'Parameter "a" for dataset',nss,':',aa(nss)
				WRITE(*,*)'Parameter "b" for dataset',nss,':',bb(nss)
			END DO
		ELSE IF (sigdep.EQ.3) THEN
			DO nss=1,nds
				WRITE(*,*)'Parameter "lambda" for dataset',nss,':',lambda(nss)
			END DO
		END IF
	END IF
        WRITE(*,*)
	WRITE(*,*)'Acceptance rate for velocity'
	WRITE(*,*)'AR1v',AR1v(3),'AR2v',AR2v(3),'ARV:',ARV(3)
	WRITE(*,*)'Acceptance rate for position'
	WRITE(*,*)'AR1p',AR1p(3),'AR2p',AR2p(3),'ARM:',ARM(3)
	WRITE(*,*)'Acceptance rate for birth and death'
	WRITE(*,*)'ARB:',ARB(3),'ARD:',ARD(3)
	WRITE(*,*)'Acceptance rate for sigma'
	WRITE(*,*)'ARS:',ARS(3)
	IF (sigdep.EQ.1.OR.sigdep.EQ.2) WRITE(*,*)'ARSa:',ARSa(3),'ARSb:',ARSb(3)
	WRITE(*,*)'Total acceptance rate',AR(3)
	WRITE(*,*)'________________________________________________________________________'
	WRITE(*,*)
END IF

!!! *****************************************************************************************
!!! Modification by EG - update raypaths
IF(urps.EQ.1) THEN
	
	node=1
	
	! Calculate Delaunay triangulation
	CALL delaun(DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
                     		worki1,worki2,worki3,eps,nv_max,&
                      		0,ldummy,0,0,0,ptsok)
	CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
        			neighbour,nnn,nnlist,ntwork)
			
	iii=0
	DO i=0,nvpr+1
		jjj=0
		iii=iii+1
		DO j=0,nvtr+1
			jjj=jjj+1					
			ii=i
			jj=j
			IF (i.EQ.0)   	ii=1
			IF (j.EQ.0)   	jj=1
			IF (i.EQ.nvpr+1) 	ii=nvpr
			IF (j.EQ.nvtr+1) 	jj=nvtr
				
			p1=longmin+(ii-1)*longr
			p2=latmax-(jj-1)*latr
			
			node=1
			
			CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro(1:2,:)),nnn,nnlist,walk)
			
			VELU(iii,jjj) = VELU(iii,jjj)+(Voro(3,node))
			
			END DO
	END DO
	
	IF(mod(ount,upthin).EQ.0) THEN
   
		WRITE(*,*)'Updating raypaths at sample',ount+ounti,'of',sampletotal,'on rank',rank+1
		
		IF(ALLOCATED(npoints)) DEALLOCATE(npoints)
		IF(ALLOCATED(raycoords)) DEALLOCATE(raycoords)
		IF(ALLOCATED(raypoints)) DEALLOCATE(raypoints)
		IF(ALLOCATED(M)) DEALLOCATE(M)
		IF(ALLOCATED(nnM)) DEALLOCATE(nnM)
	
		IF (rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
			! Save the updated geometry to file
			updout = ADJUSTL(updout)
			updvtx = trim(updout) // '_p' // trim(rankc) // '.vtx'
			!updvor = trim(updout) // '_p' // trim(rankc) // '.vor'
			
			IF(test.EQ.1) THEN
				updvtx = trim(updout) // '.vtx'
				!updvor = trim(updout) // '.vor'
			END IF	
			
			OPEN(UNIT=62,FILE=updvtx,STATUS='replace')
			WRITE(62,*)nvtr,nvpr
			WRITE(62,'(3f14.8)')latmax,longmin
			WRITE(62,'(3f14.8)')latr,longr
			WRITE(62,'(1X)')
			
			!OPEN(UNIT=64,FILE=updvor,STATUS='replace')
			!WRITE(64,*)nvtr,nvpr
			!WRITE(64,'(3f14.8)')latmax,longmin
			!WRITE(64,'(3f14.8)')latr,longr
			!WRITE(64,'(1X)')
			!WRITE(64,*)ncell		
		END IF		
		
		iii=0
		DO i=0,nvpr+1
			jjj=0
			iii=iii+1
			DO j=0,nvtr+1
				jjj=jjj+1					
				
				IF (rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
					WRITE(62,*)VELU(iii,jjj)/upthin
				END IF

			END DO

			IF (rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
				WRITE(62,'(1X)')
			END IF
		
		END DO
		
		IF (rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
			CLOSE(62)
		END IF
		
		!DO c=1,ncell
		!
		!	IF (rank.EQ.0) THEN				
		!		WRITE(64,*)Voro(2,c),Voro(1,c),Voro(3,c)
		!	END IF
		!
		!END DO
		
		!IF (rank.EQ.0) THEN
		!	CLOSE(64) ! close the file
		!END IF
		
		nru=sum(raystat(:,1),1)
		
		IF(nru.NE.nrays) THEN
			IF(rank.EQ.0) WRITE(*,*)'There is a problem with the number of active rays!!'
			CALL MPI_FINALIZE(ierror) !! comment for former MPI parallelization
			STOP
		END IF
						
		WRITE(*,*)'Calling MODRAYS to model raypaths for the accepted model on rank',rank+1


		!!------------------------------------------------------------------------------------------------------------------------------
		!! This is the last stage of the MPI parallelization 

		IF(TEST_MODRAYS.EQ.1) THEN 
			t1_rays = 0
	  	        t2_rays = 0
			CALL cpu_time(t1_rays) ! TOC - Start counting time
					
			CALL modrays(nsou,sou(:,1),sou(:,2),&
					nrec,rec(:,1),rec(:,2),&
					raystat,-1,&
					nvtr,nvpr,latmax,longmin,latr,longr,VELU/upthin,&
					gdt,gdp,&
					sgref,&
					dl,erg,&
					er,&
					fms,&
					nbs,&
					rank+1)


			CALL cpu_time(t2_rays) ! TOC - Stop counting time 
			! Calculate cumulative time spent in modrays routine 
	  		t_rays_test_base =  t_rays_test_base  +  (t2_rays-t1_rays)
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Inner parallelization starts here.  
		! We start by sending the source sets to the workers for the modrays routine for calculation.  We will send the partial  
		! information required by the modrays routines to every worker.      	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF(ALLOCATED(npoints)) DEALLOCATE(npoints)
		IF(ALLOCATED(raycoords)) DEALLOCATE(raycoords)
		IF(ALLOCATED(raypoints)) DEALLOCATE(raypoints)
		IF(ALLOCATED(Mi)) DEALLOCATE(Mi)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Time counter

		t1_rays_test = 0
 	        t2_rays_test = 0
		CALL cpu_time(t1_rays_test) ! TOC - Start counting time

		!t1_rays = 0
  	        !t2_rays = 0
		!CALL cpu_time(t1_rays) ! TOC - Start counting time
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		WRITE(*,*) 'Start parallelization'

		crazyray = 0 
		IF(ALLOCATED(offset_master)) DEALLOCATE(offset_master)
		IF(ALLOCATED(nsou_master)) DEALLOCATE(nsou_master)

		ALLOCATE(offset_master(numworkers,2))
		ALLOCATE(nsou_master(numworkers))

		! (We store data in several matrices and send to the worker tasks).  

		!WRITE(*,*) 'Number of sources', nsou
		avesou = nsou/numworkers  ! Number of sources in each subset, i.e. per worker excluding remainder.   e.g. 11/3 = 3. remainder 2.  
		extrasou = MOD(nsou,numworkers)		

                ! The offsets are the indices used for transiting along the arrays
		! In the fisrt loop we calculate the offsets that will be used by the workers.   
		! In the fisrt loop we calculate the nsources that will be send to each worker.   
  
                offset = 0        

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               DO dest=1,numworkers  
		    ! Here we need to count the number of rays corresponding to the sources sent to the workers.  Number of sources sent is nsou_sent.             
                    IF(dest.LE.extrasou) THEN     														
		        nsou_sub = avesou+1  	! nsou_subset are the number of sources sent to the worker for calculating the modrays algorithm.       		
                    ELSE 
                        nsou_sub = avesou     			
	            END IF
                    offset_end = offset + nsou_sub   ! This defines the interval we are sending.  Ex:  offset = 0, nsou_subset = 3, offset_end = 3.                      
		    nsou_master(dest) = nsou_sub
		    offset_master(dest, 1) = offset
 		    offset_master(dest, 2) = offset_end
		    ! ATTENTION offsetpts and offsetcoords are calculated and looped in the second master processor loop, but after receiving the npoints and raycoords results from the worker.
                    offset = offset_end	
  	        END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! TAG FROM_MASTER = 1              
                ! TO: DEST 
		! Send to worker ===> (dest,1)

                !Brodcast

		!nvx,nvz,goxd,gozd,dvxd,dvzd,velv
		params1(1) = nvtr
		params1(2) = nvpr
		params1(3) = latmax
		params1(4) = longmin
		params1(5) = latr
		params1(6) = longr
		params1(7) = er ! Earth radius
		params1(8) = nbs

                wrgf = -1 	                      
		params2(1) = nsou
		params2(2) = nrec
		params2(3) = wrgf          	! wrgf = write geometries to file?
		params2(4) = gdt		! Grid dicing
		params2(5) = gdp		! Grid dicing
		params2(6) = sgref		! Apply source grid refinement? (0=no,1=yes)
		params2(7) = dl			! Dicing level and extent of refined grid
		params2(8) = erg
		params2(9) = fms		! use first-order(0) or mixed-order(1) scheme
		params2(10)= rank+1		! Processor rank (in this case for embarrasingly parallel processing).
		params2(11)= ount+1		! Sample number! 


  	        counter_1 = 2*nsou*nrec               
                counter_2 = (nvpr+2)*(nvtr+2)
		!CALL MPI_BCAST(params1, 8, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 		! Brodcast parameters 1 to all processes 

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               DO dest=1,numworkers  
		    CALL MPI_SEND(0, 1, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)		 		! Send run_worker flag for running worker
		    CALL MPI_SEND(params1, 8, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror) 		! Brodcast parameters 1 to all processes 
    	            CALL MPI_SEND(params2, 11, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror) 			! Brodcast parameters 2 to all processes 
	            CALL MPI_SEND(offset_master, 2*numworkers, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)	
                    CALL MPI_SEND(nsou_master, numworkers, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)
		    CALL MPI_SEND(rec, nrec*2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)
		    CALL MPI_SEND(sou, nsou*2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)		   		   
		    CALL MPI_SEND(raystat, counter_1, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)
		    CALL MPI_SEND(VELU/upthin, counter_2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)
               END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
                !WRITE(*,*) 'All variables and arrays have been broadcasted to the workers'
		!CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

	        ! TAG = FROM_WORKER
                ! FROM: Worker
                ! Receive in master from worker ===> (from_worker,2)		

                offsetpts 	= 0  
                offsetcoords 	= 0  

		! Now we loop over all workers receiving the information for building the complete arrays.            

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	        DO  i=1,numworkers
		    ! ATTENTION offsetpts and offsetcoords are calculated and looped here, but are updated receiving the nru, npoints and raycoords results from the worker.
                    IF(ALLOCATED(offset_master)) DEALLOCATE(offset_master)
                    IF(ALLOCATED(npoints_master)) DEALLOCATE(npoints_master)
  		    ALLOCATE( offset_master(numworkers,2) )  		    
		    from_worker = i		
 
		    !TAG pair:  (i,2)		
                    CALL MPI_RECV(offset_master, 2*numworkers, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)			    		   
                    CALL MPI_RECV(nru_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)	                    
		    ALLOCATE( npoints_master( nru_master , 2 ) )			! array of length nru_subset, filled w. number of points that each ray contains and ray ID corresponding to the subset chosen.
		    CALL MPI_RECV(npoints_master, 2*nru_master, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)  

		    !numraycoords_master = sum( npoints_master( 1:nru_master , 1 ) )	
    	            CALL MPI_RECV(numraycoords_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)  
                  
                    IF(ALLOCATED(raycoords_master)) DEALLOCATE(raycoords_master)
                    ALLOCATE( raycoords_master( numraycoords_master, 2  )  )		! allocate array of ray coordinates for the sources subset.                   
		    CALL MPI_RECV(raycoords_master, numraycoords_master*2, MPI_DOUBLE_PRECISION, from_worker, 2, MPI_COMM_WORLD, status, ierror)
		    CALL MPI_RECV(crazyray_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)

		    offsetpts_end = offsetpts + nru_master
                    offsetcoords_end = offsetcoords + numraycoords_master

		    ! Next we fill in the raycoords array

		    IF (i.EQ.1) THEN		                             		    
                       ALLOCATE( npoints(nru_master, 2 ) )			! npoints is allocated here for the first time
		       npoints = npoints_master					! npoints is filled here for the first time
                       ALLOCATE( raycoords( numraycoords_master, 2  )  )      	! allocate array of ray coordinates for the sources subset for the first time
                       raycoords = raycoords_master				! raycoords is filled here for the first time				
                       IF(ALLOCATED(npoints_prev))DEALLOCATE( npoints_prev )  
                       IF(ALLOCATED(raycoords_prev))DEALLOCATE( raycoords_prev )  
		       ALLOCATE( npoints_prev( nru_master , 2 )  )
		       ALLOCATE( raycoords_prev( numraycoords_master , 2 ) )
		       npoints_prev   = npoints					  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.
		       raycoords_prev = raycoords				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.

		    ELSE IF (i.GT.1) THEN
                       DEALLOCATE( raycoords ) 
                       DEALLOCATE( npoints ) 
                       ALLOCATE( npoints( offsetpts_end, 2  )  )  	 	   ! allocate array of ray coordinates for the sources subset.		
                       ALLOCATE( raycoords( offsetcoords_end, 2  )  )      	   ! allocate array of ray coordinates for the sources subset.		
                       npoints(1:offsetpts, 1:2) = npoints_prev			   ! npoints is partially filled here with the new npoints subset that has been received recently.	                          
		       npoints(offsetpts+1: offsetpts_end, 1:2) = npoints_master   ! npoints is partially filled here with the new npoints subset that has been received recently.	                          
 		       raycoords(1:offsetcoords, 1:2 ) = raycoords_prev
                       raycoords(offsetcoords+1: offsetcoords_end, 1:2) = raycoords_master ! raycoords is partially filled here with the new raycoords subset that has been received recently.	                          

                       IF(ALLOCATED(npoints_prev))DEALLOCATE( npoints_prev )  
                       IF(ALLOCATED(raycoords_prev))DEALLOCATE( raycoords_prev )
	    	       ALLOCATE( npoints_prev( offsetpts_end , 2 )  )
		       ALLOCATE( raycoords_prev( offsetcoords_end , 2 ) )
		       raycoords_prev = raycoords				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.
		       npoints_prev   = npoints				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.

		    END IF 

                    offsetpts    = offsetpts_end 
                    offsetcoords = offsetcoords_end  !OK ALL REVISED UP TO THIS POINT.
     		    numraycoords = sum( npoints( : , 1 ) )                    ! size of raycoords_subset found by adding the first column of the corresponding subset of npoints = number of total ray coordinates in the subset.			
                    crazyray = crazyray + crazyray_master

		    ! Printing
		    !WRITE(*,*)'offset after 26 equals', offset_master(1,1), 'in processor = ',i
		    !WRITE(*,*)'offset after 26 equals', offset_master(3,1), 'in processor = ',i
		    !WRITE(*,*)'offset after 26 equals', offset_master(5,1), 'in processor = ',i
		    !WRITE(*,*)'offset after 26 equals', offset_master(7,1), 'in processor = ',i	
		    !WRITE(*,*)'nru_master is ', nru_master, 'in processor = ',i
		    !WRITE(*,*)'npoints_master is ', npoints_master(1,1), 'in processor = ',i
		    !WRITE(*,*)'numraycoords_master is ', numraycoords_master, 'in processor = ',i
		    !WRITE(*,*)'raycoords_master is ', raycoords_master(1,1), 'in processor = ',i
		    !WRITE(*,*)'crazyray_master is ', crazyray_master, 'in processor = ',i
   		    !WRITE(*,*)'raycoords is ', crazyray_master, 'in processor = ',i
		    !WRITE(*,*)'crazyray_master is ', crazyray_master, 'in processor = ',i
		    !WRITE(*,*)'Partial numraycoords at master is', numraycoords, 'processor', taskid
		    !WRITE(*,*)'Partial size of npoints at master is ', SHAPE(npoints)
		    !WRITE(*,*)'Partial size of raycoords', SHAPE(raycoords)
	        END DO

		    CALL cpu_time(t2_rays_test) ! TOC - Stop counting time 
		    ! Calculate cumulative time spent in modrays routine when parallelized
  		    t_rays_test =  t_rays_test  +  (t2_rays_test-t1_rays_test)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Inner parallelization ends here.  
		! We obtain as a result:
		! npoints, raycoords, crazyray.     	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		IF(TEST_MODRAYS.EQ.1) THEN 
  		    !Here finally we compare the results of both algorithms 
		    !WRITE(*,*)'crazyray_master is ', crazyray_master
		    !WRITE(*,*)'nru_master is ', offsetpts_end
		    !WRITE(*,*)'Size of (parallel) npoints at master is ', SHAPE(npoints)
		    !WRITE(*,*)'Size of (parallel) raycoords at master is', SHAPE(raycoords)
		    !WRITE(*,*)'original crazyray is ', crazyray_base
		    !WRITE(*,*)'original nru is ', nru_base
		    !WRITE(*,*)'Size of (original) npoints is ', SHAPE(npoints_base)
		    !WRITE(*,*)'Size of (original) raycoords is', SHAPE(raycoords_base)
		    !WRITE(*,*)

		    npoints_base = npoints_base - npoints
		    raycoords_base = raycoords_base - raycoords 
		    npoints_diff1 = abs( sum( npoints_base(:, 1) ) ) 
                    npoints_diff2 = abs( sum( npoints_base(:, 2) ) )
		    raycoords_diffx = abs( sum( raycoords_base(:, 1) ) )
		    raycoords_diffy = abs( sum( raycoords_base(:, 2) ) )

   		    WRITE(*,*) 'Results parallelization at modrays (outside), for sample', ount
		    WRITE(*,*) 'Sum of the diff sample', ount ,'is (', npoints_diff1,',', npoints_diff2, ')' 
		    WRITE(*,*) 'She sum of the diff sample', ount ,'is (', raycoords_diffx,',', raycoords_diffy, ')' 
		END IF            

		!!-------------------------------------------------------------------------------------------------------------------------------
				
		ALLOCATE(nnM(nru))
		ALLOCATE(M(sum(npoints(:,1),1)-nru,3))
		
		IF(wurf.EQ.-1.OR.wurf.EQ.rank+1.OR.rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
			rayupdout = 'rayupdate_p' // trim(rankc) // '.dat'
			IF(test.EQ.1) THEN
				rayupdout = 'rayupdate.dat'
			END IF		
			OPEN(UNIT=40,FILE=rayupdout,FORM='unformatted',STATUS='unknown')
			WRITE(40)nru
		END IF			
	
		pti=0
		ptf=0
		p=0
		DO nr=1,nru
		
			rayid=npoints(nr,2)
			raylength=npoints(nr,1)
		
			pti=pti+1
			ptf=pti+raylength-1
			
			ALLOCATE(raypoints(raylength,2))
			raypoints=raycoords(pti:ptf,:)
			
			pti=ptf
			
			lengthrayv(rayid)=raylength
			
			angle=0
				
			point(2) = raypoints(1,1)
			point(1) = raypoints(1,2)
			point(1) = point(1)*pii/180
			point(2) = point(2)*pii/180
					
			RayCel(rayid,:) = 0 
			
			nnM(nr)=p+1
			
			IF(wurf.EQ.-1.OR.wurf.EQ.rank+1.OR.rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
				WRITE(40)raylength
				WRITE(40)raypoints(1,1),raypoints(1,2)
			END IF
				
			DO i=2,raylength
			
				IF(wurf.EQ.-1.OR.wurf.EQ.rank+1.OR.rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
					WRITE(40)raypoints(i,1),raypoints(i,2)
				END IF
			
				p=p+1
				M(p,2) = raypoints(i,1)
				M(p,1) = raypoints(i,2)
				M(p,1) = M(p,1)*pii/180
				M(p,2) = M(p,2)*pii/180
				tetaS = M(p,2)			! latitude
				phiS = M(p,1)			! longitude
				IF (i.ne.2) THEN
					tetaR = M(p-1,2)	! latitude
					phiR = M(p-1,1)		! longitude
				ELSE
					tetaR = point(2)	! latitude
					phiR = point(1)		! longitude
				END IF 
				
				! Do some trigonometry on the sphere to get angle alpha between two consecutive points.
				! Use the Haversine formula (better conditioned for small distances than the spherical law of cosines)
                       		alpha=2*asin(sqrt((sin((tetaS-tetaR)/2))**2+(cos(tetaR))*(cos(tetaS))*(sin((phiS-phiR)/2))**2))
					
				! Store the distance between 2 points in M(:,3)
				M(p,3)=alpha*er
				angle = angle + alpha
				p1 = (180/pii)*M(p,1)
				p2 = (180/pii)*M(p,2)
				
				! Find in which Voronoi cell the current point is (i.e. to which node the current point is closest)
				CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro(1:2,:)),nnn,nnlist,walk)
					
				! Compute the distance of ray 'rayid' in cell 'node' and store it in RayCel(rayid,node)
				RayCel(rayid,node)=RayCel(rayid,node)+M(p,3) 
				
			END DO
				
			dis(rayid)=angle*er
			IF (sigdep.EQ.2) THEN
				srdist(rayid)=angle*180/pii !*er
			END IF
	
			DEALLOCATE(raypoints)
					
		END DO
	
		IF(wurf.EQ.-1.OR.wurf.EQ.rank+1.OR.rank+1.EQ.slmur.OR.slmur.EQ.-1) THEN
			CLOSE(40)
		END IF
		
		IF(test.EQ.1) THEN
			WRITE(ountc,*)ount
			ountc = ADJUSTL(ountc)		
			string = 'cp rayupdate.dat raysub' // trim(ountc) // '.out'
			CALL SYSTEM( string )
		END IF
		
		np=sum(npoints(:,1),1)-nru
		
		IF(np.NE.p) THEN
			IF(rank.EQ.0) WRITE(*,*)'Something went wrong when updating the rays!'
			IF(rank.EQ.0) WRITE(*,*)'Contrasting values of number of points:',np,'and',p
			!CALL MPI_FINALIZE(ierror) !! comment for former MPI parallelization
			STOP
		END IF
		
		IF(rank.EQ.0) WRITE(*,*)'Total number of ray points:',np
		IF(rank.EQ.0) WRITE(*,*)'________________________________________'
			
		DEALLOCATE(raycoords)
		DEALLOCATE(npoints)
	
		VELU=0
   		
		M(:,1)=(180/pii)*M(:,1)
		M(:,2)=(180/pii)*M(:,2)
		
		!**********************************
		! Get likelihood for the new model
		!**********************************
		
		like=0
		misfit=0
		DO nr=1,nrays
			ttime(nr)=0
			DO i=1,ncell
				ttime(nr) = ttime(nr) + RayCel(nr,i)/Voro(3,i)
			END DO
			IF(sigdep.EQ.0) THEN
				IF (uglpd.EQ.0) THEN
					like = like + (ttime(nr)-dat(nr,5))**2/(2*(sigma(1,dsstat(nr))**2))
				ELSE IF (uglpd.EQ.1) THEN
					like = like + ABS(dat(nr,5)-ttime(nr))/sigma(1,dsstat(nr))
				END IF
				misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(1,dsstat(nr)))**2
			ELSE IF(sigdep.NE.0) THEN
				IF (uglpd.EQ.0) THEN
					like = like + (ttime(nr)-dat(nr,5))**2/(2*(sigma(nr,1)**2))
				ELSE IF (uglpd.EQ.1) THEN
					like = like + ABS(dat(nr,5)-ttime(nr))/sigma(nr,1)
				END IF
				misfit = misfit + ((ttime(nr)-dat(nr,5))/sigma(nr,1))**2
			END IF	
		END DO
		WRITE(*,*)'Updated average misfit for processor',rank+1,':',misfit/nrays
	END IF
END IF
!!! *****************************************************************************************  

!IF (ount.EQ.1000) THEN
!	CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
!	EXIT
!END IF
END DO !end the Reversible Jump sampling	<<<----------------------------------- LOOP OVER SAMPLES ENDS HERE !!! -----------------------------------
!STOP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization


!--------------------------------------------------------------------------

! Deallocate arrays to clear memory
DEALLOCATE(vertices,neighbour,nnn,nnlist,ntwork,worki1,worki2,worki3)
DEALLOCATE(rec,sou,raystat,dat)
IF(ALLOCATED(lengthray)) DEALLOCATE(lengthray)
IF(ALLOCATED(RV)) DEALLOCATE(RV)
IF(ALLOCATED(ldummy)) DEALLOCATE(ldummy)
IF(ALLOCATED(cellogic)) DEALLOCATE(cellogic)
IF(ALLOCATED(RayCel)) DEALLOCATE(RayCel)
IF(ALLOCATED(RayCel_prop)) DEALLOCATE(RayCel_prop)
IF(ALLOCATED(RayCel_pprop)) DEALLOCATE(RayCel_pprop)
IF(ALLOCATED(ttime)) DEALLOCATE(ttime)
IF(ALLOCATED(ttime_prop)) DEALLOCATE(ttime_prop)
IF(ALLOCATED(ttime_pprop)) DEALLOCATE(ttime_pprop)
IF(ALLOCATED(lengthrayv)) DEALLOCATE(lengthrayv)
IF(ALLOCATED(dis)) DEALLOCATE(dis)
IF(ALLOCATED(srdist)) DEALLOCATE(srdist)
IF(ALLOCATED(srdist_prop)) DEALLOCATE(srdist_prop)
IF(ALLOCATED(Voro)) DEALLOCATE(Voro)
IF(ALLOCATED(Voro_prop)) DEALLOCATE(Voro_prop)
IF(ALLOCATED(Voro_pprop)) DEALLOCATE(Voro_pprop)
IF(ALLOCATED(VELU)) DEALLOCATE(VELU)
IF(ALLOCATED(M)) DEALLOCATE(M)
IF(ALLOCATED(nnM)) DEALLOCATE(nnM)
IF(ALLOCATED(dsstat)) DEALLOCATE(dsstat)
  
!--------------------------------------------------------------------------
 

!***************************************************************************
! Save information from all the chains
!***************************************************************************

! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
WRITE(*,*)'Finished chain on processor',rank+1,': collecting information and saving results...'

!----------------------------------------------------------------------------------------------------
! Save last accepted model for next run of the program
!----------------------------------------------------------------------------------------------------
! Save Voronoi cells
OPEN(UNIT=84,FILE=last_tmp,FORM='unformatted',STATUS='replace')
WRITE(84)nlastacc
DO c=1,nlastacc
	WRITE(84)lastacc(2,c),lastacc(1,c),lastacc(3,c)
END DO
CLOSE(84)
DEALLOCATE(lastacc)
! Save data noise parameters
OPEN(UNIT=84,FILE=sigma_tmp,FORM='unformatted',STATUS='replace')
IF (sigdep.EQ.0) THEN
	DO nss=1,nds
		WRITE(84)sigma(1,nss)
	END DO
ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
	DO nss=1,nds
		WRITE(84)aa(nss),bb(nss)
	END DO
ELSE IF (sigdep.EQ.3) THEN
	DO nss=1,nds
		WRITE(84)lambda(nss)
	END DO
END IF
CLOSE(84)
! Save acceptance ratios
OPEN(UNIT=84,FILE=arat_tmp,FORM='unformatted',STATUS='replace')
WRITE(84)totacc
WRITE(84)ounti+ount
WRITE(84)AR(1),AR(2),AR(3)
WRITE(84)ARB(1),ARB(2),ARB(3)
WRITE(84)ARD(1),ARD(2),ARD(3)
WRITE(84)ARM(1),ARM(2),ARM(3)
WRITE(84)AR1p(1),AR1p(2),AR1p(3)
WRITE(84)AR2p(1),AR2p(2),AR2p(3)
WRITE(84)ARV(1),ARV(2),ARV(3)
WRITE(84)AR1v(1),AR1v(2),AR1v(3)
WRITE(84)AR2v(1),AR2v(2),AR2v(3)
WRITE(84)ARS(1),ARS(2),ARS(3)
IF (sigdep.EQ.1.OR.sigdep.EQ.2) WRITE(84)ARSa(1),ARSa(2),ARSa(3)
IF (sigdep.EQ.1.OR.sigdep.EQ.2) WRITE(84)ARSb(1),ARSb(2),ARSb(3)
WRITE(84)randcount
CLOSE(84)
! Deallocate arrays
IF(ALLOCATED(aa)) DEALLOCATE(aa)
IF(ALLOCATED(aa_prop)) DEALLOCATE(aa_prop)
IF(ALLOCATED(bb)) DEALLOCATE(bb)
IF(ALLOCATED(bb_prop)) DEALLOCATE(bb_prop)
IF(ALLOCATED(lambda)) DEALLOCATE(lambda)
IF(ALLOCATED(lambda_prop)) DEALLOCATE(lambda_prop)
!----------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------
! Save binary or text files containing samples
!----------------------------------------------------------------------------------------------------
IF (runnum.EQ.1) THEN
	OPEN(UNIT=50,FILE=trim(sampfile) // '.out' // trim(rankc),STATUS='replace',FORM='unformatted')
	OPEN(UNIT=51,FILE=trim(ncelfile) // '.out' // trim(rankc),STATUS='replace',FORM='formatted')
	IF (sigdep.NE.4) OPEN(UNIT=52,FILE=trim(sigfile) // '.out' // trim(rankc),STATUS='replace',FORM='unformatted')
	IF (sigdep.EQ.4.AND.rank.EQ.0) OPEN(UNIT=52,FILE=trim(sigfile) // '.out',STATUS='replace',FORM='unformatted')
	OPEN(UNIT=53,FILE=trim(misffile) // '.out' // trim(rankc),STATUS='replace',FORM='formatted')
	OPEN(UNIT=54,FILE=trim(aratfile) // '.out' // trim(rankc),STATUS='replace',FORM='formatted')
	DO pp=0,nsample
		WRITE(51,'(2I10)')pp,samplesNcel(pp+1)
		IF (pp.EQ.0) THEN
			DO i=1,samplesNcel(pp+1)
				WRITE(50)firstVoro(1,i),firstVoro(2,i),firstVoro(3,i)
			END DO
		ELSE
			WRITE(50)samplesStep(pp+1,1),samplesStep(pp+1,2),samplesInd(pp+1),samplesVal(pp+1,1),samplesVal(pp+1,2),samplesVal(pp+1,3)
		END IF
		IF (sigdep.NE.4) THEN
			DO nss=1,nds
				IF (sigdep.EQ.0.OR.sigdep.EQ.3) WRITE(52)samplesSigmas(pp+1,nss)
				IF (sigdep.EQ.1.OR.sigdep.EQ.2) WRITE(52)samplesSigmas(pp+1,nss),samplesSigmas(pp+1,nss+nds)
			END DO
		ELSE IF (sigdep.EQ.4.AND.rank.EQ.0.AND.pp.EQ.0) THEN
			WRITE(52)nrays
			DO nr=1,nrays
				WRITE(52)sigma(nr,1)
			END DO
		END IF
		WRITE(53,'(1I10,1F20.6)')pp,samplesMisfit(pp+1)
		WRITE(54,'(1I10,6F14.6)')pp,samplesAratios(pp+1,1),samplesAratios(pp+1,2),samplesAratios(pp+1,3),samplesAratios(pp+1,4),samplesAratios(pp+1,5),samplesAratios(pp+1,6)
	END DO
ELSE IF (runnum.GT.1) THEN
	OPEN(UNIT=50,FILE=trim(sampfile) // '.out' // trim(rankc),STATUS='old',FORM='unformatted',ACCESS='append')
	OPEN(UNIT=51,FILE=trim(ncelfile) // '.out' // trim(rankc),STATUS='old',FORM='formatted',ACCESS='append')
	IF (sigdep.NE.4) OPEN(UNIT=52,FILE=trim(sigfile) // '.out' // trim(rankc),STATUS='old',FORM='unformatted',ACCESS='append')
	OPEN(UNIT=53,FILE=trim(misffile) // '.out' // trim(rankc),STATUS='old',FORM='formatted',ACCESS='append')
	OPEN(UNIT=54,FILE=trim(aratfile) // '.out' // trim(rankc),STATUS='old',FORM='formatted',ACCESS='append')
	DO pp=1,nsample
		WRITE(51,'(2I10)')ounti+pp,samplesNcel(pp)
		WRITE(50)samplesStep(pp,1),samplesStep(pp,2),samplesInd(pp),samplesVal(pp,1),samplesVal(pp,2),samplesVal(pp,3)
		IF (sigdep.NE.4) THEN
			DO nss=1,nds
				IF (sigdep.EQ.0.OR.sigdep.EQ.3) WRITE(52)samplesSigmas(pp,nss)
				IF (sigdep.EQ.1.OR.sigdep.EQ.2) WRITE(52)samplesSigmas(pp,nss),samplesSigmas(pp,nss+nds)
			END DO
		END IF
		WRITE(53,'(1I10,1F20.6)')ounti+pp,samplesMisfit(pp)
		WRITE(54,'(1I10,6F14.6)')ounti+pp,samplesAratios(pp,1),samplesAratios(pp,2),samplesAratios(pp,3),samplesAratios(pp,4),samplesAratios(pp,5),samplesAratios(pp,6)
	END DO
END IF
CLOSE(50)
CLOSE(51)
INQUIRE(52,OPENED=file_open)
IF (file_open) CLOSE(52)
CLOSE(53)
CLOSE(54)

IF (rank.EQ.0) THEN
	OPEN(UNIT=84,FILE='parameters.dat',FORM='unformatted',STATUS='replace')
	WRITE(84)nbproc
	WRITE(84)ounti+ount+1
	WRITE(84)ncell_min,ncell_max
	WRITE(84)mean,theta
	WRITE(84)latmin,latmax
	WRITE(84)longmin,longmax
	WRITE(84)nds
	WRITE(84)sigdep
	WRITE(84)sigma_min,sigma_max
	WRITE(84)aa_min,aa_max
	WRITE(84)bb_min,bb_max
	WRITE(84)lambda_min,lambda_max
	WRITE(84)nt_max,nv_max,nmax
	WRITE(84)nvtr,nvpr
	WRITE(84)latr,longr
	CLOSE(84)
END IF

!----------------------------------------------------------------------------------------------------

! Deallocate samples matrices
IF(ALLOCATED(samplesStep)) DEALLOCATE(samplesStep)
IF(ALLOCATED(samplesInd)) DEALLOCATE(samplesInd)
IF(ALLOCATED(samplesVal)) DEALLOCATE(samplesVal)
IF(ALLOCATED(samplesNcel)) DEALLOCATE(samplesNcel)
IF(ALLOCATED(samplesSigmas)) DEALLOCATE(samplesSigmas)
IF(ALLOCATED(samplesMisfit)) DEALLOCATE(samplesMisfit)
IF(ALLOCATED(samplesAratios)) DEALLOCATE(samplesAratios)
IF(ALLOCATED(sigma)) DEALLOCATE(sigma)
IF(ALLOCATED(sigma_prop)) DEALLOCATE(sigma_prop)


! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
WRITE(*,*)'Accepted samples:',totacc,'/',sampletotal,'on processor',rank+1
! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror) !! uncomment for former MPI parallelization
END IF ! This is the end of the master process run for the ray-routines parallelization

    IF (ount.EQ.nsample) THEN
        DO dest=1,numworkers  	
            CALL MPI_SEND(2, 1, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)		 	! Send run_worker flag for stopping all workers
        END DO
    END IF

!WRITE(*,*) 'MPI closing'


!! ---------------------------------------------- end of master section ------------------------------------------------------!!

!! -------------------------------------------------  worker section ---------------------------------------------------------!!

IF (taskid.GT.0) THEN 
   
    DO WHILE (run_worker.EQ.0)	! <<<--------------------LOOP OVER MPI PROCESSORS STARTS HERE !!! ----------------------------

        CALL MPI_RECV(run_w, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierror) 	! Receive run_worker flag for running worker
	    IF (run_w.EQ.2) THEN
                GOTO 500
   	    END IF  

	    !WRITE(*,*) 'We are processing in worker', taskid   
	    ! TAG  FROM_MASTER = 1
	    ! origin  MASTER = 0
	    ! Receive from master ===> (0,1)
	    ! CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

	    CALL MPI_RECV(params1, 8, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, status, ierror) 	! Brodcast parameters 1 to all processes 
	    CALL MPI_RECV(params2, 11, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierror) 		! Brodcast parameters 2 to all processes 

	    nvtr   = params1(1) 
	    nvpr   = params1(2) 
	    latmax = params1(3) 
	    longmin= params1(4) 
	    latr   = params1(5)  
	    longr  = params1(6)  
	    er     = params1(7) 
	    nbs    = params1(8)  

	    nsou   = params2(1) 
	    nrec   = params2(2) 
	    wrgf   = params2(3) 
	    gdt    = params2(4) 
	    gdp    = params2(5) 
	    sgref  = params2(6) 
	    dl     = params2(7) 
	    erg    = params2(8) 
	    fms    = params2(9) 
	    rank   = params2(10)-1 ! The rank of the chain (outer parallelization)
	    wcount = params2(11)

	    IF(ALLOCATED(offset_worker)) DEALLOCATE(offset_worker)   
	    IF(ALLOCATED(nsou_worker)) DEALLOCATE(nsou_worker)
	    IF(ALLOCATED(rec_worker)) DEALLOCATE(rec_worker)      
	    IF(ALLOCATED(sou_worker)) DEALLOCATE(sou_worker)      
	    IF(ALLOCATED(raystat_worker)) DEALLOCATE(raystat_worker)
	    IF(ALLOCATED(VELI_worker)) DEALLOCATE(VELI_worker)

	    ALLOCATE( offset_worker(numworkers,2) )
	    ALLOCATE( nsou_worker(numworkers) )                      
	    ALLOCATE( rec_worker(nrec, 2) )           
	    ALLOCATE( sou_worker(nrec, 2) )           
	    ALLOCATE( raystat_worker(nsou*nrec, 2) )    
	    ALLOCATE( VELI_worker(nvpr+2,nvtr+2) )    

	    counter_1 = 2*nsou*nrec
	    counter_2 = (nvpr+2)*(nvtr+2)
	    CALL MPI_RECV(offset_worker, numworkers*2, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierror) ! Brodcast parameters 2 to all processes 
	    CALL MPI_RECV(nsou_worker, numworkers, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierror)
	    CALL MPI_RECV(rec_worker, nrec*2, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, status, ierror)
	    CALL MPI_RECV(sou_worker, nsou*2, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, status, ierror) 
	    CALL MPI_RECV(raystat_worker, counter_1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierror)
	    CALL MPI_RECV(VELI_worker, counter_2, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, status, ierror)    

	    ! The processor number is taskid
	    nsou_subset = nsou_worker(taskid)
	    offset = offset_worker(taskid, 1)
	    offset_end  = offset_worker(taskid, 2)

	    IF(ALLOCATED(sou_subset)) DEALLOCATE(sou_subset) 
	    IF(ALLOCATED(raystat_subset)) DEALLOCATE(raystat_subset)   

	    ALLOCATE( sou_subset(nsou_subset, 2) )           
	    ALLOCATE( raystat_subset(nsou_subset*nrec,2) )  				! raystat has a switch for each raypath (valid = 1, non-valid =0)
	    
	    sou_subset = sou_worker( offset+1:offset_end, 1:2 )     			! sources corresponding to each worker, selected starting with the next source after the last set [offset+1] ).                    
	    raystat_subset = raystat_worker( offset*nrec+1 : offset_end*nrec,  1:2 )   	! raystat: switch in column 1, valid ray number in column 2 for the set of sources of the worker.  

	    !WRITE(*,*) 'This is samples number', wcount
	    !WRITE(*,*) '1. Number of sources for processor',taskid, 'is', nsou_subset    
	    !WRITE(*,*) '2. The offset_worker is', offset_worker(1,1), offset_worker(2,1), offset_worker(3,1), offset_worker(4,1), offset_worker(5,1), offset_worker(6,1), offset_worker(7,1)
	    !WRITE(*,*) '3. The nsou_worker is', nsou_worker(1), nsou_worker(2), nsou_worker(3), nsou_worker(4), nsou_worker(5), nsou_worker(6), nsou_worker(7)    
	    !WRITE(*,*) '4. VELI_worker last element is', VELI_worker(nvpr+2,nvtr+2), 'from the master'
	    !WRITE(*,*) 'Processor', taskid
	    !WRITE(*,*) '5. Received the receiver 1', rec_worker(1,1), 'from the master'
	    !WRITE(*,*) '6. Received the source 1', sou_worker(1,1), 'from the master'
	    !WRITE(*,*) '7. raystat 1 is', raystat_worker(1,1), raystat_worker(1,2), 'from the master'
	    !WRITE(*,*) '8. Received offset, offset_end value', offset, offset_end, 'from the master'
	    !WRITE(*,*) '9. The source_subset 1 is ', sou_subset(1,1), sou_subset(1,2), 'for the processor', taskid
	    !WRITE(*,*) '10.The raystat subset 1 end is', raystat_subset(1,1), raystat_subset(1,2) 
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    ! The routine below outputs 
	    ! nru_subset    npoints_subset    raycoords_subset    crazyray_subset
	    ! (Parallelization objective: obtaine the final arrays 
	    ! nru    npoints    raycoords    crazyray

	    nru_subset = sum(raystat_subset(:,1), 1)       			! number of active rays (switch=1) for the subset of the worker
	    IF(ALLOCATED(npoints_subset)) DEALLOCATE(npoints_subset) 
	    IF(ALLOCATED(raycoords_subset)) DEALLOCATE(raycoords_subset)   

	    CALL modrays_subset(nsou_subset,sou_subset(:,1),sou_subset(:,2),&
	    	 nrec,rec_worker(:,1),rec_worker(:,2),&
	    	 raystat_subset,-1,&
	     	 nvtr,nvpr,latmax,longmin,latr,longr,VELI_worker,&
	    	 gdt,gdp,&
	    	 sgref,&
	   	 dl,erg,&
		 er,&
		 fms,&
		 nbs,&
		 rank+1)    

	    numraycoords_subset = sum( npoints_subset( 1:nru_subset , 1 ) ) ! size of raycoords_subset found by adding the first column of the corresponding subset of npoints = number of total ray coordinates in the subset.    

	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    !CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
	    !WRITE(*,*) 'After modrays in worker'
	    !WRITE(*,*) 'Processor', taskid
	    !WRITE(*,*) 'Received offset, offset_end value', offset, offset_end, 'from the master'

	    ! TAG FROM_WORKER = 2
	    ! MASTER=0 
	    ! Send from worker to master ===> (0,2)

	    CALL MPI_SEND(offset_worker, 2*numworkers, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, ierror)  
	    CALL MPI_SEND(nru_subset, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, ierror)			
	    CALL MPI_SEND(npoints_subset, 2*nru_subset, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, ierror)      
	    CALL MPI_SEND(numraycoords_subset, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, ierror)  
	    CALL MPI_SEND(raycoords_subset, numraycoords_subset*2, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, ierror)
	    CALL MPI_SEND(crazyray_subset, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, ierror)

	    !WRITE(*,*) 'Sample number', wcount, 'of a total of', nsample, 'samples.'
	    !WRITE(*,*) 'numraycoords', numraycoords_subset, 'processor', taskid
	    !WRITE(*,*) 'size of npoints_subset', SIZE(npoints_subset)
	    !WRITE(*,*) 'size of raycoords_subset', SHAPE(raycoords_subset)
	    !WRITE(*,*) 'nru_subset', nru_subset, 'for processor', taskid
	    !WRITE(*,*) 'npoints_subset 1', npoints_subset(1,1), 'processor', taskid 
	    !WRITE(*,*) 'numraycoords', numraycoords_subset, 'processor', taskid
	    !WRITE(*,*) 'raycoords_subset 1', raycoords_subset(1,1), 'processor', taskid
	    !WRITE(*,*) 'crazyray', crazyray_subset, 'processor', taskid
	    !WRITE(*,*) 'All results sent to master.  End of processor id', taskid
	    !OK All checked!
    END DO	   
500 WRITE(*,*) 'Finished loop at worker', taskid
END IF 

!! -----------------------------------------------  end worker section -------------------------------------------------------!!


 CALL MPI_FINALIZE(ierror)! Terminate the parallelization for the ray-routines !!
 CALL cpu_time(t2) ! TOC - Stop counting time and display

IF (taskid.EQ.0) WRITE(*,*)
IF (taskid.EQ.0) WRITE(*,*)'------------------------------------------------------------------------'
IF (taskid.EQ.0) WRITE(*,*)'Total ',t2-t1, 'seconds'
IF (taskid.EQ.0) WRITE(*,*)'Modrays prll', t_rays_test, 'Orig modrays', t_rays_test_base, 'seconds'

END ! end the main program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!                SUBROUTINES                  !!!!!!!!!!!!!! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------!
!       Random number generator for a Gaussian distribution        !
!       (adapted from Numerical Recipes)                           !
! -----------------------------------------------------------------!
FUNCTION GASDEV(randcount)
      
      USE mt19937
      DOUBLE PRECISION GASDEV
      DOUBLE PRECISION v1,v2,r,fac
      INTEGER randcount
 
      !IF (idum.lt.0) iset=0
10     v1=2*grnd()-1
       randcount=randcount+1
       v2=2*grnd()-1
       randcount=randcount+1
       r=v1**2+v2**2
       IF(r.ge.1.or.r.eq.0) GOTO 10
       fac=sqrt(-2*log(r)/r)
       GASDEV=v2*fac

      RETURN
END
!------------------------------------------------------------------!  

!------------------------------------------------------------------!
!                      Update raypath geometry                     !
!------------------------------------------------------------------!
!!! Subroutine added by EG
!!! Update those raypaths that go through cells that are modified in the new proposed model
!!! Requires the use of library libuprays.a  
	
	!! This subroutine is used in case of birth, death, move, velocity change of a cell 
	SUBROUTINE updaterays(RayCel_prop,RayCel,Voro_prop,ncell,ncell_prop,cellogic,ind,pn,sn, t_upr, t_upr_test)
		
		USE mksamples_prllvars 
		USE MPI
		USE fm2dssvars
		USE fm2dout
		
		IMPLICIT NONE
		
		INTEGER			::	i,ii,iii,j,jj,jjj,node,pp,pu,nr,nc,nrr,ind,rayid,raylength,pti,ptf,pn,sn
		INTEGER			::	nsin,nrin,ncell,ncell_prop
		REAL(KIND=ii10)		::	p1,p2,point(2), t1_upr, t2_upr, t1_upr_test, t2_upr_test
		LOGICAL 		::	cellogic(ncell_max)
		REAL(KIND=ii10)	  	::	tetaS,tetaR,phiS,phiR,alpha,angle,RayCel(nrays,ncell)

		INTEGER :: status(MPI_STATUS_SIZE)
		
		REAL(KIND=ii5), DIMENSION(:,:), ALLOCATABLE 	::	raypoints
		REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 	::	Mu
		REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 	::	RAYUP
		REAL(KIND=ii10), DIMENSION(1:,1:) 		:: 	Voro_prop
		REAL(KIND=ii10), DIMENSION(1:,1:)		::	RayCel_prop
		!INTEGER, DIMENSION(:,:), ALLOCATABLE		::	switch    	This variable has been defined in the module: mksamples_prllvars
		CHARACTER (LEN=30) 				::	pnc,snc

		REAL(KIND=ii10), INTENT(out) :: t_upr, t_upr_test
						
		!WRITE(*,*)'Hello from subroutine UPDATERAYS'
				
		IF(ALLOCATED(switch)) DEALLOCATE(switch)
		ALLOCATE(switch(nsou*nrec,2))

		switch=raystat
		switch(:,1)=0
		pp=0
		node=1

                t_upr = 0
                t1_upr = 0 
                t2_upr = 0 

                t_upr_test = 0
                t1_upr_test = 0 
                t2_upr_test = 0 
		
		crazyray=0
				
		IF(ALLOCATED(RAYUP)) DEALLOCATE(RAYUP)
		IF(ALLOCATED(npoints)) DEALLOCATE(npoints)
		IF(ALLOCATED(raycoords)) DEALLOCATE(raycoords)
		
		pp=0
		! If only the raypaths going through the affected cell and its neighbours should be updated...
		IF (uar.EQ.0) THEN	
			DO nr=1,nrays
				DO nc=1,ncell
					IF ((cellogic(nc)).AND.(RayCel(nr,nc).NE.0)) THEN
					
						pp=pp+1
						
						DO nrr=1,nsou*nrec
					
							IF(switch(nrr,2).EQ.nr) THEN
								switch(nrr,1)=1
								EXIT
							END IF
					
						END DO
			
						RayCel_prop(nr,:) = 0
						
						EXIT
						
					END IF
				END DO
			END DO
		! ...otherwise, if all raypaths should be updated...
		ELSE IF (uar.EQ.1) THEN
			DO nr=1,nrays
				pp=pp+1
				DO nrr=1,nsou*nrec
					IF(switch(nrr,2).EQ.nr) THEN
						switch(nrr,1)=1
						EXIT
					END IF
				END DO
			END DO
		END IF
			
		nru=pp
			
		IF(nru.GT.0) THEN
		
			ALLOCATE(RAYUP(nvpr+2,nvtr+2))	! This is the matrix containing the proposed model velocities that will be used to update the raypaths
			
			! Calculate Delaunay triangulation
			CALL delaun(DBLE(Voro_prop(1:2,:)),ncell_prop,neighbour,vertices,nt,nt_max,&
                	      		worki1,worki2,worki3,eps,nv_max,&
                	      		0,ldummy,0,0,0,ptsok)
			CALL build_nv(ncell_prop,vertices,nt,ncell_max,nmax,&
				  neighbour,nnn,nnlist,ntwork)
				
			iii=0
			DO i=0,nvpr+1
				jjj=0
				iii=iii+1
				DO j=0,nvtr+1
					jjj=jjj+1					
					ii=i
					jj=j
					IF (i.EQ.0)   	ii=1
					IF (j.EQ.0)   	jj=1
					IF (i.EQ.nvpr+1) 	ii=nvpr
					IF (j.EQ.nvtr+1) 	jj=nvtr
						
					p1=longmin+(ii-1)*longr
					p2=latmax-(jj-1)*latr
					
					node=1
					
					CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro_prop(1:2,:)),nnn,nnlist,walk)
					
					RAYUP(iii,jjj) = (Voro_prop(3,node))
					
				END DO
			END DO
				
			
			IF(allocated(npoints)) write(*,*)'npoints allocated'
			IF(allocated(raycoords)) write(*,*)'raycoords allocated'



			! This is the third chance for MPI parallelization


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	        	! The section below corresponds to the traditional case for the modrays algorithm. Use for algorithm output comparison.

			IF(TEST_MODRAYS.EQ.1) THEN 

				t1_upr = 0
		  	        t2_upr = 0
				CALL cpu_time(t1_upr) ! TOC - Start counting time
				
				CALL modrays(nsou,sou(:,1),sou(:,2),&
					nrec,rec(:,1),rec(:,2),&
					switch,-1,&   					! switch replaces raystat 
					nvtr,nvpr,latmax,longmin,latr,longr,RAYUP,&	! RAYUP replaces VELI
					gdt,gdp,&
					sgref,&
					dl,erg,&
					er,&
					fms,&
					nbs,&
					pn)						! pn is rank+1

				CALL cpu_time(t2_upr) ! TOC - Stop counting time 
				! Calculate cumulative time spent in modrays routine 
		  		t_upr =  t2_upr-t1_upr

			        IF(ALLOCATED(npoints_base)) DEALLOCATE(npoints_base)
				IF(ALLOCATED(raycoords_base)) DEALLOCATE(raycoords_base)
			        nru_base = nru	
				crazyray_base = crazyray

				!WRITE(*,*) 'Done with test 2!', SIZE(npoints), SIZE(raycoords)

			        ALLOCATE( npoints_base( (SIZE(npoints))/2,2 ) )	
			        ALLOCATE( raycoords_base( (SIZE(raycoords)/2),2 ) )
	       		        npoints_base = npoints
		       		raycoords_base = raycoords

				IF(ALLOCATED(npoints)) DEALLOCATE(npoints)
				IF(ALLOCATED(raycoords)) DEALLOCATE(raycoords)

			END IF 
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		        ! Inner parallelization starts here.  
			! We start by sending the source sets to the workers for the modrays routine for calculation.  We will send the partial  
			! information required by the modrays routines to every worker.      	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Time counter

			t1_upr_test = 0
	  	        t2_upr_test = 0
			CALL cpu_time(t1_upr_test) ! TOC - Start counting time
		
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!WRITE(*,*) 'Start parallelization inside updaterays, sample', sn

			crazyray = 0 
			IF(ALLOCATED(offset_master)) DEALLOCATE(offset_master)
			IF(ALLOCATED(nsou_master)) DEALLOCATE(nsou_master)

			ALLOCATE(offset_master(numworkers,2))
			ALLOCATE(nsou_master(numworkers))

			! (We store data in several matrices and send to the worker tasks).  

			!WRITE(*,*) 'Number of sources inside updaterays', nsou
			!WRITE(*,*) 'Number of workers inside updaterays', numworkers
			avesou = nsou/numworkers  ! Number of sources in each subset, i.e. per worker excluding remainder.   e.g. 11/3 = 3. remainder 2.  
			extrasou = MOD(nsou,numworkers)		

		        ! The offsets are the indices used for transiting along the arrays
			! In the first loop we calculate the offsets that will be used by the workers.   
			! In the first loop we calculate the nsources that will be send to each worker.   
	 
		        offset = 0        

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		        DO dest=1,numworkers  
			    ! Here we need to count the number of rays corresponding to the sources sent to the workers.  Number of sources sent is nsou_sent.             
		            IF(dest.LE.extrasou) THEN     														
				nsou_sub = avesou+1  	! nsou_subset are the number of sources sent to the worker for calculating the modrays algorithm.       		
		            ELSE 
		                nsou_sub = avesou     			
			    END IF
		            offset_end = offset + nsou_sub   ! This defines the interval we are sending.  Ex:  offset = 0, nsou_subset = 3, offset_end = 3.                      
			    nsou_master(dest) = nsou_sub
			    offset_master(dest, 1) = offset
	 		    offset_master(dest, 2) = offset_end
			    ! ATTENTION offsetpts and offsetcoords are calculated and looped in the second master processor loop, but after receiving the npoints and raycoords results from the worker.
		            offset = offset_end	
	  	        END DO

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! TAG FROM_MASTER = 1              
		        ! TO: DEST 
			! Send to worker ===> (dest,1)

		        !Brodcast

			!nvx,nvz,goxd,gozd,dvxd,dvzd,velv
			params1(1) = nvtr
			params1(2) = nvpr
			params1(3) = latmax
			params1(4) = longmin
			params1(5) = latr
			params1(6) = longr
			params1(7) = er ! Earth radius
			params1(8) = nbs

		        wrgf = -1 	                      
			params2(1) = nsou
			params2(2) = nrec
			params2(3) = wrgf          	! wrgf = write geometries to file?
			params2(4) = gdt		! Grid dicing
			params2(5) = gdp		! Grid dicing
			params2(6) = sgref		! Apply source grid refinement? (0=no,1=yes)
			params2(7) = dl			! Dicing level and extent of refined grid
			params2(8) = erg
			params2(9) = fms		! use first-order(0) or mixed-order(1) scheme
			params2(10)= pn			! Processor rank+ (in this case for embarrasingly parallel processing).  We send pn instead of rank+1.
			params2(11)= sn 		! Sample number! 

	   	        !counter_5 = 2*nsou*nrec         ! The counter for switch.  Optional.
		        counter_2 = (nvpr+2)*(nvtr+2)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		        !WRITE(*,*) 'Before broadcasting to the workers'

		        DO dest=1,numworkers  
			    CALL MPI_SEND(0, 1, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)	   	 	        ! Send run_worker flag for running worker
			    CALL MPI_SEND(params1, 8, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror) 		! Brodcast parameters 1 to all processes 
	    	            CALL MPI_SEND(params2, 11, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror) 			! Brodcast parameters 2 to all processes 
			    CALL MPI_SEND(offset_master, 2*numworkers, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)	
		            CALL MPI_SEND(nsou_master, numworkers, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)
			    CALL MPI_SEND(rec, nrec*2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)
			    CALL MPI_SEND(sou, nsou*2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)		   		   
			    CALL MPI_SEND(switch, 2*nsou*nrec, MPI_INTEGER, dest, 1, MPI_COMM_WORLD, ierror)		! Here we send switch instead of raystat
			    CALL MPI_SEND(RAYUP, counter_2, MPI_DOUBLE_PRECISION, dest, 1, MPI_COMM_WORLD, ierror)       ! Here we send RAYUP instead of VELI 
		        END DO

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		        !WRITE(*,*) 'All variables and arrays have been broadcasted to the workers'						

			! TAG = FROM_WORKER
		        ! FROM: Worker
		        ! Receive in master from worker ===> (from_worker,2)		

		        offsetpts 	= 0  
		        offsetcoords 	= 0  

			! Now we loop over all workers receiving the information for building the complete arrays.            

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			DO  i=1,numworkers
			    ! ATTENTION offsetpts and offsetcoords are calculated and looped here, but are updated receiving the nru, npoints and raycoords results from the worker.
		            IF(ALLOCATED(offset_master)) DEALLOCATE(offset_master)
		            IF(ALLOCATED(npoints_master)) DEALLOCATE(npoints_master)
	  		    ALLOCATE( offset_master(numworkers,2) )  		    
			    from_worker = i		
	 
			    !TAG pair:  (i,2)		
		            CALL MPI_RECV(offset_master, 2*numworkers, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)
		            CALL MPI_RECV(nru_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)	                    
			    ALLOCATE( npoints_master( nru_master , 2 ) )			! array of length nru_subset, filled w. number of points that each ray contains and ray ID corresponding to the subset chosen.
			    CALL MPI_RECV(npoints_master, 2*nru_master, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)  
			    numraycoords_master = sum( npoints_master( 1:nru_master , 1 ) )	
	    	            CALL MPI_RECV(numraycoords_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)  
		          
		            IF(ALLOCATED(raycoords_master)) DEALLOCATE(raycoords_master)
		            ALLOCATE( raycoords_master( numraycoords_master, 2  )  )		! allocate array of ray coordinates for the sources subset.                   
			    CALL MPI_RECV(raycoords_master, numraycoords_master*2, MPI_DOUBLE_PRECISION, from_worker, 2, MPI_COMM_WORLD, status, ierror)
			    CALL MPI_RECV(crazyray_master, 1, MPI_INTEGER, from_worker, 2, MPI_COMM_WORLD, status, ierror)
			    offsetpts_end = offsetpts + nru_master
		            offsetcoords_end = offsetcoords + numraycoords_master
  			    ! Next we fill in the raycoords array

			    IF (i.EQ.1) THEN	
			                     		    
		                ALLOCATE( npoints(nru_master, 2 ) )			! npoints is allocated here for the first time
				npoints = npoints_master				! npoints is filled here for the first time
		                ALLOCATE( raycoords( numraycoords_master, 2  )  )      	! allocate array of ray coordinates for the sources subset for the first time
		                raycoords = raycoords_master				! raycoords is filled here for the first time				

		                IF(ALLOCATED(npoints_prev))DEALLOCATE( npoints_prev )  
		                IF(ALLOCATED(raycoords_prev))DEALLOCATE( raycoords_prev )  
				ALLOCATE( npoints_prev( nru_master , 2 )  )
				ALLOCATE( raycoords_prev( numraycoords_master , 2 ) )
				npoints_prev   = npoints				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.
				raycoords_prev = raycoords				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.


			    ELSE IF (i.GT.1) THEN
		                DEALLOCATE( raycoords ) 
		                DEALLOCATE( npoints ) 
		                ALLOCATE( npoints( offsetpts_end, 2  )  )  	 	   ! allocate array of ray coordinates for the sources subset.		
		                ALLOCATE( raycoords( offsetcoords_end, 2  )  )      	   ! allocate array of ray coordinates for the sources subset.		
		                npoints(1:offsetpts, 1:2) = npoints_prev		   ! npoints is partially filled here with the new npoints subset that has been received recently.	                          
				npoints(offsetpts+1: offsetpts_end, 1:2) = npoints_master  ! npoints is partially filled here with the new npoints subset that has been received recently.	                          
	 		        raycoords(1:offsetcoords, 1:2 ) = raycoords_prev
		                raycoords(offsetcoords+1: offsetcoords_end, 1:2) = raycoords_master ! raycoords is partially filled here with the new raycoords subset that has been received recently.	                          

		                IF(ALLOCATED(npoints_prev))DEALLOCATE( npoints_prev )  
		                IF(ALLOCATED(raycoords_prev))DEALLOCATE( raycoords_prev )
		    	        ALLOCATE( npoints_prev( offsetpts_end , 2 )  )
				ALLOCATE( raycoords_prev( offsetcoords_end , 2 ) )
				raycoords_prev = raycoords				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.
				npoints_prev   = npoints				  ! Save the previous raycoords array for saving the result in the fisrt "portion" of raycoords.

			     END IF 

		            offsetpts    = offsetpts_end 
		            offsetcoords = offsetcoords_end  !OK ALL REVISED UP TO THIS POINT.
	     		    numraycoords = sum( npoints( : , 1 ) )                    ! size of raycoords_subset found by adding the first column of the corresponding subset of npoints = number of total ray coordinates in the subset.			
		            crazyray = crazyray + crazyray_master

			    ! Printing
			    !WRITE(*,*)'offset after 26 equals', offset_master(1,1), 'in processor = ',i
			    !WRITE(*,*)'offset after 26 equals', offset_master(3,1), 'in processor = ',i
			    !WRITE(*,*)'offset after 26 equals', offset_master(5,1), 'in processor = ',i
			    !WRITE(*,*)'offset after 26 equals', offset_master(7,1), 'in processor = ',i	
			    !WRITE(*,*)'nru_master is ', nru_master, 'in processor = ',i
			    !WRITE(*,*)'npoints_master is ', npoints_master(1,1), 'in processor = ',i
			    !WRITE(*,*)'numraycoords_master is ', numraycoords_master, 'in processor = ',i
			    !WRITE(*,*)'raycoords_master is ', raycoords_master(1,1), 'in processor = ',i
			    !WRITE(*,*)'crazyray_master is ', crazyray_master, 'in processor = ',i
	   		    !WRITE(*,*)'raycoords is ', crazyray_master, 'in processor = ',i
			    !WRITE(*,*)'crazyray_master is ', crazyray_master, 'in processor = ',i
			    !WRITE(*,*)'Partial numraycoords at master is', numraycoords, 'processor', taskid
			    !WRITE(*,*)'Partial size of npoints at master is ', SHAPE(npoints)
			    !WRITE(*,*)'Partial size of raycoords', SHAPE(raycoords)
			END DO

   		        CALL cpu_time(t2_upr_test) ! TOC - Stop counting time 
			! Calculate cumulative time spent in modrays routine when parallelized
  			t_upr_test =  t2_upr_test-t1_upr_test
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		        ! Inner parallelization ends here.  
			! We obtain as a result:
			! npoints, raycoords, crazyray.     	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     		        !Here finally we compare the results of both algorithms 

			IF(TEST_MODRAYS.EQ.1) THEN 

			    !WRITE(*,*)'crazyray_master in uprays is ', crazyray_master
			    !WRITE(*,*)'nru_master in uprays is ', offsetpts_end
			    !WRITE(*,*)'Size of (parallel) npoints in uprays is ', SHAPE(npoints)
			    !WRITE(*,*)'Size of (parallel) raycoords in uprays is', SHAPE(raycoords)
			    !WRITE(*,*)'original crazyray in uprays is ', crazyray_base
			    !WRITE(*,*)'original nru in uprays is ', nru_base
			    !WRITE(*,*)'Size of (original) npoints in uprays is ', SHAPE(npoints_base)
			    !WRITE(*,*)'Size of (original) raycoords in uprays is', SHAPE(raycoords_base)
			    !WRITE(*,*)

			    npoints_base = npoints_base - npoints
	 		    raycoords_base = raycoords_base - raycoords 

			    npoints_diff1 = abs( sum( npoints_base(:, 1) ) ) 
		            npoints_diff2 = abs( sum( npoints_base(:, 2) ) )
			    raycoords_diffx = abs( sum( raycoords_base(:, 1) ) )
			    raycoords_diffy = abs( sum( raycoords_base(:, 2) ) )

			    !npoints_base_sum1 = abs( sum( npoints_base(:, 1) ) ) 
		            !npoints_base_sum2 = abs( sum( npoints_base(:, 2) ) )
			    !raycoords_base_sumx = abs( sum( raycoords_base(:, 1) ) )
			    !raycoords_base_sumy = abs( sum( raycoords_base(:, 2) ) )

			    !npoints_sum1 = abs( sum( npoints(:, 1) ) ) 
		            !npoints_sum2 = abs( sum( npoints(:, 2) ) )
			    !raycoords_sumx = abs( sum( raycoords(:, 1) ) )
			    !raycoords_sumy = abs( sum( raycoords(:, 2) ) )

			    WRITE(*,*) 'The sum of the differences of npoints is (', npoints_diff1,',', npoints_diff2, ') in uprays' 
			    WRITE(*,*) 'The sum of the differences of raycoords is (', raycoords_diffx,',', raycoords_diffy, ') in uprays' 
   			    !WRITE(*,*) 'Results parallelization inside updaterays, sample', sn
			    !WRITE(*,*) 'The base sum in npoints is (', npoints_base_sum1,',', npoints_base_sum2, ') in uprays' 
			    !WRITE(*,*) 'The base sum in raycoords is (', raycoords_base_sumx,',', raycoords_base_sumy, ') in uprays' 
			    !WRITE(*,*) 'The sum in npoints is (', npoints_sum1,',', npoints_sum2, ') in uprays' 
			    !WRITE(*,*) 'The sum in raycoords is (', raycoords_sumx,',', raycoords_sumy, ') in uprays' 

			END IF

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			IF(crazyray.GT.0) THEN
				IF(sse.EQ.1) THEN
				
					! Save the updated geometry to file
					WRITE(pnc,*)pn
					pnc = ADJUSTL(pnc)
					WRITE(snc,*)sn
					snc = ADJUSTL(snc)
					errout = ADJUSTL(errout)
					errvtx = trim(errout) // '_s' // trim(snc) // '_p' // trim(pnc) // '.vtx'
					errvor = trim(errout) // '_s' // trim(snc) // '_p' // trim(pnc) // '.vor'
			
					OPEN(UNIT=72,FILE=errvtx,STATUS='replace')
					WRITE(72,*)nvtr,nvpr
					WRITE(72,'(3f14.8)')latmax,longmin
					WRITE(72,'(3f14.8)')latr,longr
					WRITE(72,'(1X)')
			
					OPEN(UNIT=74,FILE=errvor,STATUS='replace')
					WRITE(74,*)nvtr,nvpr
					WRITE(74,'(3f14.8)')latmax,longmin
					WRITE(74,'(3f14.8)')latr,longr
					WRITE(74,'(1X)')
					WRITE(74,*)ncell				
		
		
					iii=0
					DO i=0,nvpr+1
						jjj=0
						iii=iii+1
						DO j=0,nvtr+1
							jjj=jjj+1					
						
							WRITE(72,*)RAYUP(iii,jjj)

						END DO

						WRITE(72,'(1X)')
					
					END DO
		
					CLOSE(72)
				
				
					DO nc=1,ncell
	
						WRITE(74,*)Voro_prop(2,nc),Voro_prop(1,nc),Voro_prop(3,nc)
					
					END DO
		
					CLOSE(74) ! close the file
				END IF
				WRITE(*,*)
				WRITE(*,*)'Current sample causes raypath errors.'
				RETURN
			END IF
			
			pti=0
			ptf=0
			
			IF(wurf.EQ.-1.OR.pn.EQ.wurf) THEN
				WRITE(pnc,*)pn
				pnc = ADJUSTL(pnc)
				OPEN(UNIT=40,FILE='rayupdate_p' // trim(pnc) // '.dat',FORM='unformatted',STATUS='unknown')
				WRITE(40)nru
			END IF
			
			DO nr=1,nru
			
				rayid=npoints(nr,2)
				raylength=npoints(nr,1)
				!WRITE(*,*)'Ray',rayid,'has',raylength,'points'
				
				pti=pti+1
				ptf=pti+raylength-1
				
				ALLOCATE(raypoints(raylength,2))
				raypoints=raycoords(pti:ptf,:)
				
				pti=ptf
					
				ALLOCATE(Mu(raylength-1,3))
				
				lengthrayv(rayid)=raylength
				
				angle=0
					
				point(2) = raypoints(1,1)
				point(1) = raypoints(1,2)
				point(1) = point(1)*pii/180
				point(2) = point(2)*pii/180
					
				pu = 0
				Mu = 0
					
				RayCel_prop(rayid,:) = 0 
				
				IF(wurf.EQ.-1.OR.pn.EQ.wurf) THEN
					WRITE(40)raylength
					WRITE(40)raypoints(1,1),raypoints(1,2)
				END IF
					
				DO i=2,raylength
				
					IF(wurf.EQ.-1.OR.pn.EQ.wurf) THEN
						WRITE(40)raypoints(i,1),raypoints(i,2)
					END IF
				
					pu=pu+1
					Mu(pu,2) = raypoints(i,1)
					Mu(pu,1) = raypoints(i,2)
					Mu(pu,1) = Mu(pu,1)*pii/180
					Mu(pu,2) = Mu(pu,2)*pii/180
					tetaS = Mu(pu,2)		! latitude
					phiS = Mu(pu,1)			! longitude
					IF (i.ne.2) THEN
						tetaR = Mu(pu-1,2)	! latitude
						phiR = Mu(pu-1,1)	! longitude	
					ELSE
						tetaR = point(2)	! latitude
						phiR = point(1)		! longitude
					END IF 
						
					! Do some trigonometry on the sphere to get angle alpha between two consecutive points.
					! Use the Haversine formula (better conditioned for small distances than the spherical law of cosines)
                        		alpha=2*asin(sqrt((sin((tetaS-tetaR)/2))**2+(cos(tetaR))*(cos(tetaS))*(sin((phiS-phiR)/2))**2))
				
					! Store the distance between 2 points in Mu(:,3)
					Mu(pu,3)=alpha*er
					angle = angle + alpha
					p1 = (180/pii)*Mu(pu,1)
					p2 = (180/pii)*Mu(pu,2)
				
					! Find in which Voronoi cell the current point is (i.e. to which node the current point is closest)
					CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro_prop(1:2,:)),nnn,nnlist,walk)
					
					! Compute the distance of ray 'rayid' in cell 'node' and store it in Raycel(rayid,node)
					RayCel_prop(rayid,node)=RayCel_prop(rayid,node)+Mu(pu,3) 
					
				END DO
				
				dis(rayid)=angle*er
				IF (sigdep.EQ.2) THEN
					srdist_prop(rayid)=angle*180/pii !*er
				END IF
	
				DEALLOCATE(Mu)
					
				DEALLOCATE(raypoints)
					
			END DO
			
			IF(wurf.EQ.-1.OR.pn.EQ.wurf) THEN
				CLOSE(40)
			END IF
			
			DEALLOCATE(raycoords)
			
			DEALLOCATE(RAYUP)
			
		ELSE
		
			RayCel_prop=RayCel
			
		END IF
		
	END SUBROUTINE updaterays	
!------------------------------------------------------------------!
