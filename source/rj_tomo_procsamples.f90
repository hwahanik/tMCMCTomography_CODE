!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This progranm can be used to process the outputs of the rj_tomo_mksamples.f90 code	!!
!! performing transdimensional 2D tomography with the reversible-jump algorithm.	!!
!!											!!
!! Inputs, outputs and parameters are defined in the input file procsamples.in		!!
!! 											!!
!! Erica Galetti (EG), April 2013							!!
!! Contact: erica.galetti@ed.ac.uk							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM procsamples

IMPLICIT NONE

! Define two types of real numbers (for consistency with FMST code)
INTEGER, PARAMETER :: ii5=SELECTED_REAL_KIND(5,10)
INTEGER, PARAMETER :: ii10=SELECTED_REAL_KIND(10,100)
! File names
CHARACTER (LEN=30) ::	cdum
CHARACTER (LEN=30) ::	ncellsout,evidout,noiseout,aveout,stdevout,medout,maxout,misfout
CHARACTER (LEN=30) ::	ountc,rankc
CHARACTER (LEN=30) ::	firstout,lastout,voroout,nodeout,avwtmp,postout
CHARACTER (LEN=40) ::	firstvtx,firstvor,lastvtx,lastvor
CHARACTER (LEN=40) ::	voroini,vorofin,nodeini,nodefin,nodeens
CHARACTER (LEN=30) ::	sampfile,ncelfile,sigfile,misffile,aratfile
! Markov chain parameters
INTEGER 	::	burn_in   	! Burn-in period
INTEGER 	::	sampletotal,ss	! Total number of samples (including burn-in)
INTEGER 	::	thinn,thin	! Thinning of the chain 
INTEGER 	::	nbproc,rank,rr,nchains
! Prior parameters
REAL(KIND=ii10)	::	mean,theta
REAL(KIND=ii10)	::	latmin,latmax
REAL(KIND=ii10)	::	longmin,longmax
INTEGER		::	nds,nss,sigdep
REAL(KIND=ii10)	::	sigma_min,sigma_max,aa_min,aa_max,bb_min,bb_max,lambda_min,lambda_max
INTEGER		::	ncell_min,ncell_max,ncell,c
! Parameters for displaying results and to discretize the posterior
INTEGER 	::	win	! Average window for monitoring convergence
INTEGER 	::	nvt,nvp	! Number of velocity pixels in (latitude, longitude)
INTEGER 	::	nvd,nvn	! Number of pixels in velocity and data noise
! Variables for Delaunay triangulation (construction of Voronoi cells)
INTEGER		:: nt_max,nv_max,nmax,ptsok
REAL(KIND=ii10)	:: eps = 0.000001
INTEGER		:: nt,walk
INTEGER, DIMENSION(:,:), ALLOCATABLE		::	vertices,neighbour
INTEGER, DIMENSION(:), ALLOCATABLE		::	nnn
INTEGER, DIMENSION(:), ALLOCATABLE		::	nnlist,ntwork
INTEGER, DIMENSION(:), ALLOCATABLE		::	worki1,worki2,worki3
LOGICAL, DIMENSION(:), ALLOCATABLE		::	ldummy
! Other parameters
INTEGER		::	ppf,ppl,node,maxx,oilm,nrays,i,j,ii,jj,iii,jjj,nvtr,nvpr,p,pp,nr,ind,ount,w,step,acc,nmci
REAL(KIND=ii10)	::	lat,long,p1,p2,latr,longr,u,ua,ub,misfit,latnew,longnew,velnew
REAL(KIND=ii10)	::	latp,longp,latminp,latmaxp,longminp,longmaxp
! Allocatable arrays
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE	::	Voro,sigma
REAL(KIND=ii10), DIMENSION(:), ALLOCATABLE	::	aa,bb,lambda
INTEGER, DIMENSION(:), ALLOCATABLE		::	EvidenceS	! Distribution on number of cells
INTEGER, DIMENSION(:,:), ALLOCATABLE		::	ML_sigmaS	! Distribution on data noise
INTEGER, DIMENSION(:,:,:), ALLOCATABLE		::	postS		! Distribution on velocity
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE	::	AVES,VARS	! Average and error velocity map
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE	::	maxmap,median	! Max and median velocity map
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE	:: 	area		! Area of pixels for node density
INTEGER, DIMENSION(:), ALLOCATABLE 		::	nfirsts,nlasts	! First and last number of cells
REAL(KIND=ii10), DIMENSION(:,:,:), ALLOCATABLE 	::	firsts,lasts	! First and last model of the chain
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 	::	enspS		! Ensemble of Voronoi centre density
REAL(KIND=ii10), DIMENSION(:,:), ALLOCATABLE 	::	firstpS,lastpS	! First and last Voronoi centre density
REAL(KIND=ii10), DIMENSION(:), ALLOCATABLE 	::	WA,AVW		! Windows for averaging
INTEGER, DIMENSION(:), ALLOCATABLE		::	badmc		! Bad Markov chains (to be ignored in processing)


! Read input file procsamples.in
OPEN(UNIT=20,FILE='procsamples.in',STATUS='old')
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)sampfile
READ(20,1)ncelfile
READ(20,1)sigfile
READ(20,1)misffile
READ(20,1)aratfile
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)burn_in
READ(20,*)thinn
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)win
READ(20,*)nvt,nvp
READ(20,*)latp,longp
READ(20,*)nvd,nvn
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)aveout
READ(20,1)stdevout
READ(20,1)medout
READ(20,1)maxout
READ(20,1)ncellsout
READ(20,1)misfout
READ(20,1)evidout
READ(20,1)noiseout
READ(20,1)postout
READ(20,1)voroout
READ(20,1)nodeout
READ(20,*)oilm
READ(20,1)firstout
READ(20,1)lastout
READ(20,1)cdum
READ(20,1)cdum
READ(20,1)cdum
READ(20,*)nmci
IF (nmci.GT.0) THEN
	ALLOCATE(badmc(nmci))
	DO ii=1,nmci
		READ(20,*)badmc(ii)
	END DO
END IF
1   FORMAT(a30)
CLOSE(20)

! Read parameters from file
OPEN(UNIT=84,FILE='parameters.dat',FORM='unformatted',STATUS='old')
READ(84)nbproc
!nbproc=32
READ(84)sampletotal
READ(84)ncell_min,ncell_max
READ(84)mean,theta
READ(84)latmin,latmax
READ(84)longmin,longmax
READ(84)nds
READ(84)sigdep
READ(84)sigma_min,sigma_max
IF (sigma_max.EQ.sigma_min) sigma_max=2*sigma_min
READ(84)aa_min,aa_max
READ(84)bb_min,bb_max
READ(84)lambda_min,lambda_max
READ(84)nt_max,nv_max,nmax
READ(84)nvtr,nvpr
READ(84)latr,longr
CLOSE(84)

! The number of Markov chains is equal to the total number of chains (processors) minus those that are 'bad'
IF (nmci.GE.nbproc) THEN
	WRITE(*,*)'The number of bad chains is equal/greater than the number of Markov chains available - program is terminating!!!'
	STOP
END IF
nchains = nbproc-nmci

lat = (latmax - latmin)/REAL(nvt-1,KIND=ii10)
long = (longmax - longmin)/REAL(nvp-1,KIND=ii10)

! Allocate arrays
!---------------------------------------------
! Arrays for Delaunay triangulation
ALLOCATE(vertices(3,nt_max))
ALLOCATE(neighbour(3,nt_max))
ALLOCATE(nnn(ncell_max+1))
ALLOCATE(nnlist(nmax))
ALLOCATE(ntwork(nmax))
ALLOCATE(worki1(nv_max))
ALLOCATE(worki2(nv_max))
ALLOCATE(worki3(nv_max))
ALLOCATE(ldummy(ncell_max))
vertices=0
neighbour=0
nnn=0
nnlist=0
ntwork=0
worki1=0
worki2=0
worki3=0
ldummy=0
!---------------------------------------------
! Arrays of samples
ALLOCATE(Voro(3,ncell_max))
IF (sigdep.EQ.0) THEN
	ALLOCATE(sigma(1,nds))
	sigma=0
ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
	ALLOCATE(aa(nds))
	ALLOCATE(bb(nds))
	aa=0
	bb=0
ELSE IF (sigdep.EQ.3) THEN
	ALLOCATE(lambda(nds))
	lambda=0
END IF
ALLOCATE(nfirsts(nbproc))
ALLOCATE(nlasts(nbproc))
ALLOCATE(firsts(3,ncell_max,nbproc))
ALLOCATE(lasts(3,ncell_max,nbproc))
!---------------------------------------------
! Arrays for posterior distributions and maps
IF(sigdep.EQ.0.OR.sigdep.EQ.3) THEN
	ALLOCATE(ML_sigmaS(nds,nvn))
ELSE IF(sigdep.EQ.1.OR.sigdep.EQ.2) THEN
	ALLOCATE(ML_sigmaS(2*nds,nvn))
ELSE IF(sigdep.EQ.4) THEN
	ALLOCATE(ML_sigmaS(1,nvn))
END IF
ML_sigmaS=0
ALLOCATE(AVES(nvp,nvt))
ALLOCATE(VARS(nvp,nvt))
ALLOCATE(enspS(nvp,nvt))
ALLOCATE(firstpS(nvp,nvt))
ALLOCATE(lastpS(nvp,nvt))
ALLOCATE(area(nvp,nvt))
AVES=0
VARS=0
enspS=0
firstpS=0
lastpS=0
DO ii=1,nvp
	longminp = longmin + (ii-1)*long - longp/2
	IF (longminp.LT.longmin) longminp=longmin
	longmaxp = longminp + longp
	IF (longmaxp.GT.longmax) longmaxp=longmax
	DO jj=1,nvt
		latminp = latmin + (jj-1)*lat - latp/2
		IF (latminp.LT.latmin) latminp=latmin
		latmaxp = latminp + latp
		IF (latmaxp.GT.latmax) latmaxp=latmax
		area(ii,jj)=ABS(latmaxp-latminp)*ABS(longmaxp-longminp)
	END DO
END DO
ALLOCATE(EvidenceS(ncell_max))
EvidenceS=0
ALLOCATE(postS(nvp,nvt,nvd))
postS=0
!---------------------------------------------

! Loop over all chains (i.e. over all processors on which the chains were run)
WRITE(*,*)'Looping over chains...'
DO rank=1,nbproc
	
	IF (ANY(badmc.EQ.rank)) THEN
		WRITE(*,*)'Chain',rank,'ignored'
		CYCLE
	END IF
	
	WRITE(*,*)'Processing chain',rank,'of',nbproc
	
	rr=rank
	WRITE(rankc,*)rank
	rankc = ADJUSTL(rankc)
	! Open files
	OPEN(UNIT=50,FILE=trim(sampfile) // '.out' // trim(rankc),STATUS='old',FORM='unformatted')
	OPEN(UNIT=51,FILE=trim(ncelfile) // '.out' // trim(rankc),STATUS='old',FORM='formatted')
	IF (sigdep.NE.4) OPEN(UNIT=52,FILE=trim(sigfile) // '.out' // trim(rankc),STATUS='old',FORM='unformatted')
	IF (sigdep.EQ.4.AND.rr.EQ.1) OPEN(UNIT=52,FILE=trim(sigfile) // '.out',STATUS='old',FORM='unformatted')
	
	! Read and discard burn-in samples
	ss=0
	DO WHILE (ss.LE.burn_in)
		! If we are at the first sample of the chain...
		IF (ss.EQ.0) THEN
			! ...read and save the number of cells...
			READ(51,'(2I10,2I10)')ount,ncell
			!write(*,*)'Sample',ss,'has',ncell,'cells'
			nfirsts(rr)=ncell
			! ...read and save the Voronoi geometry...
			DO c=1,ncell
				READ(50)firsts(1,c,rr),firsts(2,c,rr),firsts(3,c,rr)
				!write(*,*)firsts(1,c,rr),firsts(2,c,rr),firsts(3,c,rr)
			END DO
			Voro(:,:)=firsts(:,:,rr)
			! ...read and ignore data noise parameters...
			IF (sigdep.NE.4) THEN
				DO nss=1,nds
					READ(52)
				END DO
			END IF
			! ...evaluate Voronoi centre density...
			DO ppf=1,ncell
				p1=firsts(1,ppf,rr)
				p2=firsts(2,ppf,rr)
				DO ii=1,nvp
					longminp = longmin + (ii-1)*long - longp/2
					longmaxp = longminp + longp
					DO jj=1,nvt
						latminp = latmin + (jj-1)*lat - latp/2
						latmaxp = latminp + latp
						IF ((p2.GE.latminp).AND.(p2.LT.latmaxp).AND.(p1.GE.longminp).AND.(p1.LT.longmaxp)) THEN
							firstpS(ii,jj)=firstpS(ii,jj)+1
						END IF
					END DO
				END DO
			END DO
		ELSE
			READ(51,'(2I10)')ount,ncell
		        !write(*,*)'Sample',ss,'has',ncell,'cells'
			READ(50)step,acc,ind,longnew,latnew,velnew
			!write(*,*)step,acc,ind,longnew,latnew,velnew
			IF (step.EQ.1.AND.acc.EQ.1) THEN
				Voro(1,ncell)=longnew
				Voro(2,ncell)=latnew
				Voro(3,ncell)=velnew
			ELSE IF (step.EQ.2.AND.acc.EQ.1) THEN
				Voro(:,ind)=Voro(:,ncell+1)
				Voro(:,ncell+1:ncell_max)=0
			ELSE IF (step.EQ.3.AND.acc.EQ.1) THEN
				Voro(1,ind)=longnew
				Voro(2,ind)=latnew
			ELSE IF (step.EQ.4.AND.acc.EQ.1) THEN
				Voro(3,ind)=velnew
			ELSE IF (step.EQ.5.AND.acc.EQ.1) THEN
				Voro=Voro
			ELSE
				Voro=Voro
			END IF
			IF (sigdep.NE.4) THEN
				DO nss=1,nds
					READ(52)
				END DO
			END IF
		END IF
		ss=ss+1
	END DO
	
	! After burn-in, start collecting samples
	thin=0
	DO WHILE (ss.LE.sampletotal-1)
		!Voro=0
		READ(51,'(2I10)')ount,ncell
		write(*,*)'Sample',ss,'has',ncell,'cells'
		READ(50)step,acc,ind,longnew,latnew,velnew
		!write(*,*)step,acc,ind,longnew,latnew,velnew
		IF (step.EQ.1.AND.acc.EQ.1) THEN
			Voro(1,ncell)=longnew
			Voro(2,ncell)=latnew
			Voro(3,ncell)=velnew
		ELSE IF (step.EQ.2.AND.acc.EQ.1) THEN
			Voro(:,ind)=Voro(:,ncell+1)
			Voro(:,ncell+1:ncell_max)=0
		ELSE IF (step.EQ.3.AND.acc.EQ.1) THEN
			Voro(1,ind)=longnew
			Voro(2,ind)=latnew
		ELSE IF (step.EQ.4.AND.acc.EQ.1) THEN
			Voro(3,ind)=velnew
		ELSE IF (step.EQ.5.AND.acc.EQ.1) THEN
			Voro=Voro
		ELSE
			Voro=Voro
		END IF
		IF (sigdep.NE.4) THEN
			DO nss=1,nds
				IF (sigdep.EQ.0) READ(52)sigma(1,nss)
				IF (sigdep.EQ.1.OR.sigdep.EQ.2) READ(52)aa(nss),bb(nss)
				IF (sigdep.EQ.3) READ(52)lambda(nss)
			END DO
		ELSE IF (sigdep.EQ.4.AND.rr.EQ.1.AND.ss.EQ.sampletotal-1) THEN
			READ(52)nrays
			ALLOCATE(sigma(nrays,1))
			DO nr=1,nrays
				READ(52)sigma(nr,1)
			END DO
		END IF
		
		! Add sample to posterior
		IF (mod(ss,thinn).EQ.0) THEN
			! Number of thinned samples increases by 1
			thin=thin+1
			! Posterior on number of cells
			EvidenceS(ncell)=EvidenceS(ncell)+1
			! Calculate Delaunay triangulation for current sample
			CALL delaun (DBLE(Voro(1:2,:)),ncell,neighbour,vertices,nt,nt_max,&
					worki1,worki2,worki3,eps,nv_max,&
					0,ldummy,0,0,0,ptsok)
			CALL build_nv(ncell,vertices,nt,ncell_max,nmax,&
				neighbour,nnn,nnlist,ntwork)
			node = 1
			! Loop over latitude and lomgitude points to evaluate average and error map
			DO i=1,nvp
				DO j=1,nvt
					p1=longmin+(i-1)*long
					p2=latmin+(j-1)*lat
					CALL find_node2D(DBLE([p1,p2]),node,DBLE(Voro(1:2,:)),nnn,nnlist,walk)
					! Average
					AVES(i,j) = AVES(i,j) + DBLE(Voro(3,node))
					! Standard deviation
					VARS(i,j) = VARS(i,j) + DBLE(Voro(3,node))**2
					jj=ceiling((Voro(3,node)-(mean-theta))*nvd/(2*theta))
					IF (jj.EQ.0) jj=1
					! Posterior on velocity at each point of the grid
					postS(i,j,jj)=postS(i,j,jj)+1
				END DO
			END DO
			! Calculate Voronoi centre density
			DO ii=1,nvp
				longminp = longmin + (ii-1)*long - longp/2
				longmaxp = longminp + longp
				DO jj=1,nvt
					latminp = latmin + (jj-1)*lat - latp/2
					latmaxp = latminp + latp
					DO c=1,ncell
						p1=Voro(1,c)
						p2=Voro(2,c)
						IF ((p2.GE.latminp).AND.(p2.LT.latmaxp).AND.(p1.GE.longminp).AND.(p1.LT.longmaxp)) THEN
							enspS(ii,jj)=enspS(ii,jj)+1
						END IF
					END DO
				END DO
			END DO
			! Posterior on data noise parameters
			IF (sigdep.EQ.0) THEN
				DO nss=1,nds
					ii=ceiling((sigma(1,nss)-sigma_min)*nvn/(sigma_max-sigma_min))
					IF (ii.EQ.0) ii=1
					ML_sigmaS(nss,ii) = ML_sigmaS(nss,ii)+1
				END DO
			ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
				DO nss=1,nds
					ii=ceiling((aa(nss)-aa_min)*nvn/(aa_max-aa_min))
					IF (ii.EQ.0) ii=1	
					ML_sigmaS(nss,ii) = ML_sigmaS(nss,ii)+1
					ii=ceiling((bb(nss)-bb_min)*nvn/(bb_max-bb_min))
					IF (ii.EQ.0) ii=1	
					ML_sigmaS(nss+nds,ii) = ML_sigmaS(nss+nds,ii)+1
				END DO
			ELSE IF (sigdep.EQ.3) THEN
				DO nss=1,nds
					ii=ceiling((lambda(nss)-lambda_min)*nvn/(lambda_max-lambda_min))
					IF (ii.EQ.0) ii=1	
					ML_sigmaS(nss,ii) = ML_sigmaS(nss,ii)+1
				END DO
			ELSE IF (sigdep.EQ.4.AND.rr.EQ.1.AND.ss.EQ.sampletotal-1) THEN
				DO nr=1,nrays
					ii=ceiling((sigma(nr,1)-sigma_min)*nvn/(sigma_max-sigma_min))
					IF (ii.EQ.0) ii=1
					IF (ii.GT.nvn) ii=nvn
					ML_sigmaS(1,ii) = ML_sigmaS(1,ii)+1
				END DO
			END IF
			
		END IF
		! Evaluate Voronoi centre density for last model
		IF (ss.EQ.sampletotal-1) THEN
			nlasts(rr)=ncell
			lasts(:,:,rr)=Voro
			DO ppl=1,ncell
				p1=lasts(1,ppl,rr)
				p2=lasts(2,ppl,rr)
				DO ii=1,nvp
					longminp = longmin + (ii-1)*long - longp/2
					longmaxp = longminp + longp
					DO jj=1,nvt
						latminp = latmin + (jj-1)*lat - latp/2
						latmaxp = latminp + latp
						IF ((p2.GE.latminp).AND.(p2.LT.latmaxp).AND.(p1.GE.longminp).AND.(p1.LT.longmaxp)) THEN
							lastpS(ii,jj)=lastpS(ii,jj)+1
						END IF
					END DO
				END DO
			END DO
		END IF
		ss=ss+1
	END DO
	CLOSE(50)
	CLOSE(51)
	IF (sigdep.NE.4.OR.(sigdep.EQ.4.AND.rr.EQ.1)) CLOSE(52)
END DO

AVES=AVES/DBLE(thin*nchains)
VARS=VARS/DBLE(thin*nchains)
VARS = sqrt(VARS - AVES**2)
DO i=1,nvp
	DO j=1,nvt
		IF (ISNAN(VARS(i,j))) VARS(i,j)=0
	END DO
END DO
enspS=enspS/REAL(thin*nchains,KIND=ii10)

!----------------------------------------------------------------------!
!!!!!!!!!!                     SAVE FILES                     !!!!!!!!!!
!----------------------------------------------------------------------!
! Save files of first and last model
WRITE(*,*)'Saving first and last Voronoi cell files...'
voroout = ADJUSTL(voroout)
voroini = trim(voroout) // '.ini'
vorofin = trim(voroout) // '.fin'
OPEN(UNIT=81,FILE=voroini,STATUS='replace')
OPEN(UNIT=83,FILE=vorofin,STATUS='replace')
DO rr=1,nbproc
	
	IF (ANY(badmc.EQ.rr)) THEN
		WRITE(*,*)'   ...chain',rr,'ignored...'
		CYCLE
	END IF
	
	! Save files containing all Voronoi cells
	WRITE(81,*)nfirsts(rr)
	WRITE(83,*)nlasts(rr)
	DO c=1,nfirsts(rr)
			WRITE(81,*)firsts(2,c,rr),firsts(1,c,rr),firsts(3,c,rr)
	END DO
	WRITE(81,'(1X)')
	DO c=1,nlasts(rr)
		WRITE(83,*)lasts(2,c,rr),lasts(1,c,rr),lasts(3,c,rr)
	END DO
	WRITE(83,'(1X)')
	
	! Depending on choice made in input file, save single file of 
	! first and last model of one chain
	IF (rr.EQ.oilm.OR.oilm.EQ.-1) THEN
		rank=rr
		WRITE(rankc,*)rank
		rankc = ADJUSTL(rankc)
		
		! First model
		CALL delaun (DBLE(firsts(1:2,:,rr)),nfirsts(rr),neighbour,vertices,nt,nt_max,&
				worki1,worki2,worki3,eps,nv_max,&
				0,ldummy,0,0,0,ptsok)
		CALL build_nv(nfirsts(rr),vertices,nt,ncell_max,nmax,&
			neighbour,nnn,nnlist,ntwork)
		WRITE(*,*)'Saving geometry of first model of chain',rr
		firstout = ADJUSTL(firstout)
		firstvtx = trim(firstout) // '_p' // trim(rankc) // '.vtx'
		firstvor = trim(firstout) // '_p' // trim(rankc) // '.vor'
		OPEN(UNIT=60,FILE=firstvtx,STATUS='replace')
		WRITE(60,*)nvtr,nvpr
		WRITE(60,'(3f14.8)')latmax,longmin
		WRITE(60,'(3f14.8)')latr,longr
		WRITE(60,'(1X)')
		OPEN(UNIT=65,FILE=firstvor,STATUS='replace')
		WRITE(65,*)nvtr,nvpr
		WRITE(65,'(3f14.8)')latmax,longmin
		WRITE(65,'(3f14.8)')latr,longr
		WRITE(65,'(1X)')
		WRITE(65,*)nfirsts(rr)
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
				CALL find_node2D(DBLE([p1,p2]),node,DBLE(firsts(1:2,:,rr)),nnn,nnlist,walk)
				WRITE(60,*)firsts(3,node,rr)
			END DO
			WRITE(60,'(1X)')
		END DO
		CLOSE(60)	
		DO c=1,nfirsts(rr)
			WRITE(65,*)firsts(2,c,rr),firsts(1,c,rr),firsts(3,c,rr)
		END DO
		CLOSE(65) ! close the file
		
		! Last model
		CALL delaun (DBLE(lasts(1:2,:,rr)),nlasts(rr),neighbour,vertices,nt,nt_max,&
				worki1,worki2,worki3,eps,nv_max,&
				0,ldummy,0,0,0,ptsok)
		CALL build_nv(nlasts(rr),vertices,nt,ncell_max,nmax,&
			neighbour,nnn,nnlist,ntwork)
		WRITE(*,*)'Saving geometry of last model of chain',rr
		lastout = ADJUSTL(lastout)
		lastvtx = trim(lastout) // '_p' // trim(rankc) // '.vtx'
		lastvor = trim(lastout) // '_p' // trim(rankc) // '.vor'
		OPEN(UNIT=60,FILE=lastvtx,STATUS='replace')
		WRITE(60,*)nvtr,nvpr
		WRITE(60,'(3f14.8)')latmax,longmin
		WRITE(60,'(3f14.8)')latr,longr
		WRITE(60,'(1X)')
		OPEN(UNIT=65,FILE=lastvor,STATUS='replace')
		WRITE(65,*)nvtr,nvpr
		WRITE(65,'(3f14.8)')latmax,longmin
		WRITE(65,'(3f14.8)')latr,longr
		WRITE(65,'(1X)')
		WRITE(65,*)nlasts(rr)
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
				CALL find_node2D(DBLE([p1,p2]),node,DBLE(lasts(1:2,:,rr)),nnn,nnlist,walk)
				WRITE(60,*)lasts(3,node,rr)
			END DO
			WRITE(60,'(1X)')
		END DO
		CLOSE(60)	
		DO c=1,nlasts(rr)
			WRITE(65,*)lasts(2,c,rr),lasts(1,c,rr),lasts(3,c,rr)
		END DO
		CLOSE(65) ! close the file
	END IF
END DO
CLOSE(81) ! close the file
CLOSE(83) ! close the file
DEALLOCATE(firsts,lasts,nfirsts,nlasts)
WRITE(*,*)'...done!'
!----------------------------------------------------------------------!
! Save files of Voronoi centre densities
WRITE(*,*)'Saving files of Voronoi centre densities...'
nodeout = ADJUSTL(nodeout)
nodeens = trim(nodeout) // '.ens'
nodeini = trim(nodeout) // '.ini'
nodefin = trim(nodeout) // '.fin'
OPEN(UNIT=86,FILE=nodeens,STATUS='replace')
WRITE(86,*)nvt,nvp
WRITE(86,'(3f14.8)')latmax,longmin
WRITE(86,'(3f14.8)')lat,long
WRITE(86,'(1X)')
OPEN(UNIT=81,FILE=nodeini,STATUS='replace')
WRITE(81,*)nvt,nvp
WRITE(81,'(3f14.8)')latmax,longmin
WRITE(81,'(3f14.8)')lat,long
WRITE(81,'(1X)')
OPEN(UNIT=83,FILE=nodefin,STATUS='replace')
WRITE(83,*)nvt,nvp
WRITE(83,'(3f14.8)')latmax,longmin
WRITE(83,'(3f14.8)')lat,long
WRITE(83,'(1X)')
DO ii=0,nvp+1
	DO jj=0,nvt+1
		IF (ii.EQ.0.OR.ii.EQ.nvp+1) THEN
			WRITE(86,*) 0
		ELSE IF (jj.EQ.0.OR.jj.EQ.nvt+1) THEN
			WRITE(86,*) 0
		ELSE
			WRITE(86,*)enspS(ii,nvt+1-jj)/area(ii,nvt+1-jj)
		END IF
		
		IF (ii.EQ.0.OR.ii.EQ.nvp+1) THEN
			WRITE(81,*) 0
		ELSE IF (jj.EQ.0.OR.jj.EQ.nvt+1) THEN
			WRITE(81,*) 0
		ELSE
			WRITE(81,*)firstpS(ii,nvt+1-jj)/(area(ii,nvt+1-jj)*(nchains))
		END IF
		
		IF (ii.EQ.0.OR.ii.EQ.nvp+1) THEN
			WRITE(83,*) 0
		ELSE IF (jj.EQ.0.OR.jj.EQ.nvt+1) THEN
			WRITE(83,*) 0
		ELSE
			WRITE(83,*)lastpS(ii,nvt+1-jj)/(area(ii,nvt+1-jj)*(nchains))
		END IF
	END DO
	WRITE(86,'(1X)')
	WRITE(81,'(1X)')
	WRITE(83,'(1X)')
END DO
CLOSE(86) ! close the file
CLOSE(81) ! close the file
CLOSE(83) ! close the file
DEALLOCATE(firstpS,lastpS,area,enspS)
WRITE(*,*)'...done!'
!----------------------------------------------------------------------!
! Save files of average and error map
WRITE(*,*)'Saving average and error maps...'
! Plot the Average map
OPEN(UNIT=52,FILE=aveout,STATUS='replace')
WRITE(52,*)nvt,nvp
WRITE(52,'(3f14.8)')latmax,longmin
WRITE(52,'(3f14.8)')lat,long
WRITE(52,'(1X)')
! Plot the Variance map
OPEN(UNIT=53,FILE=stdevout,STATUS='replace')
WRITE(53,*)nvt,nvp
WRITE(53,'(3f14.8)')latmax,longmin
WRITE(53,'(3f14.8)')lat,long
WRITE(53,'(1X)')
DO i=0,nvp+1
	DO j=0,nvt+1
		ii=i
		jj=j
		IF (i.EQ.0)   ii=1
		IF (j.EQ.0)   jj=1
		IF (i.EQ.nvp+1) ii=nvp
		IF (j.EQ.nvt+1) jj=nvt
		WRITE(52,*)AVES(ii,nvt+1-jj)
		WRITE(53,*)VARS(ii,nvt+1-jj)
	END DO
	WRITE(52,'(1X)')
	WRITE(53,'(1X)')
END DO
CLOSE(52)	! close the file 
CLOSE(53)	! close the file
DEALLOCATE(AVES,VARS)
WRITE(*,*)'...done!'
!----------------------------------------------------------------------!
! Evaluate and save map of median and max, and file of posterior 
! velocity distributions
! Take the maximum and median of the posterior at each pixel
WRITE(*,*)'Evaluating maximum and median map...'
ALLOCATE(maxmap(nvp,nvt))
ALLOCATE(median(nvp,nvt))
maxmap=0
median=0
OPEN(UNIT=85,FILE=postout,FORM='unformatted',STATUS='replace')
WRITE(85)nvt,nvp
WRITE(85)latmax,longmin
WRITE(85)lat,long
WRITE(85)nvd
DO jj=1,nvd
	WRITE(85)mean-theta+(jj-0.5)*2*theta/nvd
END DO
DO i=1,nvp
	DO j=1,nvt
		maxx=0
		! Take the max of postS 
		ind=1
		p=0
		DO jj=1,nvd
			p=p+postS(i,j,jj)
			IF (postS(i,j,jj)>maxx) THEN
				ind=jj
				maxx=postS(i,j,jj)
			END IF
		END DO
		maxmap(i,j)=mean-theta+(ind-0.5)*2*theta/nvd
		! Take the median of postS
		pp=0
		DO jj=1,nvd
			pp=pp+postS(i,j,jj)
			IF (pp>(p/2)) THEN
				ind=jj
				median(i,j)=mean-theta+(ind-0.5)*2*theta/nvd
				EXIT
			END IF
		END DO
		DO jj=1,nvd
			WRITE(85)postS(i,j,jj)
		END DO
	END DO
END DO
CLOSE(85)
WRITE(*,*)'...done!'
DEALLOCATE(postS)
WRITE(*,*)'Saving maximum and median map...'
! Plot the Median map
OPEN(UNIT=54,FILE=medout,STATUS='replace')
WRITE(54,*)nvt,nvp
WRITE(54,'(3f14.8)')latmax,longmin
WRITE(54,'(3f14.8)')lat,long
WRITE(54,'(1X)')
! Plot the Maximum map
OPEN(UNIT=5,FILE=maxout,STATUS='replace')
WRITE(5,*)nvt,nvp
WRITE(5,'(3f14.8)')latmax,longmin
WRITE(5,'(3f14.8)')lat,long
WRITE(5,'(1X)')
DO i=0,nvp+1
	DO j=0,nvt+1
		ii=i
		jj=j
		IF (i.EQ.0)   ii=1
		IF (j.EQ.0)   jj=1
		IF (i.EQ.nvp+1) ii=nvp
		IF (j.EQ.nvt+1) jj=nvt
		WRITE(5,*)maxmap(ii,nvt+1-jj)
		WRITE(54,*)median(ii,nvt+1-jj)
	END DO
	WRITE(5,'(1X)')
	WRITE(54,'(1X)')
END DO
CLOSE(5)	! close the file 
CLOSE(54)	! close the file
DEALLOCATE(maxmap,median)
WRITE(*,*)'...done!'
!----------------------------------------------------------------------!
! Save file of posterior distribution on the number of cells
WRITE(*,*)'Saving posterior on number of cells...'
OPEN(UNIT=54,FILE=evidout,STATUS='replace')
DO i=1,ncell_max
	WRITE(54,*)i,EvidenceS(i)
END DO
CLOSE(54)! close the file 
DEALLOCATE(EvidenceS)
WRITE(*,*)'...done!'
!----------------------------------------------------------------------!
! Save file of posterior distribution on data noise
WRITE(*,*)'Saving posterior on data noise...'
OPEN(UNIT=59,FILE=noiseout,STATUS='replace')
IF (sigdep.EQ.0) THEN
	DO i=1,nvn
		u=sigma_min+(i-0.5)*(sigma_max-sigma_min)/nvn
		WRITE(59,'(1f14.8)',advance='no')u
		DO nss=1,nds
			IF (nss.NE.nds) WRITE(59,'(1I20)',advance='no')ML_sigmaS(nss,i)
			IF (nss.EQ.nds) WRITE(59,'(1I20)')ML_sigmaS(nss,i)
		END DO
	END DO
ELSE IF (sigdep.EQ.1.OR.sigdep.EQ.2) THEN
	DO i=1,nvn
		ua=aa_min+(i-0.5)*(aa_max-aa_min)/nvn
		ub=bb_min+(i-0.5)*(bb_max-bb_min)/nvn
		WRITE(59,'(1f14.8)',advance='no')ua
		DO nss=1,nds
			WRITE(59,'(1I20)',advance='no')ML_sigmaS(nss,i)
		END DO
		WRITE(59,'(1f14.8)',advance='no')ub
		DO nss=1,nds
			IF (nss.NE.nds) WRITE(59,'(1I20)',advance='no')ML_sigmaS(nss+nds,i)
			IF (nss.EQ.nds) WRITE(59,'(1I20)')ML_sigmaS(nss+nds,i)
		END DO
	END DO
ELSE IF (sigdep.EQ.3) THEN
	DO i=1,nvn
		u=lambda_min+(i-0.5)*(lambda_max-lambda_min)/nvn
		WRITE(59,'(1f14.8)',advance='no')u
		DO nss=1,nds
			IF (nss.NE.nds) WRITE(59,'(1I20)',advance='no')ML_sigmaS(nss,i)
			IF (nss.EQ.nds) WRITE(59,'(1I20)')ML_sigmaS(nss,i)
		END DO
	END DO
ELSE IF (sigdep.EQ.4) THEN
	DO i=1,nvn
		u=sigma_min+(i-0.5)*(sigma_max-sigma_min)/nvn
		WRITE(59,'(1f14.8,1I20)')u,ML_sigmaS(1,i)
	END DO
END IF	
CLOSE(59)! close the file with data noise
DEALLOCATE(ML_sigmaS)
WRITE(*,*)'...done!'
!----------------------------------------------------------------------!

! Deallocate arrays
DEALLOCATE(vertices)
DEALLOCATE(neighbour)
DEALLOCATE(nnn)
DEALLOCATE(nnlist)
DEALLOCATE(ntwork)
DEALLOCATE(worki1)
DEALLOCATE(worki2)
DEALLOCATE(worki3)
DEALLOCATE(ldummy)
DEALLOCATE(Voro)
IF (ALLOCATED(sigma)) DEALLOCATE(sigma)
IF (ALLOCATED(aa)) DEALLOCATE(aa)
IF (ALLOCATED(bb)) DEALLOCATE(bb)
IF (ALLOCATED(lambda)) DEALLOCATE(lambda)

!----------------------------------------------------------------------!
! Calculate window-averaged misfits and number of cells
! (to monitor convergence)
ALLOCATE(WA(win))
ALLOCATE(AVW(sampletotal-1))

! Number of cells
WRITE(*,*)'Saving window-averaged number of cells...'
WA=0
AVW=0
w=0
DO rr=1,nbproc
	
	IF (ANY(badmc.EQ.rr)) THEN
		WRITE(*,*)'   ...chain',rr,'ignored...'
		CYCLE
	END IF
	
	rank=rr
	WRITE(rankc,*)rank
	rankc = ADJUSTL(rankc)
	! Open files
	OPEN(UNIT=50,FILE=trim(ncelfile) // '.out' // trim(rankc),STATUS='old',FORM='formatted')
	READ(50,*)
	DO ss=1,sampletotal-1
		READ(50,'(2I10)')ount,ncell
		AVW(ss)=AVW(ss)+REAL(ncell,KIND=ii10)
	END DO
	CLOSE(50)
END DO
AVW=AVW/(nchains)
OPEN(UNIT=44,FILE=ncellsout,STATUS='replace')
DO ss=1,sampletotal-1
	IF (ss.LT.win) THEN
		w=w+1
		WA(w)=AVW(ss)
	ELSE IF (ss.EQ.win) THEN
		w=w+1
		WA(w)=AVW(ss)
		WRITE(44,*)SUM(WA)/REAL(win,KIND=ii10)
	ELSE IF (ss.GT.win) THEN 
		DO i=1,win-1
			WA(i)=WA(i+1)
		END DO
		WA(win)=AVW(ss)
		WRITE(44,*)SUM(WA)/REAL(win,KIND=ii10)
	END IF
END DO
CLOSE(44)

! Misfit
WRITE(*,*)'...and misfit...'
WA=0
AVW=0
w=0
DO rr=1,nbproc
	
	IF (ANY(badmc.EQ.rr)) THEN
		WRITE(*,*)'   ...chain',rr,'ignored...'
		CYCLE
	END IF
	
	rank=rr
	WRITE(rankc,*)rank
	rankc = ADJUSTL(rankc)
	! Open files
	OPEN(UNIT=50,FILE=trim(misffile) // '.out' // trim(rankc),STATUS='old',FORM='formatted')
	READ(50,*)
	DO ss=1,sampletotal-1
		READ(50,'(1I10,1F20.6)')ount,misfit
		AVW(ss)=AVW(ss)+misfit
	END DO
	CLOSE(50)
END DO
AVW=AVW/(nchains)
OPEN(UNIT=44,FILE=misfout,STATUS='replace')
DO ss=1,sampletotal-1
	IF (ss.LT.win) THEN
		w=w+1
		WA(w)=AVW(ss)
	ELSE IF (ss.EQ.win) THEN
		w=w+1
		WA(w)=AVW(ss)
		WRITE(44,*)SUM(WA)/REAL(win,KIND=ii10)
	ELSE IF (ss.GT.win) THEN 
		DO i=1,win-1
			WA(i)=WA(i+1)
		END DO
		WA(win)=AVW(ss)
		WRITE(44,*)SUM(WA)/REAL(win,KIND=ii10)
	END IF
END DO
CLOSE(44)
DEALLOCATE(WA,AVW)
WRITE(*,*)'...done!'


END PROGRAM
