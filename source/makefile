# COMPILATION COMMANDS REWRITTEN FOR RUNNING AT SGR CLUSTER, CAMBRIDGE, UK 
# 
# THIS VERSION HAS BEEN USED FOR EXPERIMENTS FOR THE PARALLELIZATION OF THE RAY ROUTINES
#------------------------------------------------------------------------
#
#		Makefile for RJ-MCMC code driver program 
#
#  Includes:
#	
#	driver		compiles mksamples and procsamples programs only
#	nn, qhull 	contain the subroutines necessary to compute the
#			Delaunay triangulation (Voronoi cells) and 
#			find to which Voronoi cell a certain point on a
#			grid belongs to
#	MT		contains the subroutine that generates random
#			numbers using the Mersenne-Twister algorithm
#	fm2d		contains the subroutines necessary to model 
#			raypaths (an adapted version of the fm2dss code 
#			in FMST) and update raypath geometries in the
#			rj_tomo code
#	EGcodes		contains the code voroslice (which produces a file
#			of Voronoi cell velocities compatible with GMT -
#			adapted from tslicess in FMST) and the code 
#			fm2dssRJ (which produces a raypath input file 
#			compatible with rj_tomo - adapted from fm2dss
#			in FMST)
#	fmst		contains the FMST code (bug fixed in fm2dss)
#	input		all input files
#
#------------------------------------------------------------------------


# Initialize the environment	
# Use Faramir compiler with MPI 

# Here we have two options:  use GNU compiler or Intel compiler for Fortran and Fortran-MPI 

# Run this module for GNU compiler.
# module load openmpi/1.4/gnu_4.1.2

# Run this module for Intel compilers
# module load openmpi/1.4/intel_11.1
# module load intel/11.1

#F90 = gfortran  #GNU Fortran
#F90mpi = mpif90 #mpif90 GNU Fortran
#CC = gcc	 # GNU cc 

# Intel compilers

# Remove comment for static library compilation 
# LDFLAGS=–static –disable-shared –enable-static

F90 = ifort
F90mpi = mpif90  #mpi ifort
CCmpi = mpicc
CC = icc 

rj_dir=../RJ_MCMC_TOMO
bin_dir = $(rj_dir)/bin
plot_dir = $(bin_dir)/gmtplot

#-----------------------------------------------------------------------
# HW: IMPORTANT NOTES:  
# Line LDFLAGS=–static –disable-shared –enable-static can be inserted above for static library compilation.
# Line -fast -O5 xprofile=use -unroll=8 was included in rj_tomo_mksamples.f90 and fm2dssRJu_traveltime.f90
# Remove for normal compilation


all  : nn qhull fm2d EGcodes driver fmst input

all16  : nn qhull fm2d EGcodes driver16 fmst input

exec : nn qhull fm2d EGcodes driver fmst

driver : fm2d MT rjtomo input

other : EGcodes fmst

driver16: fm2d MT rjtomo16 input


# TESTS

driverold : fm2d MT rjtomoa input

old : nn qhull fm2d EGcodes driverold fmst input

hello: hellow


rjtomo  :: 
	$(F90mpi) -o mks_fastrays rj_tomo_mksamples_fastrays.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90) -o procsamples rj_tomo_procsamples.f90 -traceback -O1 -L./ -lnn2d -lqh 
	\rm globalp.mod traveltime.mod fm2dout.mod fm2dssvars.mod mt19937.mod  
	\mv mks_fastrays $(bin_dir)
	\mv procsamples $(bin_dir)

rjtomo16  :: 
	$(F90mpi) -o mks_fastrays1 rj_tomo_mksamples_fastrays1.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays2 rj_tomo_mksamples_fastrays2.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays3 rj_tomo_mksamples_fastrays3.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays4 rj_tomo_mksamples_fastrays4.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays5 rj_tomo_mksamples_fastrays5.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays6 rj_tomo_mksamples_fastrays6.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays7 rj_tomo_mksamples_fastrays7.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays8 rj_tomo_mksamples_fastrays8.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays9 rj_tomo_mksamples_fastrays9.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays10 rj_tomo_mksamples_fastrays10.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays11 rj_tomo_mksamples_fastrays11.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays12 rj_tomo_mksamples_fastrays12.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays13 rj_tomo_mksamples_fastrays13.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays14 rj_tomo_mksamples_fastrays14.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays15 rj_tomo_mksamples_fastrays15.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90mpi) -o mks_fastrays16 rj_tomo_mksamples_fastrays16.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90) -o procsamples rj_tomo_procsamples.f90 -traceback -O1 -L./ -lnn2d -lqh 
	\rm globalp.mod traveltime.mod fm2dout.mod fm2dssvars.mod mt19937.mod  
	\mv mks_fastrays1 $(bin_dir)
	\mv mks_fastrays2 $(bin_dir)
	\mv mks_fastrays3 $(bin_dir)
	\mv mks_fastrays4 $(bin_dir)
	\mv mks_fastrays5 $(bin_dir)
	\mv mks_fastrays6 $(bin_dir)
	\mv mks_fastrays7 $(bin_dir)
	\mv mks_fastrays8 $(bin_dir)
	\mv mks_fastrays9 $(bin_dir)
	\mv mks_fastrays10 $(bin_dir)
	\mv mks_fastrays11 $(bin_dir)
	\mv mks_fastrays12 $(bin_dir)
	\mv mks_fastrays13 $(bin_dir)
	\mv mks_fastrays14 $(bin_dir)
	\mv mks_fastrays15 $(bin_dir)
	\mv mks_fastrays16 $(bin_dir)
	\mv procsamples $(bin_dir)
	
MT	::
	$(F90) -c MT/mt19937.f90 
	ar -r libranmt.a mt19937.o
	\rm mt19937.o

nn	:: 
	$(F90) -assume no2underscore -c nn/del_sub.f 
	$(F90) -assume no2underscore -c nn/nn.f 
	$(F90) -assume no2underscore -c nn/delaun.f 
	$(F90) -assume no2underscore -c nn/nnplot.f 
	$(CC) -c nn/stack.c 
	$(CC) -c nn/stackpair.c 
	$(CC) -c nn/volume.c 
	$(CC) -c nn/utils.c 
	ar -r libnn2d.a  del_sub.o nn.o delaun.o nnplot.o stack.o stackpair.o volume.o utils.o
	\rm del_sub.o nn.o delaun.o nnplot.o stack.o stackpair.o volume.o utils.o


qhull	:: 
	$(CC) -c qhull/geom.c 
	$(CC) -c qhull/globals.c 
	$(CC) -c qhull/io.c 
	$(CC) -c qhull/poly.c 
	$(CC) -c qhull/qhull.c 
	$(CC) -c qhull/set.c 
	$(CC) -c qhull/mem.c 
	ar -r libqh.a  geom.o globals.o io.o poly.o qhull.o set.o mem.o
	\rm geom.o globals.o io.o poly.o qhull.o set.o mem.o

fm2d	::
	$(F90) -c fm2d/fm2dssRJu_globalp.f90 
	$(F90) -c fm2d/fm2dssRJu_traveltime.f90
	$(F90) -c fm2d/fm2dssRJu_fm2dout.f90  
	$(F90) -c fm2d/fm2dssRJu_main.f90 
	$(F90) -c fm2d/fm2dssRJu_main_subset.f90 
	ar -r libuprays.a fm2dssRJu_globalp.o fm2dssRJu_traveltime.o fm2dssRJu_fm2dout.o fm2dssRJu_main.o fm2dssRJu_main_subset.o  
	\rm fm2dssRJu_globalp.o fm2dssRJu_traveltime.o fm2dssRJu_fm2dout.o fm2dssRJu_main.o fm2dssRJu_main_subset.o   

EGcodes	::
	$(F90) -o bin2txt EGcodes/bin2txt.f90 
	$(F90) -o fm2dssRJ EGcodes/fm2dssRJ.f90 
	$(F90) -o gauserrdist EGcodes/gauserrdist.f90 
	$(F90) -o gauserr EGcodes/gauserr.f90 
	$(F90) -o mkstraightrays EGcodes/mkstraightrays.f90 
	$(F90) -o raylength EGcodes/raylength.f90 
	$(F90) -o velpost EGcodes/velpost.f90 
	$(F90) -o vorogrid EGcodes/vorogrid.f90 -L./ -lnn2d -lqh 
	$(F90) -o vororay EGcodes/vororay.f90 -L./ -lnn2d -lqh -luprays 
	$(F90) -o voroslice EGcodes/voroslice.f90 -L./ -lnn2d -lqh  
	\rm fm2dssvars.mod
	\mv bin2txt fm2dssRJ gauserr gauserrdist mkstraightrays raylength vorogrid vororay $(bin_dir)
	\mv velpost voroslice $(plot_dir)
	
fmst	::
	$(F90) -o fm2dss fmst/fm2dss.f90 
	$(F90) -o grid2dss fmst/grid2dss.f90 
	$(F90) -o misfitss fmst/misfitss.f90 
	$(F90) -o residualss fmst/residualss.f90 
	$(F90) -o resplotss fmst/resplotss.f90 
	$(F90) -o subinvss fmst/subinvss.f90 
	$(F90) -o synthtss fmst/synthtss.f90 
	$(F90) -o tslicess fmst/tslicess.f90 
	\mv fm2dss grid2dss misfitss residualss subinvss synthtss $(bin_dir)
	\cp fmst/ttomoss $(bin_dir)
	\mv resplotss tslicess $(plot_dir)
	
input	:: 
	\cp mksamples.in $(bin_dir)
	\cp procsamples.in $(bin_dir)
	\cp fmst/inputfiles/fm2dss.in $(bin_dir)
	\cp fmst/inputfiles/grid2dss.in $(bin_dir)
	\cp fmst/inputfiles/misfitss.in $(bin_dir)
	\cp fmst/inputfiles/residualss.in $(bin_dir)
	\cp fmst/inputfiles/resplotss.in $(plot_dir)
	\cp fmst/inputfiles/subinvss.in $(bin_dir)
	\cp fmst/inputfiles/synthtss.in $(bin_dir)
	\cp fmst/inputfiles/tslicess.in $(plot_dir)
	\cp fmst/inputfiles/ttomoss.in $(bin_dir)
	\cp EGcodes/inputfiles/bin2txt.in $(bin_dir)
	\cp EGcodes/inputfiles/fm2dssRJ.in $(bin_dir)
	\cp EGcodes/inputfiles/gauserr.in $(bin_dir)
	\cp EGcodes/inputfiles/gauserrdist.in $(bin_dir)
	\cp EGcodes/inputfiles/mkstraightrays.in $(bin_dir)
	\cp EGcodes/inputfiles/raylength.in $(bin_dir)
	\cp EGcodes/inputfiles/velpost.in $(plot_dir)
	\cp EGcodes/inputfiles/vorogrid.in $(bin_dir)
	\cp EGcodes/inputfiles/vororay.in $(bin_dir)
	\cp EGcodes/inputfiles/voroslice.in $(plot_dir)

rjtomoa :: 
	$(F90mpi) -o mksamples rj_tomo_mksamples.f90 -L./ -lnn2d -lqh -luprays -lranmt 
	$(F90) -o procsamples rj_tomo_procsamples.f90 -traceback -O1 -L./ -lnn2d -lqh 
	\rm globalp.mod traveltime.mod fm2dout.mod fm2dssvars.mod mt19937.mod  
	\mv mksamples $(bin_dir)make all
	\mv procsamples $(bin_dir)

hellow  :: 
	$(F90mpi) -o hellowfortranmpi hellowfortranmpi.f90
	\mv hellowfortranmpi $(bin_dir)

clean   ::
	\rm -f *.o *.a *.mod nn/*.o qhull/*.o EGcodes/*.o EGcodes/*.mod fm2d/*.o fm2d/*.mod 

