!!! Output module !!!

MODULE fm2dout

IMPLICIT NONE
!INTEGER, PARAMETER :: iii5=SELECTED_REAL_KIND(5,10)
INTEGER, PARAMETER :: iii10=SELECTED_REAL_KIND(10,100)

INTEGER 					:: 	nru    		!Defined in module fm2dout
INTEGER, DIMENSION(:,:), ALLOCATABLE		::	npoints  	!Defined in module fm2dout
INTEGER, DIMENSION(:,:), ALLOCATABLE		::	npoints_prev
INTEGER, DIMENSION(:,:), ALLOCATABLE		::	npoints_subset
REAL(KIND=iii10), DIMENSION(:,:), ALLOCATABLE	:: 	raycoords  	!Defined in module fm2dout
REAL(KIND=iii10), DIMENSION(:,:), ALLOCATABLE	::	raycoords_subset
REAL(KIND=iii10), DIMENSION(:,:), ALLOCATABLE	::	raycoords_prev

INTEGER						:: 	crazyray  	!Defined in module fm2dout
INTEGER						:: 	crazyray_subset

INTEGER 					:: 	nru_subset
INTEGER 					::	numraycoords, numraycoords_prev, numraycoords_subset


!For testing the parallelization

INTEGER, DIMENSION(:,:), ALLOCATABLE		::	npoints_base
REAL(KIND=iii10), DIMENSION(:,:), ALLOCATABLE	::      raycoords_base
INTEGER 					:: 	crazyray_base, nru_base


END MODULE fm2dout



