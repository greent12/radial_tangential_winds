MODULE maskmod

USE kinds
USE constants

IMPLICIT NONE

CONTAINS

SUBROUTINE masked_mean_2d(array,mask,mean)
!Calculates the mean of 'array' for all entries where 'mask' is true
!Array is assumed to be of the type real and be a 'single' as defined in kinds.f90

REAL(KIND=single),DIMENSION(:,:),INTENT(IN) :: array
LOGICAL,DIMENSION(:,:), INTENT(IN) :: mask
REAL(KIND=single),INTENT(OUT) :: mean


IF (COUNT(mask) > 0) THEN
  mean = SUM(array,mask)/COUNT(mask)
ELSE
  mean = FILL_VAL
ENDIF

END SUBROUTINE

SUBROUTINE masked_mean_1d(array,mask,mean)
!Calculates the mean of 'array' for all entries where 'mask' is true
!Array is assumed to be of the type real and be a 'single' as defined in kinds.f90

REAL(KIND=single),DIMENSION(:),INTENT(IN) :: array
LOGICAL,DIMENSION(:), INTENT(IN) :: mask
REAL(KIND=single),INTENT(OUT) :: mean


IF (COUNT(mask) > 0) THEN
  mean = SUM(array,mask)/COUNT(mask)
ELSE
  mean = FILL_VAL
ENDIF

END SUBROUTINE
END MODULE
