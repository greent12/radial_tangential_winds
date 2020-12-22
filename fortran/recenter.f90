PROGRAM recenter

USE read_test_data  !Has subroutine that reads in test data
USE kinds
USE constants
USE ideal_angle
USE calc_distance_from_point
USE uv_to_rt
USE maskmod
USE omp_lib

IMPLICIT NONE
INTEGER :: ios,n,ri,thi,i
INTEGER,PARAMETER :: GRIDDIMY=150,GRIDDIMX=150
REAL(KIND=single), DIMENSION(GRIDDIMY,GRIDDIMX) :: lons,lats,uwind,vwind
LOGICAL :: converged = .FALSE.

!Variables that come from the namelist
INTEGER :: num_sectors,num_iterations,spad,grf
REAL(KIND=single) :: dist_coeff,wind_coeff,rxpad,olon,olat,curr_delta,gr0

INTEGER :: pnxi,pnyi,xloc,yloc,ybi,xbi
INTEGER, DIMENSION(1) :: rmw_indx
REAL(KIND=single),ALLOCATABLE :: angle_thresh(:),rad_gaussian(:),&
                                 vt_ann_azi(:),sector_mean_error(:,:,:)
REAL(KIND=single),DIMENSION(GRIDDIMY,GRIDDIMX) :: ws,XD,YD, &
                  angle, curr_dist, curr_rw, curr_vt, dist_weight, &
                  curr_weight,curr_angle_dif, curr_weighted_dif
REAL(KIND=single) :: prev_mean_dif=FILL_VAL,tc_center_lon, &
                     tc_center_lat,vt_azi_max,curr_rmw,& 
                     ws_max,angle_upper,angle_lower,curr_mean_dif
REAL(KIND=single),DIMENSION(GRIDDIMY,GRIDDIMX,GRIDDIMY,GRIDDIMX) :: the_ideal_angle
LOGICAL,DIMENSION(GRIDDIMY,GRIDDIMX) :: mask_ann,mask_sector

!Read namelist variables
NAMELIST/parm/num_sectors,dist_coeff,wind_coeff,rxpad,olon,olat,&
              num_iterations,spad,curr_delta,gr0,grf
OPEN(UNIT=11,FILE="parm.nml",IOSTAT=ios)
READ(11,parm,IOSTAT=ios)
CLOSE(11)

!Print namelist variables to make sure they are correct
PRINT*, "***********PARAMETERS************"
PRINT*, "num_sectors",num_sectors
PRINT*, "dist_coeff",dist_coeff
PRINT*, "wind_coeff",wind_coeff
PRINT*, "rxpad",rxpad
PRINT*, "olon",olon
PRINT*, "olat",olat
PRINT*, "num_iterations",num_iterations
PRINT*, "spad",spad
PRINT*, "curr_delta",curr_delta
PRINT*, "gr0",gr0
PRINT*, "grf",grf
PRINT*, "*********************************"

!Initialize and allocate arrays
ws = FILL_VAL
XD = 0.0
YD = 0.0

ALLOCATE(angle_thresh(num_sectors+1),rad_gaussian(grf),&
         sector_mean_error(num_sectors,GRIDDIMY,GRIDDIMX),&
         vt_ann_azi(grf))

!Read in a simple case 150 x 150 grids from Hurricane Matthew 19z analysis HWRF run
CALL load_data(lats,lons,uwind,vwind)

!Calculate ideal angle array
CALL calc_ideal_angle(GRIDDIMY,GRIDDIMX,lats,lons,the_ideal_angle)

!Bounds of azimuthal sectors to average angle errors
angle_thresh = (/(-pi+((2*pi)/num_sectors)*i,i=0,num_sectors)/)

!Array of radii to use for radius of maximum wind computations
rad_gaussian = (/(gr0+i*curr_delta,i=0,grf-1)/)

!Find indicies of closest lat/lon to first guesses
pnyi=MINLOC(ABS(lats(:,1)-olat),1)
pnxi=MINLOC(ABS(lons(1,:)-olon),1)

!Compute wind speed
ws = SQRT(uwind**2+vwind**2)
ws_max = MAXVAL(ws,ws .NE. FILL_VAL)

!Begin TC center search
OUTER: DO n =1,num_iterations  !Replace with num_iterations

   PRINT*,'  recenter_tc: current iteration is:',n
   
   !If first iteration has been completed, copy center estimate from previous iteration:
   IF ( n .GT. 1 ) THEN
      IF ( prev_mean_dif .NE. FILL_VAL ) THEN
         pnyi = yloc
         pnxi = xloc
      ELSE
         tc_center_lon = FILL_VAL
         tc_center_lat = FILL_VAL
         EXIT OUTER 
      ENDIF
   ENDIF
 
   !Loop over meridional grid points
   YPOINT: DO ybi=pnyi-spad,pnyi+spad
      !Loop over zonal grid points
      XPOINT: DO xbi=pnxi-spad,pnxi+spad

         !Check to make sure were searching within the grid bounds
         IF (ybi .LT. 1 .OR. &
             ybi .GT. GRIDDIMY .OR. &
             xbi .LT. 1 .OR. &
             xbi .GT. GRIDDIMX) CYCLE XPOINT
          
         !Compute zonal displacement and relative angles
         XD = lons - lons(ybi,xbi)
         YD = lats - lats(ybi,xbi)
         angle = ATAN2(YD,XD)
     
         !Calculate distances and then tangential wind
         CALL calc_distance_g2p(lons(ybi,xbi),lats(ybi,xbi),&
                                GRIDDIMX,GRIDDIMY,&
                                lons,lats,curr_dist)
         CALL calc_uv_to_rt(GRIDDIMX,GRIDDIMY,lons,lats,lons(ybi,xbi),&
                            lats(ybi,xbi),uwind,vwind,curr_rw,curr_vt)

         !Identify the annulus where the tangential wind is maximized
         vt_ann_azi=FILL_VAL
         DO ri =1,grf
           mask_ann = (curr_dist .GE. rad_gaussian(ri)-0.5*curr_delta &
                   .AND. curr_dist < rad_gaussian(ri) + 0.5*curr_delta)
           CALL masked_mean_2d(curr_vt,mask_ann,vt_ann_azi(ri))
         ENDDO
         vt_azi_max = MAXVAL(vt_ann_azi,vt_ann_azi .NE. FILL_VAL)
         rmw_indx=MAXLOC(vt_ann_azi,vt_ann_azi .NE. FILL_VAL)
         curr_rmw = rad_gaussian(rmw_indx(1))

         !Establish the radial distance weighting for errors:
         dist_weight = EXP(-1.*(((curr_dist - 0.)**2)/(2.*(curr_rmw**2))))

         !Compute sum of weights 
         curr_weight = dist_coeff*dist_weight + wind_coeff*(ws/ws_max)

         !Compute the angle difference between idealized vortex and observed flow
         !This Where/Elsewhere block is acting as a array-array that might be
         ! masked or don's support nan math
         WHERE (the_ideal_angle(ybi,xbi,:,:) .NE. FILL_VAL)
          curr_angle_dif = ATAN2(vwind,uwind) - the_ideal_angle(ybi,xbi,:,:)
         ELSEWHERE
          curr_angle_dif = FILL_VAL
         ENDWHERE

         !Correct angle difference to account for values outside -PI<=x<=PI
         WHERE ( curr_angle_dif .EQ. FILL_VAL )
           curr_angle_dif = FILL_VAL
         ELSEWHERE ( curr_angle_dif .GT. PI )
           curr_angle_dif = curr_angle_dif - 2.*PI
         ELSEWHERE ( curr_angle_dif .LT. -1.*PI )
           curr_angle_dif = curr_angle_dif + 2.*PI
         ENDWHERE

         !Compute finalized errors at each grid point
         WHERE ( curr_angle_dif .EQ. FILL_VAL )
           curr_weighted_dif = FILL_VAL
         ELSEWHERE
           curr_weighted_dif = curr_weight*curr_angle_dif
         ENDWHERE

         !Set al values outside of RMW+rxpad to FILL_VAL
         WHERE ( curr_dist .GT. curr_rmw + rxpad )
           curr_weighted_dif = FILL_VAL
         ENDWHERE

         !Seperate grid into azimuthal sectors
         SECTORS: DO thi =1,num_sectors
           angle_lower=angle_thresh(thi)
           angle_upper=angle_thresh(thi+1)
           mask_sector = (angle .GE. angle_lower .AND. &
                          angle .LE. angle_upper .AND. &
                          curr_weighted_dif .NE. FILL_VAL )

           CALL masked_mean_2d(ABS(curr_weighted_dif),&
                               mask_sector,&
                               sector_mean_error(thi,ybi,xbi))
         ENDDO SECTORS

         !Compute mean error of quadrants  
         CALL masked_mean_1d(sector_mean_error(:,ybi,xbi),&
                             sector_mean_error(:,ybi,xbi) &
                             .NE. fill_val,curr_mean_dif)

         !See if current location yields a better center estimate
         ! than the previous iteration
         IF (curr_mean_dif .LT. prev_mean_dif) THEN
           !Important that FILL_VAL is set to positive (large number)
           ! so that on first iteration prev_mean_dif will be set to 
           ! curr_mean_dif
           prev_mean_dif = curr_mean_dif
           tc_center_lon = lons(ybi,xbi)
           tc_center_lat = lats(ybi,xbi)
           yloc = ybi
           xloc = xbi
         ENDIF
         
      ENDDO XPOINT
   ENDDO YPOINT

   !If the estimate has converged, break from OUTER loop and keep estimate
   IF (n .GT. 1) THEN
     IF (yloc .EQ. pnyi .AND. xloc .EQ. pnxi) THEN
        converged = .true.
        EXIT OUTER 
     ENDIF
   ENDIF
ENDDO OUTER


print*, converged
print*, tc_center_lat,tc_center_lon,xloc,yloc

DEALLOCATE(angle_thresh,rad_gaussian,&
         sector_mean_error,vt_ann_azi)

END PROGRAM
