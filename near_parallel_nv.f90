program near_parallel
   ! use cudafor
    use openacc
    use omp_lib

!!    This program was developed to explore the use of GPUs in a calculation that is similar to other
!!    that are part of the PRISM model.  It reads a prism grid and a prism station file, uses a lat/long box 
!!    to identify stations that are near each grid cell, calculates the actual distances between each gridcell
!!    and those stations, sorts the result, then adds the offsets to the closest stations to another array.
!!
!!    PRISM would then use that array of station offsets to identify, for each gridcell, the stations that would be
!!    used in each of the regressions (one regression for each grid cell)
!!    
!!    have used -stdpar=gpu, -stdpar=multicore and -mp=gpu, -mp=multicore.    
!!        -Minfo=accel produces a rundown of parallelization at compile    
!!
!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!         ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!                |              |  
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!                |              | 
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~

!!    haversine
!!    write_bil
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~


    implicit none

    integer, parameter :: num_stations = 26911
    integer, parameter :: ncols = 7025
    integer, parameter :: nrows = 3105
    integer, parameter :: NDISTANCES = 50  !The number of closest stations to identify.  This is the size of the 2 arrays that are placed in the register.
                                           !H100 has 255 32bit registers per thread. NDISTANCE > ~50 requires 100 registers. Anything above that degrades performance 

    integer, dimension(nrows, ncols) :: theElevGrid  ! Array to store the elevationraster data
    real(4), dimension(NDISTANCES) :: local_distance
    integer, dimension(NDISTANCES) :: local_station
    real(4), dimension(:,:), allocatable :: lats  ! Array to store latitude for the raster
    real(4), dimension(:,:), allocatable :: lons  ! Array to store longitude for the raster
    real(4), dimension(:), allocatable :: latitudes !latitude for the station file
    real(4), dimension(:), allocatable :: longitudes !longitude for the station file
    integer, dimension(:,:,:), allocatable :: near_stations  ! Array to store the local station list (Those within the threshold). The kth dim will be the fffset into the station arrays 
    real(4), dimension(:,:,:), allocatable :: distances  !Distance from cell to each of the stations within threshold 
    real(4), dimension(:,:,:), allocatable :: elev_weight  !Distance from cell to each of the stations within threshold 
    real(4), dimension(:,:,:), allocatable :: distance_weight  !Distance from cell to each of the stations within threshold 
    real(4), dimension(:,:,:), allocatable :: proximity_weight  !Coastal proximity weight 



    character(len=100) :: filename, filename2
    integer :: ios
    integer :: i, j, k
    real(4) :: ul_lat, ul_lon, cell_size, no_data
 
    !!!!! Arrays to store indivdual columns from the stations list
    character(len=30), dimension(num_stations) :: station_names
    integer, dimension(num_stations) :: station_ids, elevations, other_int
    character(len=30), dimension(num_stations) :: additional_info
    character(len=100) :: header_line
    !!!!!!!!!!
  
    integer :: start_clock, end_clock, clock_rate, indx
    real(4) :: es1, distance, lats_scalar, lons_scalar
    real(4)  :: distb, dz

    real(4), parameter :: rm = 10 !minimum distance, in km
    real(4), parameter :: a = 2    !distance weighting exponent.  2 means inverse distance squared
    real(4), parameter :: b = 1  ! elevation weighting exponent
    real(4), parameter ::dzm=200 ! minimum elevation difference meters
    real(4), parameter ::dzx=1000 ! max elevation difference meters
    
    !These should be read from the bil, but are consistent across all prism grids
    ul_lat = 49.933333333333    ! Upper-left corner latitude
    ul_lon = -125.016666666667    ! Upper-left corner longitude
    cell_size = 0.008333333333    ! Cell size in degrees

    distance=2 ! 'Distance' in degrees lon/lat that defines the coarse seach rectangle
    no_data=-9999
    
    allocate( distances(nrows,ncols,NDISTANCES),lats(nrows,ncols),lons(nrows,ncols),latitudes(num_stations),longitudes(num_stations),  near_stations(nrows,ncols,NDISTANCES))
    allocate( elev_weight(nrows,ncols,NDISTANCES),distance_weight(nrows,ncols,NDISTANCES))

    distances(:,:,:)=-9999
    local_distance(:)=-9999
    ! File containing the station data
    filename = 'cai_tmax_us_us_30s_20240318.stn'

    ! Open the text file for reading
    open(unit=10, file=filename, status='old', action='read',form='formatted' ,iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(filename)
        stop
    end if

    ! Skip the header lines
    do i = 1, 2
        read(10, '(A)', iostat=ios) header_line
        if (ios /= 0) then
            print *, 'Error reading header line: ', i
            stop
        end if
    end do

    ! Read the data into arrays
    i = 0
    do
        i = i + 1
        read(10, '(I8,  A13, I5, F13.4,  F12.4, I6,  A10)', iostat=ios) &
            station_ids(i), station_names(i), elevations(i), longitudes(i), latitudes(i), other_int(i), additional_info(i)
        if (ios < 0) exit
    end do
    close(10)

    filename2 = 'PRISM_us_dem_800m_bil.bil'
    ! Open the .bil file for reading
    open(unit=10, file=filename2, status='old', access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', trim(filename2)
        stop
    end if
    do i = 1, nrows
        read(10, rec=i) theElevGrid(i,:)
    end do
   ! read(10) theElevGrid
    close(10)

    ! Initialize to no_data so that the output bils display correctly in arc
    !$acc parallel loop gang vector collapse(2) 
    do j = 1, ncols  ! Across CUDA blocks, for example
        do i = 1, nrows  ! Threads in a block, for example
            !$acc loop seq
            do k = 1, NDISTANCES
                near_stations(i,j,k)=-9999
                elev_weight(i,j,k)=-9999
                distance_weight(i,j,k)=-9999
            end do
        end do    
    end do    

! Find the closest stations, building arrays to store distances and station offsets
    call system_clock(start_clock, clock_rate)
    ! Find and sort the closest NDISTANCES stations.
    !$acc parallel loop gang vector collapse(2) private(local_distance, local_station)
        do j = 1, ncols  ! Across CUDA blocks, for example
            do i = 1, nrows  ! Threads in a block, for example
                if (theElevGrid(i,j) .ne. -9999) then  ! A quick check??? - this check produces vertical lines of 0s.  Not sure why.  
                    lats_scalar = ul_lat - (i - 1) * cell_size  ! scalar
                    lons_scalar = ul_lon + (j - 1) * cell_size  ! scalar            
                    local_distance(:) = 250
                    local_station(:)  = 0
                    !$acc loop seq
                    do k = 1, num_stations
                        ! Am I within the minimal distance to even consider
                        if ((abs(lats_scalar-latitudes(k))  < distance) .and. (abs(lons_scalar-longitudes(k)) < distance)) then
                            ! Yes, so compute a better distance
                         !   distb = sqrt((lats_scalar-latitudes(k))**2 + (lons_scalar-longitudes(k))**2)
                            distb = haversine(lats_scalar, lons_scalar, latitudes(k), longitudes(k))
                            ! Now put it in the proper place of distances
                            if (distb .lt. local_distance(NDISTANCES)) then
                            ! At least one in the list will change
                                local_distance(NDISTANCES) = distb
                                local_station(NDISTANCES) = k
                                !$acc loop seq
                                do indx = NDISTANCES-1, 1, -1
                                    if (distb .ge. local_distance(indx)) exit
                                    ! This might not be optimal, even correct.  But move bigger things up.
                                    local_distance(indx+1) = local_distance(indx)
                                    local_station(indx+1)  = local_station(indx)
                                    local_distance(indx) = distb
                                    local_station(indx) = k
                                end do
                            end if
                        end if
                    end do
                    distances(i,j,:) = local_distance(:) ! Finally write them to global memory, just once
                    near_stations(i,j,:) = local_station(:)
               end if
            end do
        end do    


            ! Calculate the elevation and distance weight
        !$acc parallel loop gang vector collapse(2) 
        do j = 1, ncols  ! Across CUDA blocks, for example
            do i = 1, nrows  ! Threads in a block, for example
                if (theElevGrid(i,j) .ne. -9999) then 
                    !$acc loop seq
                    do k = 1, NDISTANCES
                        ! Elevation Weight
                        dz=abs(theElevGrid(i,j) - elevations(near_stations(i,j,k)))
                        if (dz .le. dzm) then
                            elev_weight(i,j,k)=1
                        else if (dz .gt. dzm .and. dz .le. dzx) then
                            elev_weight(i,j,k)=1/(dz**b)
                        else if (dz .gt. dzx) then
                            elev_weight(i,j,k)=0
                        end if
                        !Distance Weigth
                        if ((distances(i,j,k) - rm) .le. 0) then
                            distance_weight(i,j,k)=1
                        else
                            distance_weight(i,j,k)=1/((distances(i,j,k)-rm)**a)
                        end if
                        !Coastal Proximity Weight


                    end do
                end if
            end do    
        end do

    call system_clock(end_clock)

    es1 = real(end_clock - start_clock) / real(clock_rate)
    print '(A, F10.4, A)', 'Finding nearby stations AND complete distances AND weights with OpenACC: Elapsed time : ', es1, ' seconds. '
  
    call write_bil(distance_weight(:,:,1)) 
    deallocate( distances,lats,lons,latitudes,longitudes, near_stations,distance_weight, elev_weight)

    contains

            ! Pure function to calculate the Haversine distance between two points.  The accounts for the curvature of the earth and returns distance in m
        pure real(4) function haversine(lat1, lon1, lat2, lon2)
            real(4), intent(in) :: lat1, lon1, lat2, lon2
            real(4) :: dlat, dlon, a, c, r
            real(4), parameter :: pi = 3.141592653589793
            real(4) :: lat1_rad, lat2_rad
            r = 6371.0  ! Radius of the Earth in kilometers

            lat1_rad = lat1 * pi / 180.0
            lat2_rad = lat2 * pi / 180.0
            dlat = (lat2 - lat1) * pi / 180.0
            dlon = (lon2 - lon1) * pi / 180.0

            a = sin(dlat / 2.0)**2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon / 2.0)**2
            c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
            haversine=r * c          
        end function haversine  



!  Write the data to a (currently hardcoded!).bil format. Note that this includes only the data.  Additional files needed to describe data are not written.

        subroutine write_bil(array)
            integer :: nx, ny, ii
            real(4), dimension(:,:), intent(in) :: array
            !integer, dimension(:,:), intent(inout) :: array
            ny = SIZE(array, 1)
            nx = SIZE(array, 2)
            open(unit=10, file='near_d.bil', status='replace', access='stream', form='unformatted')
        
            ! Write the array to the file in BIL format
            do ii = 1, ny
                write(10, rec=ii) array(ii,:)
            end do
          !  write(10) array
            close(10)
        end subroutine write_bil

end program near_parallel
