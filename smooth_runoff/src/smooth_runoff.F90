program smooth_runoff

  use netcdf

  implicit none

  integer, parameter :: &
       RKIND = selected_real_kind(13)

  ! input parameters
  character(len=1000) :: &
       filename_mapping_in, &
       filename_mapping_out, &
       filename_mpas_mesh, &
       filename_diagnostic

  real(kind=RKIND) :: &
       distanceLimit, & ! (km)
       stdDeviation     ! (km)

!      if (iCell_in == 95036) then  ! Amazon
!      if (iCell_in == 29684) then  ! Yukon
!      if (iCell_in == 29883) then  ! Macenzie
  integer, parameter :: &
       singleOceanCell = 0

  integer :: &
       n_s_buffer, rec_length

  logical :: &
       outputDiagnostic

  integer :: &
       diagnosticPointIndex

  namelist /smooth/ &
       filename_mapping_in, &
       filename_mapping_out, &
       filename_mpas_mesh, &
       filename_diagnostic, &
       distanceLimit, &
       stdDeviation, &
       n_s_buffer, &
       outputDiagnostic, &
       diagnosticPointIndex

  ! input mapping file
  ! dimensions
  integer :: &
       n_a, &
       ni_a, &
       nj_a, &
       nv_a, &
       src_grid_rank, &
       n_b, &
       ni_b, &
       nj_b, &
       nv_b, &
       dst_grid_rank, &
       n_s

  ! variables
  real(kind=RKIND), dimension(:), allocatable :: &
       xc_a, &
       yc_a, &
       area_a, &
       frac_a, &
       xc_b, &
       yc_b, &
       area_b, &
       frac_b, &
       S

  real(kind=RKIND), dimension(:,:), allocatable :: &
       xv_a, &
       yv_a, &
       xv_b, &  
       yv_b

  integer, dimension(:), allocatable :: &
       mask_a, &
       src_grid_dims, &
       mask_b, &
       dst_grid_dims, &
       col, & ! RTM combined index
       row    ! MPAS iCell

  ! output mapping file
  integer*8 :: &
       n_s_out

  real(kind=RKIND), dimension(:), allocatable :: &
       S_out

  integer, dimension(:), allocatable :: &
       col_out, &
       row_out

  ! mpas grid
  integer :: &
       nCells, &
       nEdges, &
       maxEdges

  real(kind=RKIND), dimension(:), allocatable :: &
       xCell, &
       yCell, &
       zCell, &
       areaCell, &
       dcEdge

  integer, dimension(:), allocatable :: &
       nEdgesOnCell

  integer, dimension(:,:), allocatable :: &
       cellsOnCell, &
       edgesOnCell

  ! temporary variables
  integer :: &
       is, &
       ia, &
       ib, &
       iBuffer, &
       nBuffers, &
       start, &
       count, &
       ncid, &
       varid_S, &
       varid_col, &
       varid_row
  
  real(kind=RKIND), dimension(:), allocatable :: &
       weights, &
       cellDistance, &
       distanceSpread

  integer, dimension(:), allocatable :: &
       iCellsSpread, &
       cellStatus

  ! diagnostic variables
  real(kind=RKIND), dimension(:), allocatable :: &
       weightsDiagnostic

  real(kind=RKIND) :: &
       sumWeightsAreaOceanOrig, sumWeightsAreaOceanNew, &
       sumWeightsAreaRunoffOrig, sumWeightsAreaRunoffNew, &
       areaRunoffCell

  ! read in namelist
  write(*,*) "Reading namelist..."
  open(11,file="smooth_runoff_in",status='old')
  read(11,nml=smooth)
  close(11)

  ! load the nearest neighbour mapping file
  call load_mapping_file(filename_mapping_in)

  ! load the MPAS mesh file
  call load_mesh_file()

  ! allocate temporary arrays
  write(*,*) "Allocate temporary variables..."
  allocate(weights(nCells))
  allocate(cellDistance(nCells))
  allocate(distanceSpread(nCells))
  allocate(iCellsSpread(nCells))
  allocate(cellStatus(nCells))

  write(*,*) "Perform smoothing..."

  ! open a temporary file to store the output sparse matrix
!maltrud debug
! open(11,status="scratch")
! open(11,name="scratch",status="unknown")
  is = 0
  inquire(iolength = rec_length) is
  open(11,file="intScratch",access='direct',form='unformatted',  &
                recl=2*rec_length,status='unknown')

  inquire(iolength = rec_length) weights(1)
  open(12,file="r8Scratch",access='direct',form='unformatted',  &
                recl=rec_length,status='unknown')

  !  initial check
  sumWeightsAreaOceanOrig = 0.
  sumWeightsAreaRunoffOrig = 0. 
  
  !$OMP parallel do private(areaRunoffCell) reduction(+:sumWeightsAreaOceanOrig,sumWeightsAreaRunoffOrig)
  do is = 1, n_s
!    if (mod(is,1000) == 0) write(*,*) "INITIAL CHECK: Doing ", is, "of", n_s, " (", nint(100.0_RKIND*(real(is,RKIND)/real(n_s,RKIND))) ,"%)..."
     if (singleOceanCell /= 0) then
      if (row(is) == singleOceanCell) then
       sumWeightsAreaOceanOrig = sumWeightsAreaOceanOrig + S(is)*areaCell(row(is))
       areaRunoffCell = S(is)*areaCell(row(is))
       sumWeightsAreaRunoffOrig = sumWeightsAreaRunoffOrig + areaRunoffCell
      endif
     else
      sumWeightsAreaOceanOrig = sumWeightsAreaOceanOrig + S(is)*areaCell(row(is))
      areaRunoffCell = S(is)*areaCell(row(is))
      sumWeightsAreaRunoffOrig = sumWeightsAreaRunoffOrig + areaRunoffCell
     endif
!maltrud debug
! if (col(is) == 1000) write(*,*)col(is),row(is),areaRunoffCell,S(is),sumWeightsAreaRunoffOrig
  enddo ! is
  !$OMP end parallel do

  ! perform the smoothing
  n_s_out = 0
  do is = 1, n_s
     if (mod(is,1000) == 0) write(*,*) "Doing ", is, "of", n_s, " (", nint(100.0_RKIND*(real(is,RKIND)/real(n_s,RKIND))) ,"%)..."
     call perform_smoothing(col(is), row(is), S(is), 11, distanceLimit, stdDeviation)
  enddo ! is

! rewind(11)
  close(11)
  close(12)

  inquire(iolength = rec_length) is
  open(11,file="intScratch",access='direct',form='unformatted',  &
                recl=2*rec_length,status='old')

  inquire(iolength = rec_length) weights(1)
  open(12,file="r8Scratch",access='direct',form='unformatted',  &
                recl=rec_length,status='old')

  ! output the mapping file
  call create_mapping_file(filename_mapping_out, ncid, varid_S, varid_col, varid_row)

  write(*,*) "Read temporary file in..."

  ! allocate the output sparse matrix arrays
  allocate(S_out(n_s_buffer))
  allocate(col_out(n_s_buffer))
  allocate(row_out(n_s_buffer))

  ! diagnostics
  if (outputDiagnostic) then
     allocate(weightsDiagnostic(nCells))
     weightsDiagnostic = 0.0_RKIND
  endif

  !  final check
  sumWeightsAreaOceanNew = 0.
  sumWeightsAreaRunoffNew = 0. 

  nBuffers = ceiling(real(n_s_out, RKIND) / real(n_s_buffer, RKIND))

  do iBuffer = 1, nBuffers

     start = n_s_buffer*(iBuffer-1) + 1
     count = min(n_s_buffer,n_s_out-start+1)

     write(*,*) "Buffer start: ", start, ", end: ", start+count-1, ", total: ", n_s_out, " (", nint(100.0_RKIND*(real(start,RKIND)/real(n_s_out,RKIND))) ,"%)..."

     ! read from the temporary file the sparse matrix 
     do is = 1, count
        read(11,rec=start+is-1) col_out(is), row_out(is)
        read(12,rec=start+is-1) S_out(is)
!maltrud debug
!       if(S_out(is) /= 0.0) write(*,*)'link,col,row,S = ',is,col_out(is),row_out(is),S_out(is)

!maltrud debug
! if(row_out(is) == 29684) write(*,*)'link,col,row,S = ',is,col_out(is),row_out(is),S_out(is)

        sumWeightsAreaOceanNew = sumWeightsAreaOceanNew + S_out(is)*areaCell(row_out(is))
        do ib = 1, n_s
          if (col_out(is) == col(ib)) then
            areaRunoffCell = S(ib)*areaCell(row(ib))
            sumWeightsAreaRunoffNew = sumWeightsAreaRunoffNew + areaRunoffCell
!maltrud debug
! if (col(ib) == 1000) write(*,*)ib,is,col(ib),row(ib),areaRunoffCell,S(ib)
          endif
        enddo
     enddo

     if (outputDiagnostic) then
        do is = 1, n_s_buffer
           if (col_out(is) == diagnosticPointIndex) then
              weightsDiagnostic(row_out(is)) = S_out(is)
           endif
        enddo ! is
     endif

     ! write out 
     call write_mapping_file(ncid, varid_S, varid_col, varid_row, start, count, S_out, col_out, row_out)

  enddo ! iBuffer

  ! output the mapping file
  call close_mapping_file(ncid)

  write(*,*)'sumWeightsAreaOceanOrig, sumWeightsAreaOceanNew = ',  &
       sumWeightsAreaOceanOrig, sumWeightsAreaOceanNew
  write(*,*)'sumWeightsAreaRunoffOrig = ',  &
       sumWeightsAreaRunoffOrig

  ! output the diagnostics file
  if (outputDiagnostic) &
       call output_weights_diagnostic(filename_diagnostic, weightsDiagnostic)
  
contains

  !----------------------------------------------------------------

  subroutine perform_smoothing(col, iCell_in, S_col, unitOut, distanceLimit, stdDeviation)

    integer, intent(in) :: &
         col, iCell_in, unitOut

    real(kind=RKIND), intent(in) :: &
         S_col,  &
         distanceLimit, &
         stdDeviation

    integer :: &
         nSpread, &
         iCell

    real(kind=RKIND) :: &
         distance, &
         sum_weights, &
         areaRunoffCell

    real(kind=RKIND), parameter :: &
         km_to_m = 1000.0_RKIND

    call get_active_distances(&
         iCell_in, &
         distanceLimit * km_to_m, &
         nEdgesOnCell, &
         cellsOnCell, &
         edgesOnCell, &
         dcEdge, &
         cellStatus, &
         cellDistance, &
         nSpread, &
         iCellsSpread, &
         distanceSpread)

!maltrud debug
!   if (nSpread > 1) write(*,*)' nSpread > 1!!! ', iCell_in
    do iCell = 1, nSpread
       weights(iCell) = gaussian_weight(distanceSpread(iCell), stdDeviation * km_to_m)
    enddo ! iCell
       
    ! renormalize weights
    sum_weights = sum(weights(1:nSpread))
    weights(1:nSpread) = weights(1:nSpread)/sum_weights

    areaRunoffCell = S_col*areaCell(iCell_in)
    do iCell = 1, nSpread
       weights(iCell) = weights(iCell) * areaRunoffCell/areaCell(iCellsSpread(iCell))
    enddo ! iCell

    ! write weights to temporary file
    do iCell = 1, nSpread

    !print*,'n_s_out ',icell,n_s_out
       write(unitOut,rec=n_s_out+1) col, iCellsSpread(iCell)
       if (iCell_in == singleOceanCell .or. singleOceanCell == 0) then
!        write(*,*)'col,row,iCellspread,weight = ',col,iCell_in,iCellsSpread(iCell),weights(iCell)
         write(unitOut+1,rec=n_s_out+1) weights(iCell)
       else
         write(unitOut+1,rec=n_s_out+1) weights(iCell) * 0.0
       endif
       n_s_out = n_s_out + 1

    enddo ! iCell

  end subroutine perform_smoothing

  !----------------------------------------------------------------

  function gaussian_weight(distance, stdDeviation) result(weight)

    real(kind=RKIND), intent(in) :: &
         distance, &
         stdDeviation

    real(kind=RKIND) :: weight

    weight = exp(-1.0_RKIND * (distance**2 / (2.0_RKIND * stdDeviation**2)))
    
  end function gaussian_weight

  !----------------------------------------------------------------

  subroutine get_active_distances(&
       iCell_Start, &
       distanceLimit, &
       nEdgesOnCell, &
       cellsOnCell, &
       edgesOnCell, &
       dcEdge, &
       cellStatus, &
       cellDistance, &
       nCellsOut, &
       iCellsOut, &
       distanceOut)

    integer, intent(in) :: &
         iCell_Start

    real(kind=RKIND), intent(in) :: &
         distanceLimit

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell

    integer, dimension(:,:), intent(in) :: &
         cellsOnCell, &
         edgesOnCell

    real(kind=RKIND), dimension(:), intent(in) :: &
         dcEdge

    integer, dimension(:), intent(out) :: &
         cellStatus, &
         iCellsOut

    real(kind=RKIND), dimension(:), intent(out) :: &
         cellDistance, &
         distanceOut

    integer, intent(out) :: &
         nCellsOut

    integer, dimension(:), allocatable :: &
         iCellsPrev, &
         iCellsNext

    integer :: &
         nCellsPrev, &
         nCellsNext, &
         nLayer, &
         iCellPrev, &
         iCellOnCell, &
         iCellNext, &
         iEdge

    real(kind=RKIND) :: &
         cellDistanceToAdd      
    
    allocate(iCellsPrev(nCells))
    allocate(iCellsNext(nCells))

    ! set the internal variables
    cellStatus   = -1
    cellDistance = -1.0_RKIND

    ! set first cell 
    cellStatus(iCell_Start) = 0
    nCellsPrev = 1
    iCellsPrev(1) = iCell_Start

    ! first output cell is the starting cell
    nCellsOut = 1
    iCellsOut(1) = iCell_Start
    distanceOut(1) = 0.0_RKIND

    ! loop over cell layers from original cell
    nLayer = 0
    do while (nCellsPrev > 0)

       ! reset number of cells found this layer
       nCellsNext = 0

       ! loop over cells defined in the previous iteration
       do iCellPrev = 1, nCellsPrev

          ! loop over neighbours of these previous cells
          do iCellOnCell = 1, nEdgesOnCell(iCellsPrev(iCellPrev))

             ! get the iCell of the next potential cell in the next cells
             iCellNext = cellsOnCell(iCellOnCell,iCellsPrev(iCellPrev))

             ! get the edge index for this crossing
             iEdge = edgesOnCell(iCellOnCell,iCellsPrev(iCellPrev))

             cellDistanceToAdd = cellDistance(iCellsPrev(iCellPrev)) + dcEdge(iEdge)

             ! check to see if we need to add it to the next array
             if (icellnext .ne. 0) then
             if (cellStatus(iCellNext) == -1 .and. cellDistanceToAdd < distanceLimit) then

                ! count how many on the next list
                nCellsNext = nCellsNext + 1

                ! add this new cell to the next list
                iCellsNext(nCellsNext) = iCellNext

                ! update the status of the cell
                cellStatus(iCellNext) = nLayer

                ! calculate the distance to this cell
                cellDistance(iCellNext) = cellDistanceToAdd

                ! output
                nCellsOut = nCellsOut + 1
                iCellsOut(nCellsOut) = iCellNext
                distanceOut(nCellsOut) = cellDistanceToAdd

             endif ! cellStatus(iCellNext) == -1
             endif

          enddo ! iCellOnCell

       enddo ! iCellPrev

       ! swap next and prev
       nCellsPrev = nCellsNext

       iCellsPrev(1:nCellsNext) = iCellsnext(1:nCellsNext)

       ! increment the layer number
       nLayer = nLayer + 1

    enddo ! nCellsNext > 0

    deallocate(iCellsPrev)
    deallocate(iCellsNext)

  end subroutine get_active_distances

  !----------------------------------------------------------------

  subroutine create_mapping_file(filename, ncid, varid_S, varid_col, varid_row)

    character(len=*), intent(in) :: &
         filename

    integer, intent(out) :: &
         ncid, &
         varid_S, &
         varid_col, &
         varid_row
    
    integer :: &
         status

    integer :: &
         dimid_n_a, &
         dimid_ni_a, &
         dimid_nj_a, &
         dimid_nv_a, &
         dimid_src_grid_rank, &
         dimid_n_b, &
         dimid_ni_b, &
         dimid_nj_b, &
         dimid_nv_b, &
         dimid_dst_grid_rank, &
         dimid_n_s

    integer :: &
         varid_xc_a, &
         varid_yc_a, &
         varid_xv_a, &
         varid_yv_a, &
         varid_mask_a, &
         varid_area_a, &
         varid_frac_a, &
         varid_src_grid_dims, &
         varid_xc_b, &
         varid_yc_b, &
         varid_xv_b, &  
         varid_yv_b, &
         varid_mask_b, &
         varid_area_b, &
         varid_frac_b, &
         varid_dst_grid_dims

    character(len=1000) :: &
         date, &
         time, &
         zone, &
         datetime

    write(*,*) "Create mapping file...", n_s_out

    ! create
!   status = nf90_create(trim(filename), NF90_CLOBBER, ncid)
    status = nf90_create(trim(filename), NF90_64BIT_OFFSET, ncid)
    call netcdf_error(status, "create_mapping_file: nf90_open")

    ! define dimensions
    ! n_a
    status = nf90_def_dim(ncid, "n_a", n_a, dimid_n_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_dim n_a")

    !! ni_a
    !status = nf90_def_dim(ncid, "ni_a", ni_a, dimid_ni_a)
    !call netcdf_error(status, "create_mapping_file: nf90_def_dim ni_a")

    !! nj_a
    !status = nf90_def_dim(ncid, "nj_a", nj_a, dimid_nj_a)
    !call netcdf_error(status, "create_mapping_file: nf90_def_dim nj_a")

    ! nv_a
    status = nf90_def_dim(ncid, "nv_a", nv_a, dimid_nv_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_dim nv_a")

    ! src_grid_rank
    status = nf90_def_dim(ncid, "src_grid_rank", src_grid_rank, dimid_src_grid_rank)
    call netcdf_error(status, "create_mapping_file: nf90_def_dim src_grid_rank")

    ! n_b
    status = nf90_def_dim(ncid, "n_b", n_b, dimid_n_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_dim n_b")

    !! ni_b
    !status = nf90_def_dim(ncid, "ni_b", ni_b, dimid_ni_b)
    !call netcdf_error(status, "create_mapping_file: nf90_def_dim ni_b")

    !! nj_b
    !status = nf90_def_dim(ncid, "nj_b", nj_b, dimid_nj_b)
    !call netcdf_error(status, "create_mapping_file: nf90_def_dim nj_b")

    ! nv_b
    status = nf90_def_dim(ncid, "nv_b", nv_b, dimid_nv_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_dim nv_b")

    ! dst_grid_rank
    status = nf90_def_dim(ncid, "dst_grid_rank", dst_grid_rank, dimid_dst_grid_rank)
    call netcdf_error(status, "create_mapping_file: nf90_def_dim dst_grid_rank")
    
    ! n_s
    status = nf90_def_dim(ncid, "n_s", int(n_s_out,4), dimid_n_s)
    call netcdf_error(status, "create_mapping_file: nf90_def_dim n_s")

    ! define variables
    ! xc_a
    status = nf90_def_var(ncid, "xc_a", NF90_DOUBLE, (/dimid_n_a/), varid_xc_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_var xc_a")

    status = nf90_put_att(ncid, varid_xc_a, "long_name", "longitude of grid cell center (input)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xc_a")

    status = nf90_put_att(ncid, varid_xc_a, "units", "degrees east")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xc_a")

    ! yc_a
    status = nf90_def_var(ncid, "yc_a", NF90_DOUBLE, (/dimid_n_a/), varid_yc_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_var yc_a")

    status = nf90_put_att(ncid, varid_yc_a, "long_name", "latitude of grid cell center (input)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yc_a")

    status = nf90_put_att(ncid, varid_yc_a, "units", "degrees north")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yc_a")
 
    ! xv_a
    status = nf90_def_var(ncid, "xv_a", NF90_DOUBLE, (/dimid_nv_a,dimid_n_a/), varid_xv_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_var xv_a")

    status = nf90_put_att(ncid, varid_xv_a, "long_name", "longitude of grid cell verticies (input)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xv_a")

    status = nf90_put_att(ncid, varid_xv_a, "units", "degrees east")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xv_a")

    ! yv_a
    status = nf90_def_var(ncid, "yv_a", NF90_DOUBLE, (/dimid_nv_a,dimid_n_a/), varid_yv_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_var yv_a")

    status = nf90_put_att(ncid, varid_yv_a, "long_name", "latitude of grid cell verticies (input)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yv_a")

    status = nf90_put_att(ncid, varid_yv_a, "units", "degrees north")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yv_a")

    ! mask_a
    status = nf90_def_var(ncid, "mask_a", NF90_INT, (/dimid_n_a/), varid_mask_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_var mask_a")

    status = nf90_put_att(ncid, varid_mask_a, "long_name", "domain mask (input)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att mask_a")

    ! area_a
    status = nf90_def_var(ncid, "area_a", NF90_DOUBLE, (/dimid_n_a/), varid_area_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_var area_a")

    status = nf90_put_att(ncid, varid_area_a, "long_name", "area of cell (input)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att area_a")

    ! frac_a
    status = nf90_def_var(ncid, "frac_a", NF90_DOUBLE, (/dimid_n_a/), varid_frac_a)
    call netcdf_error(status, "create_mapping_file: nf90_def_var frac_a")

    status = nf90_put_att(ncid, varid_frac_a, "long_name", "fraction of domain intersection (input)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att frac_a")

    ! src_grid_dims
    status = nf90_def_var(ncid, "src_grid_dims", NF90_INT, (/dimid_src_grid_rank/), varid_src_grid_dims)
    call netcdf_error(status, "create_mapping_file: nf90_def_var src_grid_dims")

    ! xc_b
    status = nf90_def_var(ncid, "xc_b", NF90_DOUBLE, (/dimid_n_b/), varid_xc_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_var xc_b")

    status = nf90_put_att(ncid, varid_xc_b, "long_name", "longitude of grid cell center (output)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xc_b")

    status = nf90_put_att(ncid, varid_xc_b, "units", "degrees east")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xc_b")

    ! yc_b
    status = nf90_def_var(ncid, "yc_b", NF90_DOUBLE, (/dimid_n_b/), varid_yc_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_var yc_b")

    status = nf90_put_att(ncid, varid_yc_b, "long_name", "latitude of grid cell center (output)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yc_b")

    status = nf90_put_att(ncid, varid_yc_b, "units", "degrees north")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yc_b")
 
    ! xv_b
    status = nf90_def_var(ncid, "xv_b", NF90_DOUBLE, (/dimid_nv_b,dimid_n_b/), varid_xv_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_var xv_b")

    status = nf90_put_att(ncid, varid_xv_b, "long_name", "longitude of grid cell verticies (output)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xv_b")

    status = nf90_put_att(ncid, varid_xv_b, "units", "degrees east")
    call netcdf_error(status, "create_mapping_file: nf90_put_att xv_b")

    ! yv_b
    status = nf90_def_var(ncid, "yv_b", NF90_DOUBLE, (/dimid_nv_b,dimid_n_b/), varid_yv_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_var yv_b")

    status = nf90_put_att(ncid, varid_yv_b, "long_name", "latitude of grid cell verticies (output)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yv_b")

    status = nf90_put_att(ncid, varid_yv_b, "units", "degrees north")
    call netcdf_error(status, "create_mapping_file: nf90_put_att yv_b")

    ! mask_b
    status = nf90_def_var(ncid, "mask_b", NF90_INT, (/dimid_n_b/), varid_mask_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_var mask_b")

    status = nf90_put_att(ncid, varid_mask_b, "long_name", "domain mask (output)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att mask_b")

    ! area_b
    status = nf90_def_var(ncid, "area_b", NF90_DOUBLE, (/dimid_n_b/), varid_area_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_var area_b")

    status = nf90_put_att(ncid, varid_area_b, "long_name", "area of cell (output)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att area_b")

    ! frac_b
    status = nf90_def_var(ncid, "frac_b", NF90_DOUBLE, (/dimid_n_b/), varid_frac_b)
    call netcdf_error(status, "create_mapping_file: nf90_def_var frac_b")

    status = nf90_put_att(ncid, varid_frac_b, "long_name", "fraction of domain intersection (output)")
    call netcdf_error(status, "create_mapping_file: nf90_put_att frac_b")

    ! dst_grid_dims
    status = nf90_def_var(ncid, "dst_grid_dims", NF90_INT, (/dimid_dst_grid_rank/), varid_dst_grid_dims)
    call netcdf_error(status, "create_mapping_file: nf90_def_var dst_grid_dims")

    ! S
    status = nf90_def_var(ncid, "S", NF90_DOUBLE, (/dimid_n_s/), varid_S)
    call netcdf_error(status, "create_mapping_file: nf90_def_var S")

    status = nf90_put_att(ncid, varid_S, "long_name", "sparse matrix for mapping S:a->b")
    call netcdf_error(status, "create_mapping_file: nf90_put_att S")

    ! col
    status = nf90_def_var(ncid, "col", NF90_INT, (/dimid_n_s/), varid_col)
    call netcdf_error(status, "create_mapping_file: nf90_def_var col")

    status = nf90_put_att(ncid, varid_col, "long_name", "column corresponding to matrix elements")
    call netcdf_error(status, "create_mapping_file: nf90_put_att col")

    ! row
    status = nf90_def_var(ncid, "row", NF90_INT, (/dimid_n_s/), varid_row)
    call netcdf_error(status, "create_mapping_file: nf90_def_var row")

    status = nf90_put_att(ncid, varid_row, "long_name", "row corresponding to matrix elements")
    call netcdf_error(status, "create_mapping_file: nf90_put_att row")
    
    ! global attributes
    status = nf90_put_att(ncid, NF90_GLOBAL, "input_mapping_file", trim(filename_mapping_in))
    call netcdf_error(status, "create_mapping_file: nf90_put_att input_mapping_file")

    status = nf90_put_att(ncid, NF90_GLOBAL, "input_mesh_file", trim(filename_mpas_mesh))
    call netcdf_error(status, "create_mapping_file: nf90_put_att input_mesh_file")

    status = nf90_put_att(ncid, NF90_GLOBAL, "smoothing_method", "2D Gaussian smoothing")
    call netcdf_error(status, "create_mapping_file: nf90_put_att smoothing_method")

    status = nf90_put_att(ncid, NF90_GLOBAL, "created_by", "smooth_runoff.exe")
    call netcdf_error(status, "create_mapping_file: nf90_put_att created_by")

    call date_and_time(date, time, zone)
    datetime = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"_"//time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//trim(zone)
    status = nf90_put_att(ncid, NF90_GLOBAL, "created_at", trim(datetime))
    call netcdf_error(status, "create_mapping_file: nf90_put_att created_at")

    status = nf90_put_att(ncid, NF90_GLOBAL, "distanceLimit", distanceLimit)
    call netcdf_error(status, "create_mapping_file: nf90_put_att distanceLimit")

    status = nf90_put_att(ncid, NF90_GLOBAL, "stdDeviation", stdDeviation)
    call netcdf_error(status, "create_mapping_file: nf90_put_att stdDeviation")

    ! end definition phase
    status = nf90_enddef(ncid)
    call netcdf_error(status, "create_mapping_file: nf90_enddef")

    ! write variables
    ! xc_a
    status = nf90_put_var(ncid, varid_xc_a, xc_a)
    call netcdf_error(status, "create_mapping_file: nf90_put_var xc_a")
    
    ! yc_a
    status = nf90_put_var(ncid, varid_yc_a, yc_a)
    call netcdf_error(status, "create_mapping_file: nf90_put_var yc_a")
    
    ! xv_a
    status = nf90_put_var(ncid, varid_xv_a, xv_a)
    call netcdf_error(status, "create_mapping_file: nf90_put_var xv_a")
        
    ! yv_a
    status = nf90_put_var(ncid, varid_yv_a, yv_a)
    call netcdf_error(status, "create_mapping_file: nf90_put_var yv_a")
    
    ! mask_a
    status = nf90_put_var(ncid, varid_mask_a, mask_a)
    call netcdf_error(status, "create_mapping_file: nf90_put_var mask_a")
    
    ! area_a
    status = nf90_put_var(ncid, varid_area_a, area_a)
    call netcdf_error(status, "create_mapping_file: nf90_put_var area_a")
    
    ! frac_a
    status = nf90_put_var(ncid, varid_frac_a, frac_a)
    call netcdf_error(status, "create_mapping_file: nf90_put_var frac_a")
    
    ! src_grid_dims
    status = nf90_put_var(ncid, varid_src_grid_dims, src_grid_dims)
    call netcdf_error(status, "create_mapping_file: nf90_put_var src_grid_dims")
    
    ! xc_b
    status = nf90_put_var(ncid, varid_xc_b, xc_b)
    call netcdf_error(status, "create_mapping_file: nf90_put_var xc_b")
    
    ! yc_b
    status = nf90_put_var(ncid, varid_yc_b, yc_b)
    call netcdf_error(status, "create_mapping_file: nf90_put_var yc_b")
    
    ! xv_b
    status = nf90_put_var(ncid, varid_xv_b, xv_b)
    call netcdf_error(status, "create_mapping_file: nf90_put_var xv_b")
    
    ! yv_b
    status = nf90_put_var(ncid, varid_yv_b, yv_b)
    call netcdf_error(status, "create_mapping_file: nf90_put_var yv_b")
    
    ! mask_b
    status = nf90_put_var(ncid, varid_mask_b, mask_b)
    call netcdf_error(status, "create_mapping_file: nf90_put_var mask_b")
    
    ! area_b
    status = nf90_put_var(ncid, varid_area_b, area_b)
    call netcdf_error(status, "create_mapping_file: nf90_put_var area_b")
    
    ! frac_b
    status = nf90_put_var(ncid, varid_frac_b, frac_b)
    call netcdf_error(status, "create_mapping_file: nf90_put_var frac_b")
    
    ! dst_grid_dims
    status = nf90_put_var(ncid, varid_dst_grid_dims, dst_grid_dims)
    call netcdf_error(status, "create_mapping_file: nf90_put_var dst_grid_dims")

  end subroutine create_mapping_file

  !----------------------------------------------------------------

  subroutine write_mapping_file(ncid, varid_S, varid_col, varid_row, start, count, S_out, col_out, row_out)

    integer, intent(in) :: &
         ncid, &
         varid_S, &
         varid_col, &
         varid_row, &
         start, &
         count

    integer :: &
         status

    integer, dimension(:), intent(in) :: &
         col_out, &
         row_out

    real(kind=RKIND), dimension(:), intent(inout) :: &
         S_out

    integer, dimension(1) :: &
         start1D, &
         count1D

    integer :: n

    start1D(1) = start
    count1D(1) = count

! area weight the remap weights
!    do n = 1,count             
!       s_out(n) = s_out(n) * area_a(col_out(n)) / area_b(row_out(n))
!    enddo

    ! S
    status = nf90_put_var(ncid, varid_S, S_out, start1D, count1D)
    call netcdf_error(status, "create_mapping_file: nf90_put_var S")
    
    ! col
    status = nf90_put_var(ncid, varid_col, col_out, start1D, count1D)
    call netcdf_error(status, "create_mapping_file: nf90_put_var col")
    
    ! row
    status = nf90_put_var(ncid, varid_row, row_out, start1D, count1D)
    call netcdf_error(status, "create_mapping_file: nf90_put_var row")

  end subroutine write_mapping_file
    
  !----------------------------------------------------------------

  subroutine close_mapping_file(ncid)

    integer, intent(in) :: &
         ncid

    integer :: &
         status

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "close_mapping_file: nf90_close ")

  end subroutine close_mapping_file
 
  !----------------------------------------------------------------

  subroutine output_weights_diagnostic(filename, weights)

    character(len=*), intent(in) :: &
         filename

    real(kind=RKIND), dimension(:), intent(in) :: &
         weights

    integer :: &
         status, &
         ncid, &
         dimid, &
         varid

    write(*,*) "Output weights diagnostic..."

    status = nf90_open(trim(filename), NF90_WRITE, ncid)
    call netcdf_error(status, "add_weights_variable: nf90_open")

    status = nf90_inq_dimid(ncid, "nCells", dimid)
    call netcdf_error(status, "add_weights_variable: nf90_inq_dimid")
    
    status = nf90_redef(ncid)
    call netcdf_error(status, "add_weights_variable: nf90_redef")
    
    status = nf90_def_var(ncid, "weights", NF90_DOUBLE, (/dimid/), varid)
    call netcdf_error(status, "add_weights_variable: nf90_def_var weights")
    
    status = nf90_enddef(ncid)
    call netcdf_error(status, "add_weights_variable: nf90_enddef")
    
    status = nf90_put_var(ncid, varid, weights)
    call netcdf_error(status, "add_weights_variable: nf90_put_var weights")

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "add_weights_variable: nf90_close")
    
  end subroutine output_weights_diagnostic
  
  !----------------------------------------------------------------

  subroutine load_mesh_file()

    integer :: &
         status, &
         ncid, &
         dimid, &
         varid

    write(*,*) "Load mesh file..."

    ! open file
    status = nf90_open(trim(filename_mpas_mesh), NF90_NOWRITE, ncid)
    call netcdf_error(status, "load_mesh_file: nf90_open")

    ! nCells
    status = nf90_inq_dimid(ncid, "nCells", dimid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_dimid nCells")
    
    status = nf90_inquire_dimension(ncid, dimid, len=nCells)
    call netcdf_error(status, "load_mesh_file: nf90_inquire_dimension nCells")

    ! nEdges
    status = nf90_inq_dimid(ncid, "nEdges", dimid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_dimid nEdges")
    
    status = nf90_inquire_dimension(ncid, dimid, len=nEdges)
    call netcdf_error(status, "load_mesh_file: nf90_inquire_dimension nEdges")

    ! maxEdges
    status = nf90_inq_dimid(ncid, "maxEdges", dimid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_dimid maxEdges")
    
    status = nf90_inquire_dimension(ncid, dimid, len=maxEdges)
    call netcdf_error(status, "load_mesh_file: nf90_inquire_dimension maxEdges")

    allocate(xCell(nCells))
    allocate(yCell(nCells))
    allocate(zCell(nCells))
    allocate(areaCell(nCells))
    allocate(dcEdge(nEdges))
    allocate(nEdgesOnCell(nCells))
    allocate(cellsOnCell(maxEdges,nCells))
    allocate(edgesOnCell(maxEdges,nCells))

    ! xCell
    status = nf90_inq_varid(ncid, "xCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid xCell")
    
    status = nf90_get_var(ncid, varid, xCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var xCell")

    ! yCell
    status = nf90_inq_varid(ncid, "yCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid yCell")
    
    status = nf90_get_var(ncid, varid, yCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var yCell")
    
    ! zCell
    status = nf90_inq_varid(ncid, "zCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid zCell")
    
    status = nf90_get_var(ncid, varid, zCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var zCell")
        
    ! areaCell
    status = nf90_inq_varid(ncid, "areaCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid areaCell")
    
    status = nf90_get_var(ncid, varid, areaCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var areaCell")

    ! dcEdge
    status = nf90_inq_varid(ncid, "dcEdge", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid dcEdge")
    
    status = nf90_get_var(ncid, varid, dcEdge)
    call netcdf_error(status, "load_mesh_file: nf90_get_var dcEdge")

    ! nEdgesOnCell
    status = nf90_inq_varid(ncid, "nEdgesOnCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid nEdgesOnCell")
    
    status = nf90_get_var(ncid, varid, nEdgesOnCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var nEdgesOnCell")

    ! cellsOnCell
    status = nf90_inq_varid(ncid, "cellsOnCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid cellsOnCell")
    
    status = nf90_get_var(ncid, varid, cellsOnCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var cellsOnCell")

    ! edgesOnCell
    status = nf90_inq_varid(ncid, "edgesOnCell", varid)
    call netcdf_error(status, "load_mesh_file: nf90_inq_varid edgesOnCell")
    
    status = nf90_get_var(ncid, varid, edgesOnCell)
    call netcdf_error(status, "load_mesh_file: nf90_get_var edgesOnCell")

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "load_mesh_file: nf90_close")
    
  end subroutine load_mesh_file
  
  !----------------------------------------------------------------

  subroutine load_mapping_file(filename)

    character(len=*), intent(in) :: &
         filename

    integer :: &
         status, &
         ncid, &
         dimid, &
         varid

    integer :: n
  real(kind=RKIND) :: dx, dy, dtor
  dtor = (4._RKIND * atan(1._RKIND)) / 180._RKIND
  dx = 0.5_RKIND * dtor
  dy = 0.25_RKIND * dtor

    write(*,*) "Load mapping file..."

    ! open file
    status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
    call netcdf_error(status, "load_mapping_file: nf90_open")

    ! n_a
    status = nf90_inq_dimid(ncid, "n_a", dimid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_dimid n_a")

    status = nf90_inquire_dimension(ncid, dimid, len=n_a)
    call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension n_a")

    !! ni_a
    !status = nf90_inq_dimid(ncid, "ni_a", dimid)
    !call netcdf_error(status, "load_mapping_file: nf90_inq_dimid ni_a")

    !status = nf90_inquire_dimension(ncid, dimid, len=ni_a)
    !call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension ni_a")

    !! nj_a
    !status = nf90_inq_dimid(ncid, "nj_a", dimid)
    !call netcdf_error(status, "load_mapping_file: nf90_inq_dimid nj_a")

    !status = nf90_inquire_dimension(ncid, dimid, len=nj_a)
    !call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension nj_a")

    ! nv_a
    status = nf90_inq_dimid(ncid, "nv_a", dimid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_dimid nv_a")

    status = nf90_inquire_dimension(ncid, dimid, len=nv_a)
    call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension nv_a")

    ! src_grid_rank
    status = nf90_inq_dimid(ncid, "src_grid_rank", dimid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_dimid src_grid_rank")

    status = nf90_inquire_dimension(ncid, dimid, len=src_grid_rank)
    call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension src_grid_rank")
    
    ! n_b
    status = nf90_inq_dimid(ncid, "n_b", dimid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_dimid n_b")

    status = nf90_inquire_dimension(ncid, dimid, len=n_b)
    call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension n_b")

    !! ni_b
    !status = nf90_inq_dimid(ncid, "ni_b", dimid)
    !call netcdf_error(status, "load_mapping_file: nf90_inq_dimid ni_b")

    !status = nf90_inquire_dimension(ncid, dimid, len=ni_b)
    !call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension ni_b")

    !! nj_b
    !status = nf90_inq_dimid(ncid, "nj_b", dimid)
    !call netcdf_error(status, "load_mapping_file: nf90_inq_dimid nj_b")

    !status = nf90_inquire_dimension(ncid, dimid, len=nj_b)
    !call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension nj_b")

    ! nv_b
    status = nf90_inq_dimid(ncid, "nv_b", dimid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_dimid nv_b")

    status = nf90_inquire_dimension(ncid, dimid, len=nv_b)
    call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension nv_b")

    ! dst_grid_rank
    status = nf90_inq_dimid(ncid, "dst_grid_rank", dimid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_dimid dst_grid_rank")

    status = nf90_inquire_dimension(ncid, dimid, len=dst_grid_rank)
    call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension dst_grid_rank")
    
    ! n_s
    status = nf90_inq_dimid(ncid, "n_s", dimid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_dimid n_s")

    status = nf90_inquire_dimension(ncid, dimid, len=n_s)
    call netcdf_error(status, "load_mapping_file: nf90_inquire_dimension n_s")

    ! allocations
    allocate(xc_a(n_a))
    allocate(yc_a(n_a))
    allocate(xv_a(nv_a,n_a))
    allocate(yv_a(nv_a,n_a))
    allocate(mask_a(n_a))
    allocate(area_a(n_a))
    allocate(frac_a(n_a))
    allocate(src_grid_dims(src_grid_rank))
    allocate(xc_b(n_b))
    allocate(yc_b(n_b))
    allocate(xv_b(nv_b,n_b))
    allocate(yv_b(nv_b,n_b))
    allocate(mask_b(n_b))
    allocate(area_b(n_b))
    allocate(frac_b(n_b))
    allocate(dst_grid_dims(dst_grid_rank))
    allocate(S(n_s))
    allocate(col(n_s))
    allocate(row(n_s))

    ! xc_a
    status = nf90_inq_varid(ncid, "xc_a", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid xc_a")
    
    status = nf90_get_var(ncid, varid, xc_a)
    call netcdf_error(status, "load_mapping_file: nf90_get_var xc_a")

    ! yc_a
    status = nf90_inq_varid(ncid, "yc_a", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid yc_a")
    
    status = nf90_get_var(ncid, varid, yc_a)
    call netcdf_error(status, "load_mapping_file: nf90_get_var yc_a")

    ! xv_a
    status = nf90_inq_varid(ncid, "xv_a", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid xv_a")
    
    status = nf90_get_var(ncid, varid, xv_a)
    call netcdf_error(status, "load_mapping_file: nf90_get_var xv_a")

    ! yv_a
    status = nf90_inq_varid(ncid, "yv_a", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid yv_a")
    
    status = nf90_get_var(ncid, varid, yv_a)
    call netcdf_error(status, "load_mapping_file: nf90_get_var yv_a")

    ! mask_a
    status = nf90_inq_varid(ncid, "mask_a", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid mask_a")
    
    status = nf90_get_var(ncid, varid, mask_a)
    call netcdf_error(status, "load_mapping_file: nf90_get_var mask_a")

    ! area_a
    status = nf90_inq_varid(ncid, "area_a", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid area_a")
    
    status = nf90_get_var(ncid, varid, area_a)
    call netcdf_error(status, "load_mapping_file: nf90_get_var area_a")

!!hardwire define area_a for rx5 grid
!    do n = 1,n_a
!       area_a(n) = dx * (sin(dtor*yc_a(n)+dy)-sin(dtor*yc_a(n)-dy))
!    enddo

    ! frac_a
    status = nf90_inq_varid(ncid, "frac_a", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid frac_a")
    
    status = nf90_get_var(ncid, varid, frac_a)
    call netcdf_error(status, "load_mapping_file: nf90_get_var frac_a")

    ! src_grid_dims
    status = nf90_inq_varid(ncid, "src_grid_dims", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid src_grid_dims")
    
    status = nf90_get_var(ncid, varid, src_grid_dims)
    call netcdf_error(status, "load_mapping_file: nf90_get_var src_grid_dims")

    ! xc_b
    status = nf90_inq_varid(ncid, "xc_b", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid xc_b")
    
    status = nf90_get_var(ncid, varid, xc_b)
    call netcdf_error(status, "load_mapping_file: nf90_get_var xc_b")

    ! yc_b
    status = nf90_inq_varid(ncid, "yc_b", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid yc_b")
    
    status = nf90_get_var(ncid, varid, yc_b)
    call netcdf_error(status, "load_mapping_file: nf90_get_var yc_b")

    ! xv_b
    status = nf90_inq_varid(ncid, "xv_b", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid xv_b")
    
    status = nf90_get_var(ncid, varid, xv_b)
    call netcdf_error(status, "load_mapping_file: nf90_get_var xv_b")

    ! yv_b
    status = nf90_inq_varid(ncid, "yv_b", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid yv_b")
    
    status = nf90_get_var(ncid, varid, yv_b)
    call netcdf_error(status, "load_mapping_file: nf90_get_var yv_b")

    ! mask_b
    status = nf90_inq_varid(ncid, "mask_b", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid mask_b")
    
    status = nf90_get_var(ncid, varid, mask_b)
    call netcdf_error(status, "load_mapping_file: nf90_get_var mask_b")

    ! area_b
    status = nf90_inq_varid(ncid, "area_b", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid area_b")
    
    status = nf90_get_var(ncid, varid, area_b)
    call netcdf_error(status, "load_mapping_file: nf90_get_var area_b")

    ! frac_b
    status = nf90_inq_varid(ncid, "frac_b", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid frac_b")
    
    status = nf90_get_var(ncid, varid, frac_b)
    call netcdf_error(status, "load_mapping_file: nf90_get_var frac_b")

    ! dst_grid_dims
    status = nf90_inq_varid(ncid, "dst_grid_dims", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid dst_grid_dims")
    
    status = nf90_get_var(ncid, varid, dst_grid_dims)
    call netcdf_error(status, "load_mapping_file: nf90_get_var dst_grid_dims")

    ! S
    status = nf90_inq_varid(ncid, "S", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid S")
    
    status = nf90_get_var(ncid, varid, S)
    call netcdf_error(status, "load_mapping_file: nf90_get_var S")

    ! col
    status = nf90_inq_varid(ncid, "col", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid col")
    
    status = nf90_get_var(ncid, varid, col)
    call netcdf_error(status, "load_mapping_file: nf90_get_var col")

    ! row
    status = nf90_inq_varid(ncid, "row", varid)
    call netcdf_error(status, "load_mapping_file: nf90_inq_varid row")
    
    status = nf90_get_var(ncid, varid, row)
    call netcdf_error(status, "load_mapping_file: nf90_get_var row")

    ! close
    status = nf90_close(ncid)
    call netcdf_error(status, "load_mapping_file: nf90_close")

print*,'col',maxval(col),maxloc(col),size(col)
print*,'col',minval(col),minloc(col)
! area weight the remap weights
    do n = 1,n_s
       s(n) = s(n) * area_a(col(n)) / area_b(row(n))
    enddo
    
  end subroutine load_mapping_file

  !----------------------------------------------------------------

  subroutine netcdf_error(status, message)

    integer, intent(in) :: &
         status

    character(len=*), intent(in) :: &
         message

    if (status /= 0) then
       write(*,*) "Netcdf error: ", status, nf90_strerror(status)
       write(*,*) trim(message)
       stop
    endif
    
  end subroutine netcdf_error
    
  !----------------------------------------------------------------
  
end program smooth_runoff

! gfortran -ffree-line-length-none -I/Users/akt/Work/libraries/netcdf/netcdf/include -L/Users/akt/Work/libraries/netcdf/netcdf/lib -lnetcdff -lnetcdf smooth_runoff.F90 -o smooth_runoff.exe
