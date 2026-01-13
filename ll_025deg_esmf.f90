program ll_025deg_esmf
! This program creates an ESMF compliant grid file for 1 degree longitude/latitude JRA55 ROF
!

use netcdf

implicit none

! variable declarations

   ! netcdf variables
   integer ncidin, ncidout, rc
   integer dimid, dimid_nodecount, dimid_elementcount, &
           dimid_maxnodepelement, dimid_coorddim
   integer varid
   character(len=80) :: fileout

   ! dimensions
   integer, parameter :: nlon = 1440, nlat = 720
   integer, parameter :: nodeCount = nlon * (nlat+1)
   integer, parameter :: elementCount = nlon * nlat 
   integer, parameter :: maxNodePElement = 4
   integer, parameter :: coordDim = 2

   ! node lats and lons
   REAL*8, dimension(nlon+1) :: lonv
   REAL*8, dimension(nlat+1) :: latv

   ! output arrays
   real*8,  dimension(coorddim,nodecount) :: nodeCoords
   real*8,  dimension(coorddim,elementcount) :: centerCoords
   real*8,  dimension(elementcount) :: elementArea
   integer, dimension(elementcount) :: elementMask
   integer, dimension(maxnodepelement,elementcount) :: elementConn
   integer*1, dimension(elementcount) :: numElementConn

   ! work variables
   integer i, j, n, ncon
   real*8 r2d, dy

   integer m, ne, nlen, k

! begin executable code

   ! assign filenames, open the input and create the output files
   fileout = 'll_025deg_ESMFmesh.nc'

   ! define the output file
   rc = nf90_create(trim(fileout),ior(nf90_clobber,nf90_netcdf4),ncidout)
   print*,'create ',rc

   ! dimensions
   rc = nf90_def_dim(ncidout, 'nodeCount', nodeCount, dimid_nodeCount)
   print*,'def_dim ',rc
   rc = nf90_def_dim(ncidout, 'elementCount', elementCount, dimid_elementCount)
   print*,'def_dim ',rc
   rc = nf90_def_dim(ncidout, 'maxNodePElement', maxNodePElement, dimid_maxNodePElement)
   print*,'def_dim ',rc
   rc = nf90_def_dim(ncidout, 'coordDim', coordDim, dimid_coordDim)
   print*,'def_dim ',rc

   ! variables
   rc = nf90_def_var(ncidout, 'nodeCoords', nf90_double, (/dimid_coordDim,dimid_nodeCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'units', 'degrees')
   rc = nf90_def_var(ncidout, 'elementConn', nf90_int, (/dimid_maxNodePElement,dimid_elementCount/), varid )
      rc = nf90_put_att(ncidout, varid, 'long_name', 'Node Indices that define the element connectivity')
      rc = nf90_put_att(ncidout, varid, '_FillValue', -1)
   rc = nf90_def_var(ncidout, 'numElementConn', nf90_byte, (/dimid_elementCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'long_name', 'Number of nodes per element')
   rc = nf90_def_var(ncidout, 'centerCoords', nf90_double, (/dimid_coordDim,dimid_elementCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'units', 'degrees')
   rc = nf90_def_var(ncidout, 'elementArea', nf90_double, (/dimid_elementCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'units', 'radians^2')
      rc = nf90_put_att(ncidout, varid, 'long_name', 'area weights')
   rc = nf90_def_var(ncidout, 'elementMask', nf90_int, (/dimid_elementCount/), varid )
   print*,'def_var ',rc

   ! global attributes
   rc = nf90_put_att(ncidout, nf90_global, 'gridType', 'unstructured')
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'version', '1.0')
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'inputFile', 'none' )
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'timeGenerated', '08/26/22')
   print*,'put_att ',rc
   rc = nf90_put_att(ncidout, nf90_global, 'history', 'created by Dazlich ll_0p25_esmfmesh program')
   print*,'put_att ',rc

   rc = nf90_enddef(ncidout)

   r2d = 180._8 / (4._8 * atan(1._8))
   ! compute nodes lats and lon
   do j =1,nlat+1
      latv(j) = -90. + 0.25*(j-1)
   enddo
   do i =1,nlon+1
      lonv(i) = 0.25*(i-1)
   enddo

   ! define the output variables
   elementmask = 1
   numelementconn = 4

   n=0
   do j = 1,nlat+1
      do i = 1,nlon
         n = n + 1
         nodecoords(1,n) = lonv(i)
         nodecoords(2,n) = latv(j)
      enddo
   enddo
   n = 0
   do j = 1,nlat
      dy = sin(latv(j+1)/r2d) - sin(latv(j)/r2d)
      do i = 1,nlon
         n = n + 1
         centercoords(1,n) = 0.5 * (lonv(i) + lonv(i+1))
         centercoords(2,n) = 0.5 * (latv(j) + latv(j+1))
         elementarea(n) = 0.25 * dy / ( r2d )
         elementconn(1,n) = n ! sw vertex
         elementconn(2,n) = nlon * (j-1) + mod(i,nlon) + 1 ! se vertex
         elementconn(3,n) = nlon * j + mod(i,nlon) + 1 ! ne vertex
         elementconn(4,n) = n + nlon ! nw vertex
      enddo
   enddo

   ! write
   rc = nf90_inq_varid(ncidout,'nodeCoords',varid)
      rc = nf90_put_var(ncidout,varid, nodeCoords)
   rc = nf90_inq_varid(ncidout,'elementConn',varid)
      rc = nf90_put_var(ncidout,varid, elementConn)
   rc = nf90_inq_varid(ncidout,'numElementConn',varid)
      rc = nf90_put_var(ncidout,varid, numElementConn)
   rc = nf90_inq_varid(ncidout,'centerCoords',varid)
      rc = nf90_put_var(ncidout,varid, centerCoords)
   rc = nf90_inq_varid(ncidout,'elementArea',varid)
      rc = nf90_put_var(ncidout,varid, elementArea)
   rc = nf90_inq_varid(ncidout,'elementMask',varid)
      rc = nf90_put_var(ncidout,varid, elementMask)
 
   ! close the files
   rc = nf90_close(ncidin)
   rc = nf90_close(ncidout)

   print*,'4*pi, sum(elementArea)',16.*atan(1.),sum(elementArea)

end program ll_025deg_esmf
