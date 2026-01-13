program ll_0p5deg_esmf
! This program creates an ESMF compliant grid file for 0.5 degree longitude/latitude MOSART
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
   integer, parameter :: nlon = 720, nlat = 360
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
   fileout = 'll_0p5deg_ESMFmesh.nc'

   ! define the output file
   rc = nf90_open(trim(fileout),ior(nf90_write,nf90_netcdf4),ncidout)
   print*,'create ',rc

   ! dimensions
   rc = nf90_inq_dimid(ncidout, 'elementCount', dimid_elementCount)
   print*,'inq_dimid ',rc

   ! variables
   rc = nf90_inq_varid(ncidout, 'centerCoords', varid )
   print*,'inq_varid ',rc
   rc = nf90_get_var(ncidout, varid, centercoords )
   print*,'get_var ',rc

   rc = nf90_redef(ncidout)
   print*,'redef ',rc
   rc = nf90_def_var(ncidout, 'elementArea', nf90_double, (/dimid_elementCount/), varid )
   print*,'def_var ',rc
      rc = nf90_put_att(ncidout, varid, 'units', 'radians^2')
      rc = nf90_put_att(ncidout, varid, 'long_name', 'area weights')

   rc = nf90_enddef(ncidout)
   print*,'enddef ',rc

   r2d = 180._8 / (4._8 * atan(1._8))

   ! define the output variables

   n = 0
   do j = 1,nlat
      do i = 1,nlon
         n = n + 1
         dy = sin((centercoords(2,n)+0.25)/r2d) - sin((centercoords(2,n)-0.25)/r2d)
         elementarea(n) = 0.5 * dy / ( r2d )
      enddo
   enddo

   ! write
      rc = nf90_put_var(ncidout,varid, elementArea)
   print*,'put_var ',rc
 
   ! close the files
   rc = nf90_close(ncidout)

   print*,'4*pi, sum(elementArea)',16.*atan(1.),sum(elementArea)

end program ll_0p5deg_esmf
