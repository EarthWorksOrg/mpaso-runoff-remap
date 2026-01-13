program modify_r025_mesh
! this program will modify the input mesh file by defining a mask for all cells
!  for which runoff is always zero in all data files.
! A copy of the mesh file must already reside in the current working directory.
! This version is for the 0.25 degree mesh used by JRA55 data roff.

use netcdf
implicit none

   ! filenames
   character(len=80), parameter :: meshfile =    &
      './ll_025deg_ESMFmesh.nc'
   character(len=128) :: file =    &
      '/glade/p/cesmdata/cseg/inputdata/lnd/dlnd7/JRA55/JRA.v1.1.runoff.xxxx.170807.nc' 
   character(len=4) :: timevar

   ! grid variables
   integer, parameter :: nlon = 1440, nlat = 720
   integer n, i, j, m, nt, year
   integer, dimension(nlon*nlat) :: mask

   real*4, dimension(nlon,nlat) :: rofl, rofi

   ! netcdf variables
   integer :: rc, ncid, dimid, varid, varidl, varidi

   mask = 0
         timevar = 'time'

   do year = 1958,2016
      write(unit=file(66:69),fmt='(i4)')year
      ! open the data file and get the time dimension and varid
      rc = nf90_open(trim(file),nf90_nowrite,ncid)
print*,year,'open data',rc
      rc = nf90_inq_dimid(ncid,trim(timevar),dimid)
print*,year,'inq_dimid data',rc
      rc = nf90_inquire_dimension(ncid,dimid,len=nt)
print*,year,'inq_dimension data',rc
      rc = nf90_inq_varid(ncid,'rofl',varidl)
print*,year,'inq_varid data',rc
      rc = nf90_inq_varid(ncid,'rofi',varidi)
print*,year,'inq_varid data',rc

      ! time loop
      do m = 1,nt
         rc = nf90_get_var(ncid,varidl,rofl,(/1,1,m/),(/nlon,nlat,1/))
         rc = nf90_get_var(ncid,varidi,rofi,(/1,1,m/),(/nlon,nlat,1/))
         rofl = rofl + rofi
print*,year,m,'get_var data',rc
         do j = 1,nlat
            do i = 1,nlon
               if(rofl(i,j) > 0.0) mask(i+(j-1)*nlon) = 1
            enddo
         enddo
      enddo

      ! close the file
      rc = nf90_close(ncid)

   enddo

   ! overwrite mask in meshfile
   rc = nf90_open(trim(meshfile),nf90_write,ncid)
   print*,'open mesh',rc
   rc = nf90_inq_varid(ncid,'elementMask',varid)
   print*,'inqvarid mesh',rc
   rc = nf90_put_var(ncid,varid,mask)
   print*,'putvar mesh',rc
   rc = nf90_close(ncid)
   print*,'close mesh',rc

end program modify_r025_mesh

