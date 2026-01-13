program modify_rx1_mesh
! this program will modify the input mesh file by defining a mask for all cells
!  for which runoff is always zero in all data files.
! A copy of the mesh file must already reside in the current working directory.
! This version is for the 1x1 degree mesh

use netcdf
implicit none

   ! filenames
   integer, parameter :: nf = 2
   character(len=80), parameter :: meshfile =    &
      './ll_1deg_ESMFmesh.nc'
   character(len=128) :: datafile 
   character(len=128), parameter :: file_nyf =    &
      '/glade/p/cesmdata/cseg/inputdata/lnd/dlnd7/RX1/runoff.daitren.annual.20190226.nc'
   character(len=128), parameter :: file_iaf =    &
      '/glade/p/cesmdata/cseg/inputdata/lnd/dlnd7/RX1/runoff.daitren.iaf.20120419.nc'
   character(len=4) :: timevar

   ! grid variables
   integer, parameter :: nlon = 360, nlat = 180
   integer n, i, j, m, nt
   integer, dimension(nlon*nlat) :: mask

   real*8, dimension(nlon,nlat) :: roff

   ! netcdf variables
   integer :: rc, ncid, dimid, varid

   mask = 0
   do n = 1,nf
      ! open a data file
      select case(n)
      case(1)
         datafile = file_nyf  
         timevar = 'nt'
      case(2)
         datafile = file_iaf  
         timevar = 'time'
      end select  

      ! open the data file and get the time dimension and varid
      rc = nf90_open(trim(datafile),nf90_nowrite,ncid)
print*,n,'open data',rc
      rc = nf90_inq_dimid(ncid,trim(timevar),dimid)
print*,n,'inq_dimid data',rc
      rc = nf90_inquire_dimension(ncid,dimid,len=nt)
print*,n,'inq_dimension data',rc
      rc = nf90_inq_varid(ncid,'runoff',varid)
print*,n,'inq_varid data',rc

      ! time loop
      do m = 1,nt
         rc = nf90_get_var(ncid,varid,roff,(/1,1,m/),(/nlon,nlat,1/))
print*,n,m,'get_var data',rc
         do j = 1,nlat
            do i = 1,nlon
               if(roff(i,j) > 0.0) mask(i+(j-1)*nlon) = 1
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

end program modify_rx1_mesh

