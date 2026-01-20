PROGRAM modify_dst

   USE netcdf

   IMPLICIT NONE

   ! remap file contents
   CHARACTER(len=128) :: file, filempas
   REAL*8, DIMENSION(:), ALLOCATABLE :: wgt
   INTEGER, DIMENSION(:), ALLOCATABLE :: dst, src
   INTEGER :: n_s, ncells

   ! netcdf variables
   INTEGER :: rc, ncid, varid, dimid, ncidmpas

   ! mpas grid variables
   real*8 dist, distmin, lon0,lat0,coslat0
   real*8, allocatable, dimension(:) :: lon,lat

   ! work and loop variables
   INTEGER n, i, j, nsrc, iii, nlat, nlon, dst_black

   integer, dimension(4) :: slim, nlim, elim, wlim

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   ! 	begin executable code
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   print*,'enter 1 for 1x1 drof, 2 for mosart (0.5 degree), 3 for JRA drof (0.25 degree)'
   read(5,*)iii
   
   print*,'enter the initial remap weights file name'
   read(5,'(a128)')file
   
   ! Define the black sea box on the runoff grid
   IF(iii==1)THEN
      slim(1)=131
      nlim(1) = 139
      wlim(1) = 207
      elim(1) = 222
   ELSEIF(iii==2) THEN
      slim(1)=262
      nlim(1) = 278
      wlim(1) = 414
      elim(1) = 444
   ELSEIF(iii==3) THEN
      slim(1)=524
      nlim(1) = 556
      wlim(1) = 108
      elim(1) = 168
   ELSE
      STOP ' bad input'
   ENDIF
      nlat = 180 * (2**(iii-1))
      nlon = 360 * (2**(iii-1))

   ! find the destination index for Black Sea runoff on the mpas grid (mouth of Dardenelles)
   print*,'enter the mpas mesh variables file name'
   read(5,'(a128)')filempas
   rc = nf90_open(TRIM(filempas),nf90_nowrite,ncidmpas)
   rc = nf90_inq_dimid(ncidmpas,'nCells',dimid)
   rc = nf90_inquire_dimension(ncidmpas,dimid,len=ncells)
   print*,'ncells',ncells

   allocate(lon(ncells))
   allocate(lat(ncells))

   lat0=40.02
   lon0=26.14
   coslat0 = cos(lat0*3.1415926/180.)

   rc = nf90_inq_varid(ncidmpas,'lonCell',varid)
   rc = nf90_get_var(ncidmpas,varid,lon)
   rc = nf90_inq_varid(ncidmpas,'latCell',varid)
   rc = nf90_get_var(ncidmpas,varid,lat)

   lon = lon * 180./ 3.1415926
   lat = lat * 180./ 3.1415926

   distmin = 999.
   do n = 1,ncells
      dist = ((lon0-lon(n))*coslat0)**2 + (lat0-lat(n))**2
      if(dist < distmin) then
         dst_black = n
         distmin = dist
     endif
   enddo
   rc = nf90_close(ncidmpas)
   
   ! now modify the destination vector in the remap file
   rc = nf90_open(TRIM(file),nf90_write,ncid)
   rc = nf90_inq_dimid(ncid,'n_s',dimid)
   rc = nf90_inquire_dimension(ncid,dimid,len=n_s)
   ALLOCATE(wgt(n_s))
   ALLOCATE(src(n_s))
   ALLOCATE(dst(n_s))
   rc = nf90_inq_varid(ncid,'S',varid)
   rc = nf90_get_var(ncid,varid,wgt)
   rc = nf90_inq_varid(ncid,'col',varid)
   rc = nf90_get_var(ncid,varid,src)
   rc = nf90_inq_varid(ncid,'row',varid)
   rc = nf90_get_var(ncid,varid,dst)

   DO n = 1,n_s
      nsrc = src(n)
      j = (nsrc-1) / nlon + 1
      i = nsrc - (j-1) * nlon

         IF(slim(1) .LE. j .AND. j .LE. nlim(1)) THEN
            IF(wlim(1) .LE. i .AND. i .LE. elim(1)) THEN
               dst(n) = dst_black
            ENDIF
         ENDIF

   ENDDO

   rc = nf90_inq_varid(ncid,'row',varid)
   rc = nf90_put_var(ncid,varid,dst)
   print*,'putvar ',rc

END PROGRAM modify_dst
