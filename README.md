# mpaso-runoff-remap
Programs and documentation for creating runoff remap files used to remap terrestrial runoff to the MPAS ocean in EarthWorks

There are three possible runoff models for use in EarthWorks - 0.5 degree prognostic MOSART; and two data runoff models - 1 degree DROF for used with DATM%NYF and DATM%IAF, and 0.25 JRA55 runoff.

Do these steps on a machine with access to the runoff data and ESMF command line functions:

1) Build a local ESMF mesh file that mimics the runoff source grid. The ll_*.f90 programs do this. Edit paths as necessary.

2) Modify the local ESMF mesh file by masking those cells that contribute no runoff. The modify_*.f90 programs do this for the data runoff models. Again, edit paths as necessary.

3) Generate an initial remap file. the interp_*.pbs scripts do this. These call the ESMF_RegridWeightGen executable and generates a nearest neighbor remapping to the MPAS-ocean mesh. The user specifies the MPAS-ocean mesh and the remap file to be created.

The next steps you will want to do on a local machine as some of the steps can consume a lot of wall clock at higher resolutions and you want to avoid batch wall clock limits.

Transfer to your local machine the mpas-ocean esmf mesh file, the initial generated remap file, and the ocean domain mesh file (referred to in the mpas-ocean-inputdata repository)