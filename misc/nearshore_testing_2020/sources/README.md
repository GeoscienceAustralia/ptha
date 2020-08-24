# How to download a version of this folder with all the source-inversions
-------------------------------------------------------------------------

A compressed version of these files has been posted to the NCI THREDDS Server here: http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/sources.zip

The above link also contains some data that is needed to create the rasters and plots, but is too large to place in this repository. However we keep the codes + basic data here.

Below is the regular documentation.

# Earthquake source inversions
------------------------------

The directories here contain code and data to create vertical deformation
rasters from published earthquake source inversions, and (in most cases) smooth
the result with a Kajiura filter.

*Beware in some cases the folder or file names are misleading -- e.g. the years
are wrong, or we incorrectly spell the author name. For instance the source
model referred to as "RomanoEtAl2015" is actually based on a paper "Romano et
al. (2014)". In our manuscript the naming conventions are correct (e.g. for the
aforementioned case we refer to the "R14" source, not R15). The name
corrections are enforced as required in the plotting code; but GD didn't fix
the file names because the models had already been run and such changes would 
risk breaking the code.* 

# How to run it

## Step 0: Install R and the rptha package

Follow the instructions [here](https://github.com/GeoscienceAustralia/ptha/tree/master/R)

## Step 1: Compute the unfiltered vertical deformation from the earthquake-source-inversion data

To create the un-filtered vertical deformation rasters, use `Rscript` to run each of the following codes
(called from inside their own directory).

```
    Chile1960/FujiSatake2013/Okada_vertical_component.R  
    Chile1960/HoEtAl2019/reconstruct_free_surface.R      

    Chile2010/FujiSatake2013/Okada_vertical_component.R  
    Chile2010/LoritoEtAl2011/Okada_vertical_component.R  

    Chile2015/RomanoEtAl2016/convert_for_TFD.R           
    Chile2015/WilliamsonEtAl2017/Okada_vertical_deformation.R

    Sumatra2004/FujiSatake2007/Okada_vertical_component.R
    Sumatra2004/LoritoEtAl2010/Okada_vertical_component.R
    Sumatra2004/PiatanesiLorito2007/Okada_vertical_component.R

    Tohoku2011/SatakeEtAl2013/Okada_vertical_component.R
    Tohoku2011/YamakaziEtAl2018/Okada_vertical_component.R
    (not required for Tohoku2011/RomanoEtAl2015)
```

For example, to run the file `Chile1960/FujiSatake2013/Okada_vertical_component.R`, start a terminal and do:

```
    # Move into the directory
    cd Chile1960/FujiSatake2013/
    # Run the script -- requires that rptha is installed.
    Rscript Okada_vertical_component.R
```

This will create a number of output files, including a file containing the
vertical deformation `Fuji_chile1960_sources_SUM.tif`. 

For other inversions that have a `Okada_vertical_component.R` script, analogous files are created. In some other cases different approaches apply (i.e. some inversions do not have a `Okada_vertical_component.R` script). Thoses cases are:
* Chile1960/HoEtAl2019/  - the script makes a linear combination of gaussian free-surface perturbations.
* Chile2015/RomanoEtAl2016/ - this inversions uses triangular elements, and relies on the Fortran version of the TFD code (Meade 2007) provided by Fabrizio Romano (not currently included herein -- but if that is installed, then the provided script can produce the deformation raster).
* Tohoku2011/RomanoEtAl2015 - this inversion is a sum of unit-sources derived from a 3D finite-element model. In this case only the raster is provided. (*As mentioned above the name is misleading, I should have called it RomanoEtAl2014, but the names persist throughout this project code so I have chosen not to correct it.*).

## Step 2: Apply Kajiura filter to the vertical deformations (in most cases)

Once all the rasters have been created, the script `apply_kajiura_to_rasters.R` can be run
to apply a Kajiura filter to most of the rasters. 

    Rscript apply_kajiura_to_rasters.R

This is not applicable for Chile1960/HoEtAl2019, which is already a
water-surface deformation - so that is skipped.

# Further background on the source-inversions
---------------------------------------------

Here we give a brief overview of the source inversion techniques used in the original studies, in terms of the following categories:
1. Data
2. Geometry
3. Earthquake Green's functions details
4. Tsunami Green's functions details (e.g. SW --> use of shallow water equations to model tsunami)
5. Inversion (i.e. which one was used in our paper)


## Chile 1960

F13
1. Tsunami waveforms from tide-gauges, coseismic deformation from land-level changes
2. Planar fault.
3. Okada, horizontal displacement.
4. SW.
5. Best slip model from joint inversion (tsunami + geodetic), instantaneous rupture,

H19
1. Tsunami waveforms from tide-gauges, coseismic deformation from land-level changes
2. Planar fault.
3. Okada.
4. SW. 
5. Best slip model from joint inversion (tsunami + geodetic), instantaneous rupture, optimal time alignment + correction.

## Sumatra 2004

P07
1. Tsunami waveforms from tide-gauges.
2. Variable strike and dip.
3. Okada.
4. SW.
5. Best slip model from nonlinear inversion, constant rupture velocity.

L10
1. Tsunami waveforms from tide-gauges and satellite altimetry, coseismic deformation from GPS.
2. Variable strike and dip.
3. Okada+FEM (GPS, ridigity as free parameter).
4. SW.
5. Best slip model from nonlinear inversion, constant rupture velocity.

F07
1. Tsunami waveforms from tide-gauges and satellite altimetry.
2. Variable strike and dip.
3. Okada, horizontal displacement.
4. SW.
5. Best slip model from linear joint inversion of tide-gauges and satellite altimetry, instantaneous rupture.

## Chile 2010

L11:
1. Tsunami waveforms from DART and tide-gauges, coseismic deformation from InSAR, GPS, and land-level changes.
2. Variable strike and dip (rake as free parameter).
3. Okada.
4. SW.
5. Best slip model from nonlinear joint inversion, constant rupture velocity, smoothness regularization, moment minimization.

F13:
1. Tsunami waveforms from DART and tide-gauges, coseismic deformation from land-level changes.
2. Planar.
3. Okada, horizontal displacement.
4. SW.
5. Best slip model from joint inversion (tsunami + geodetic), instantaneous rupture

## Tohoku 2011

S13:
1. Tsunami waveforms from ocean-bottom pressure, DART, GPS buoys, coastal wave gauges, tsunameters and tide gauges.
2. Variable strike and dip.
3. Okada, horizontal displacement. 
4. SW. 
5. Best slip model from Linear, multiple time-window inversion.

R14:
1. Tsunami waveforms from ocean-bottom pressure, DART, GPS buoys, coastal wave gauges, tsunameters, coseismic deformation from inland GPS and sea bottom data.
2. 3D faults.
3. Heterogeneous FEM, horizontal displacement.
4. Non-hydrostatic SW. 
5. Slip model from nonlinear joint inversion (tsunami + geodetic), constant rupture velocity, smoothness regularization, moment minimization.

Y18:
1. Teleseismic waveforms, tsunami waveforms from ocean-bottom pressure, DART, GPS buoys, tsunameters, runup, coseismic deformation from inland GPS and sea bottom data.
2. Variable dip.
3. Okada.
4. Non-hydrostatic SW. 
5. Best slip model, teleseismic inversion, smoothness regularization + near-field iterative, kinematic.

## Chile 2015

R16
1. Tsunami waveforms from tide-gauges and DART, coseismic deformation from InSAR.
2. Variable strike and dip.
3. Meade (a triangular variant of Okada), horizontal displacement.
4. Non-hydrostatic SW.
5. Slip model from nonlinear joint inversion (tsunami + geodetic), constant rupture velocity, smoothness regularization, moment minimization.

W17
1. Tsunami waveforms from tide-gauges and DART buoys, coseismic deformation from GPS and InSAR.
2. Variable strike and dip.
3. Okada.
4. Kajiura, SW.
5. Best slip model from linear joint inversion (tsunami + geodetic), smoothness regularization, instantaneous rupture.
