## **Installation**

The ptha package is most often used on Linux. More recently (2023/05) the
author has successfully installed it on Windows following the standard instructions.

Alternatively, users running Windows or Mac could try using a virtual machine
to run Linux (e.g. via [VirtualBox](https://www.virtualbox.org) or similar).
Before 2023 this was thought to be the only approach that worked.

### **Installing R**
If you don't already have R installed, the you need to get it by following the
instructions on the [R website](https://www.r-project.org/). Use the most recent
version. 

If you are running Windows then you'll also need to 
[install the appropriate version of Rtools](https://cran.r-project.org/bin/windows/Rtools/)

### **Installing netcdf**

You need to install the R package `ncdf4` with OPeNDAP support.

On Linux, this requires pre-installing the netcdf libraries. 
```
sudo apt-get install netcdf-bin libnetcdf-dev libnetcdff-dev
```

On all platforms it should then be possible to install netcdf4 from inside R.
```r
# Run this from inside R.
# Ubuntu users may need to start R with "sudo R" to get permissions to install
install.packages('ncdf4')
```

#### **Troubleshooting netcdf**

For many years it was difficult to get opendap support in Windows (preventing
installation on that platform), although this issue seems to be resolved as of
2023/05. An earlier version of Ubuntu also shipped netcdf with an opendap bug which
caused problems with ptha data access, but this is also resolevd nowadays. 

To confirm that your `ncdf4` installation is using a suitably recent netcdf-c
library, please run the following code:


```r
library(ncdf4)

# This is a file from the PTHA, describing earthquake events on the Kermadec-Tonga
# source-zone. Note I pre-pend [stringlength=4096], which prevents truncation of
# long character strings. This functionality was broken in older netcdf-c versions
test_file = paste0('[stringlength=4096]http://dapds00.nci.org.au/thredds/dodsC/fj6/',
    'PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/',
    'all_stochastic_slip_earthquake_events_kermadectonga2.nc')

# Open it (this will not read everything)
fid = nc_open(test_file, readunlim=FALSE)

# Try to read the event_index_string. This will be artificially truncated if 
# using an old version of netcdf-c
event_index_string = ncvar_get(fid, 'event_index_string')

#
# Report whether it worked.
#

if(max(nchar(event_index_string)) == 756){
    print('Success! Your ncdf4 install can read large character strings remotely')
}else{
    print('FAIL. Perhaps ncdf4 is linking to an older netcdf-c version?')
}

# Shut down the connection 
nc_close(fid)
```
If the above code leads to the `Success! ...` message, then the install is
working. Otherwise you will have to troubleshoot your netcdf-c install (or
your internet connection, or check for a change to the NCI THREDDS server,
etc).

If your internet is not working perfectly, or the NCI server is down, you will see an
message like this:


```r
#    Error in Rsx_nc4_get_vara_double: NetCDF: DAP failure
#    Var: gaugeID  Ndims: 1   Start: 0 Count: 20185
#    Error in ncvar_get_inner(ncid2use, varid2use, nc$var[[li]]$missval, addOffset,  :  
#      C function R_nc4_get_vara_double returned error
```
In this case, just try again -- after a few attempts it usually works. If not,
then check if your internet is working. Also check whether the NCI THREDDS
server is running (occasionally it goes down for maintainence or technical
problems).


### **Installing rptha**

Finally you need to install the `rptha` package. This must be built from
source, after obtaining the code from Geoscience Australia's the `PTHA` github
repository. [See instructions here](https://github.com/GeoscienceAustralia/ptha/blob/master/R/README.md).

### **Installing mapview**
You may optionally install the `mapview` package. This is not critical but can enable some interactive plots shown in the main README.

```r
install.packages('mapview')
```

## **Unit tests**
Assuming you have installed all the above dependencies, you can run the unit
tests with:

```r
source('test_all.R')
```
This should print a number of 'PASS' statements, and no 'FAIL' statements. It might 
take a minute or more depending on your internet connection, because it reads datasets
from the NCI THREDDS server as part of the test.

Some of the tests will fail if you haven't installed the dependencies to read
time-series. Beware some failures can occur if your internet connection is not
performing well (or if there are issues with the NCI THREDDS server). Thus, if
the tests fail it is worth re-trying a few times to confirm it is not due to a
transient network issue.


