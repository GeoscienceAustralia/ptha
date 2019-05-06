#
# To get the PTHA package running with R-3.5.1 on NCI, we run the following sequence
# of commands (note the comment where some manual software patching is required)
#

install.packages('sp')
install.packages('Rcpp')
install.packages('raster')
install.packages('rgdal')
install.packages('rgeos')
install.packages('FNN')
install.packages('minpack.lm')
#
## This failed
# install.packages('geometry')
## because of a bug in that package with the intel compiler
## The work-around is to download the package source, patch it as described below, then
## do a local install (i.e. R CMD INSTALL geometry_XX.tar.gz, where you have
## patched the source in the latter .gz file)
##
##    ## HOW TO PATCH THE GEOMETRY PACKAGE FOR INTEL-COMPILER
##    ## THIS IS A COPY OF A BUG REPORT I FOUND ON THE INTERNET
##    Date:
##    2017-05-11 16:41  Priority:
##    3
##    State:
##    Open
##    Submitted by:
##    Miron Kursa (mbq)         Assigned to:
##    Nobody (None)
##    Hardware:
##    None      Product:
##    None
##    Operating System:
##    None      Component:
##    None
##    Version:
##    None      Severity:
##    None
##    Resolution:
##    None
##    URL:
##    Summary:
##    Geometry does not compile on modern Intel compilers (2016+)
##
##    Detailed description
##    Geometry does not compile on modern Intel compilers (2016+) due to a qhull compilation override which does not apply any more. The problem is in qhull_a.h, lines 105+
##
##    #if defined(__INTEL_COMPILER) && !defined(QHULL_OS_WIN)
##    template <typename T>
##    inline void qhullUnused(T &x) { (void)x; }
##    # define QHULL_UNUSED(x) qhullUnused(x);
##    #else
##    # define QHULL_UNUSED(x) (void)x;
##    #endif
##
##    throws
##
##    In file included from Rconvhulln.c(31):
##    qhull_a.h(106): warning #77: this declaration has no storage class or type specifier
##    template <typename T>
##    ^
##
##    In file included from Rconvhulln.c(31):
##    qhull_a.h(106): error: expected a ";"
##    template <typename T>
##    ^
##
##    In file included from Rconvhulln.c(31):
##    qhull_a.h(115): warning #12: parsing restarts here after previous syntax error
##    void qh_qhull(void);
##    ^
##
##    compilation aborted for Rconvhulln.c (code 2)
##    make: *** [Rconvhulln.o] Error 2
##    ERROR: compilation failed for package �~@~Xgeometry�~@~Y
##
##
##    replacing this with just
##
##    # define QHULL_UNUSED(x) (void)x;
##
##    fixes the problem.
##
##    Probably related qhull issue: https://github.com/qhull/qhull/issues/16

install.packages('geosphere')
install.packages('ncdf4')
install.packages('testthat')

#
# I had problems getting roxygen2 and devtools to build. However, they are
# not required to build the package, if you already have a rptha_XXXX.tar.gz
# file that has been created on another machine. That is a practical compromise
#

## This one is taking a long time and fails because of dependency failures
# install.packages('devtools')
## This is also taking ages, and failed because of dependency issues
# install.packages('roxygen2')


