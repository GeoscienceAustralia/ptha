suppressPackageStartupMessages(library(geosphere))
# This shows how to replace a function in a package with something else

.antipodal_modified<-function (p1, p2, tol = 1e-09) 
{
    p1 <- suppressWarnings(geosphere:::.pointsToMatrix(p1))
    p2 <- suppressWarnings(geosphere:::.pointsToMatrix(p2))
    p <- cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2])
    p[, c(1, 3)] <- suppressWarnings(geosphere:::.normalizeLonDeg(p[, c(1, 3)]))
    diflon <- abs(p[, 1] - p[, 3])
    diflat <- abs(p[, 2] + p[, 4])
    (diflat < tol) & (abs(diflon%%360 - 180) < tol)
}

.onLoad<-function(lib, pkg){
    library(utils)
    suppressPackageStartupMessages(library(geosphere))

    unlockBinding('antipodal', as.environment('package:geosphere'))

    assignInNamespace('antipodal', .antipodal_modified, 
        ns = 'geosphere', envir = as.environment('package:geosphere'))

    lockBinding('antipodal', as.environment('package:geosphere'))
}


