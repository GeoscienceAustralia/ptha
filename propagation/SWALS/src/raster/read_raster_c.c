#include "gdal.h"
#include "cpl_conv.h" /* for CPLMalloc() */

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Routines here for object based fortran interface
//

void open_gdal_raster(char inputFile[], GDALDatasetH *hDataset){

    // Give 30MB cache (note this is for GDAL as a whole, not just this dataset.)
    int cache_size_bytes = 1024*1024*30;
    // Try to avoid too much cache use
    GDALSetCacheMax(cache_size_bytes);

    GDALAllRegister();

    *hDataset = GDALOpen( inputFile , GA_ReadOnly); 

    if( *hDataset == NULL )
    {
        printf("Failed to read input file \n");
        exit(0);
    }

    //*hDriver = NULL; //GDALGetDatasetDriver( *hDataset );
    
    return;
}

void close_gdal_raster(GDALDatasetH *hDataset){

    GDALClose(*hDataset);

    // FIXME: How to delete these properly? Is it even required?
    *hDataset = NULL;
    return;

}

// Get basic metadata of raster
void get_gdal_raster_dimensions(GDALDatasetH *hDataset, 
    int xydim[], double lowerleft[], double upperright[], double adfGeoTransform[], double dx[],
    double* nodata_value){

    GDALRasterBandH hBand = GDALGetRasterBand( *hDataset, 1);
    double xleft, ytop;
    int* pbSuccess;

    xydim[0] = GDALGetRasterXSize(*hDataset);
    xydim[1] = GDALGetRasterYSize(*hDataset);

    *nodata_value = GDALGetRasterNoDataValue(hBand, pbSuccess);
    //printf("%e", *nodata_value);

    if( GDALGetGeoTransform( *hDataset, adfGeoTransform ) == CE_None )
    {
        xleft = adfGeoTransform[0];
        ytop = adfGeoTransform[3];

        dx[0] = adfGeoTransform[1];
        dx[1] = adfGeoTransform[5];

        lowerleft[0] = xleft;
        lowerleft[1] = ytop + dx[1]*xydim[1];

        upperright[0] = xleft + dx[0]*xydim[0];
        upperright[1] = ytop;
    }

    return;
}

void get_values_at_xy(GDALDatasetH *hDataset, double adfGeoTransform[], 
    double x[], double y[], double z[], int N, int verbose, int bilinear, int band){

    GDALRasterBandH hBand = GDALGetRasterBand( *hDataset, band);
    
    double *pafScanline;
    int   nXSize = GDALGetRasterBandXSize( hBand );
    int   nYSize = GDALGetRasterBandYSize( hBand );
    pafScanline = (double *) CPLMalloc(sizeof(double));

    int* pbSuccess;
    double nodata_value;

    double xleft, ytop, dx, dy;
    xleft = adfGeoTransform[0];
    ytop = adfGeoTransform[3];
    dx = adfGeoTransform[1];
    dy = adfGeoTransform[5];

    int j, i, xindex, yindex, xoff, yoff;
    double xweight, yweight;
    double z00, z01, z10, z11;
    double z_interp_x0, z_interp_x1;
    double x_local, y_local;
    CPLErr Err; 

    // Beware of hitting nodata
    nodata_value = GDALGetRasterNoDataValue(hBand, pbSuccess);

    //printf("Bilinear interpolation C %d \n", bilinear); 
    // j is the index of x,y,z
    for ( j = 0; j < N; j++){

        x_local = x[j];
        y_local = y[j];

        // Deliberately convert float to integer (rounds towards 0)
        xindex = (x_local - xleft)/dx ;
        yindex = (y_local - ytop)/dy ;

        // Read a line at xindex, yindex
        Err = GDALRasterIO(hBand, GF_Read, xindex, yindex, 1, 1, 
                      pafScanline, 1, 1, GDT_Float64, 
                      0, 0);

        if(bilinear == 1){
            // Bilinear interpolation
            z00 = *(pafScanline);
        
            // For bilinear, we need to read 4 cells, going in either the + or - x
            // and y directions. Make 'xoff' and 'yoff' record these offsets
            if( (x_local - xleft) > (xindex + 0.5) * dx){
                xoff = 1;
            }else{
                xoff = -1;
            }
    
            // Deal with boundary cases where we can't interpolate
            if((xindex == 0) & xoff < 0){
                xoff = 0;
            }
            if((xindex == nXSize-1) & xoff > 0){
                xoff = 0;
            }

            if( (y_local - ytop) < (yindex + 0.5) * dy){
                yoff = 1;
            }else{
                yoff = -1;
            }

            // Deal with boundary cases where we can't interpolate
            if((yindex == 0) & yoff < 0){
                yoff = 0;
            }
            if((yindex == nYSize-1) & yoff > 0){
                yoff = 0;
            }

            // Get weights for bilinear interpolation. 
            xweight = (x_local - xleft - (xindex * dx))/dx ;
            if(xoff == 1){
                xweight = xweight - 0.5;
            }else if(xoff == -1){
                xweight = 0.5 - xweight;
            }else{
                xweight = 1.0;
            }

            yweight = (y_local - ytop  - (yindex * dy))/dy ;
            if(yoff == 1){
                yweight = yweight - 0.5;
            }else if(yoff == -1){
                yweight = 0.5 - yweight;
            }else{
                yweight = 1.0;
            }

            //if(xweight > 1.0) xweight = 1.0 ;
            //if(xweight < 0.0) xweight = 0.0 ;
            //if(yweight > 1.0) yweight = 1.0 ;
            //if(yweight < 0.0) yweight = 0.0 ;

            if( (xweight > 1.0 + 1.0e-10) | (xweight < 0.0 - 1.0e-10) | (yweight > 1.0 + 1.0e-10) | (yweight < 0.0 - 1.0e-10)){
                printf("Bilinear interpolation C error. xw: %e, yw: %e, x_local: %e, y_local: %e, xleft: %e, ytop: %e \n", 
                    xweight, yweight, x_local, y_local, xleft, ytop);
                fflush(stdout);
            }

            // Read the line at xindex + 1, yindex
            Err = GDALRasterIO(hBand, GF_Read, xindex+xoff, yindex, 1, 1, 
                          pafScanline, 1, 1, GDT_Float64, 
                          0, 0);
            z10 = *(pafScanline);

            // Read the line at xindex, yindex + 1
            Err = GDALRasterIO(hBand, GF_Read, xindex, yindex+yoff, 1, 1, 
                          pafScanline, 1, 1, GDT_Float64, 
                          0, 0);
            z01 = *(pafScanline);
            
            // Read the line at xindex + 1, yindex + 1
            Err = GDALRasterIO(hBand, GF_Read, xindex+xoff, yindex+yoff, 1, 1, 
                          pafScanline, 1, 1, GDT_Float64, 
                          0, 0);
            z11 = *(pafScanline);

            // Treat nodata cleanly
            if(z00 == nodata_value | z01 == nodata_value | z10 == nodata_value | z11 == nodata_value){
                // If this is not 'nodata' we might get something useful
                z[j] = z00 ;

            }else{
                z_interp_x0 = z00 * (1.0-xweight) + z10 * xweight;
                z_interp_x1 = z01 * (1.0-xweight) + z11 * xweight;

                z[j] = z_interp_x0 * (1.0 - yweight) + z_interp_x1 * yweight;
            }

            //printf("xleft: %e, ytop: %e, dx: %e, dy: %e, x: %e, y: %e, xR: %e, yR: %e, xw: %e, y: %e, xI: %d, yI: %d\n", xleft, ytop, dx, dy, x[j], y[j], xindex*dx + xleft, yindex*dy + ytop, xweight, yweight, xindex, yindex);
        }else{
            // Nearest cell interpolation
            z[j] = *(pafScanline);
        }


        if(z[j] != z[j]){
            if(verbose){
                printf("NA values found in raster \n");
            }
            exit(0);
        }

    }

    CPLFree(pafScanline) ;

    return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Routines here for non-object based interface
//

//
// Get the raster dimensions, upper left, and lower right
//
void get_raster_dimensions(char inputFile[], int xydim[], double lowerleft[], double upperright[]){

    GDALDatasetH  hDataset;
  
    GDALAllRegister();

    // Open the raster
    hDataset = GDALOpen( inputFile , GA_ReadOnly );

    if( hDataset == NULL )
    {
        printf("Failed to read input file \n");
        exit(0);
    }

    //GDALDriverH   hDriver;
    double        adfGeoTransform[6], xleft, ytop, dx, dy;

    //hDriver = GDALGetDatasetDriver( hDataset );

    xydim[0] = GDALGetRasterXSize( hDataset );
    xydim[1] = GDALGetRasterYSize( hDataset );

    if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None )
    {
        xleft = adfGeoTransform[0];
        ytop = adfGeoTransform[3];

        dx = adfGeoTransform[1];
        dy = adfGeoTransform[5];

        lowerleft[0] = xleft;
        lowerleft[1] = ytop + dy*xydim[1];

        upperright[0] = xleft + dx*xydim[0];
        upperright[1] = ytop;
    }

    GDALClose(hDataset);
    return;
}

//
// Find values of 'z' at 'x,y' points in the raster, either using bilinear interpolation or nearest neighbours
//
void read_gdal_raster(char inputFile[], double x[], double y[], double z[], int N, int verbose, int bilinear){

    GDALDatasetH  hDataset;
  
    GDALAllRegister();

    // Open the raster
    hDataset = GDALOpen( inputFile , GA_ReadOnly );

    if( hDataset == NULL )
    {
        printf("Failed to read input file \n");
        exit(0);
    }

   
    // Tutorial part 2
    // Get driver, size, projection information
 
    GDALDriverH   hDriver;
    double        adfGeoTransform[6], xleft, ytop, dx, dy;

    hDriver = GDALGetDatasetDriver( hDataset );

    if( verbose == 1){
        printf( "Driver: %s/%s\n",
                GDALGetDriverShortName( hDriver ),
                GDALGetDriverLongName( hDriver ) );

        printf( "Size is %dx%dx%d\n",
                GDALGetRasterXSize( hDataset ), 
                GDALGetRasterYSize( hDataset ),
                GDALGetRasterCount( hDataset ) );

        if( GDALGetProjectionRef( hDataset ) != NULL )
            printf( "Projection is `%s'\n", GDALGetProjectionRef( hDataset ) );
    }

    if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None )
    {
        xleft = adfGeoTransform[0];
        ytop = adfGeoTransform[3];

        dx = adfGeoTransform[1];
        dy = adfGeoTransform[5];

        if(verbose){
            printf( "Origin = (%.6f,%.6f)\n", xleft, ytop );
            printf( "Pixel Size = (%.6f,%.6f)\n", dx, dy);
        }
    }

    

    // Tutorial part 3 
    // Get some info on the raster size , min, max

    GDALRasterBandH hBand;
    int             nBlockXSize, nBlockYSize;
    int             bGotMin, bGotMax;
    double          adfMinMax[2];
   
    // Block size information 
    hBand = GDALGetRasterBand( hDataset, 1 );
    GDALGetBlockSize( hBand, &nBlockXSize, &nBlockYSize );

    if(verbose){
        printf( "Block=%dx%d Type=%s, ColorInterp=%s\n",
            nBlockXSize, nBlockYSize,
            GDALGetDataTypeName(GDALGetRasterDataType(hBand)),
            GDALGetColorInterpretationName(
                GDALGetRasterColorInterpretation(hBand)) );
    }

    // Min/max info
    adfMinMax[0] = GDALGetRasterMinimum( hBand, &bGotMin );
    adfMinMax[1] = GDALGetRasterMaximum( hBand, &bGotMax );
    if( ! (bGotMin && bGotMax) )
        GDALComputeRasterMinMax( hBand, TRUE, adfMinMax );

    if(verbose){
        printf( "Min=%.3fd, Max=%.3f\n", adfMinMax[0], adfMinMax[1] );
   
        // Overviews?? 
        if( GDALGetOverviewCount(hBand) > 0 )
            printf( "Band has %d overviews.\n", GDALGetOverviewCount(hBand));

        // Color Tables ??
        if( GDALGetRasterColorTable( hBand ) != NULL )
            printf( "Band has a color table with %d entries.\n", 
                     GDALGetColorEntryCount(
                         GDALGetRasterColorTable( hBand ) ) );

    }
    // Tutorial part 4, with some modifications by GD
    // Read raster values
    // We loop over the entire raster and print non-na values + indices
    // Note we are reading entire blocks (efficient)

    double *pafScanline;
    int   nXSize = GDALGetRasterBandXSize( hBand );
    int   nYSize = GDALGetRasterBandYSize( hBand );
    
    if(verbose){
        printf("nXSize %d \n", nXSize);
        printf("nYSize %d \n", nYSize);
    }

    //pafScanline = (float *) CPLMalloc(sizeof(float)*nXSize);
    //pafScanline = (float *) malloc(sizeof(float)*nXSize);
    pafScanline = (double *) CPLMalloc(sizeof(double));

    int j, i, xindex, yindex, xoff, yoff;
    double xweight, yweight;
    double z00, z01, z10, z11;
    double z_interp_x0, z_interp_x1;
    double x_local, y_local;
    CPLErr Err; 
    // j = 0 corresponds to North
    // i = 0 corresponds to West
    //for ( j = 0; j < nYSize; j++){
        //GDALRasterIO(hBand, GF_Read, 0, j, nXSize, 1, 
        //              pafScanline, nXSize, 1, GDT_Float32, 
        //              0, 0);

    // Try to avoid too much cache use
    GDALSetCacheMax(1048576);

    // j is the index of x,y,z
    for ( j = 0; j < N; j++){

        x_local = x[j];
        y_local = y[j];

        // Deliberately convert float to integer (rounds towards 0)
        xindex = (x_local - xleft)/dx ;
        yindex = (y_local - ytop)/dy ;

        // Read a line at xindex, yindex
        Err = GDALRasterIO(hBand, GF_Read, xindex, yindex, 1, 1, 
                      pafScanline, 1, 1, GDT_Float64, 
                      0, 0);

        if(bilinear == 1){
            // Bilinear interpolation
            z00 = *(pafScanline);
        
            // For bilinear, we need to read 4 cells, going in either the + or - x
            // and y directions. Make 'xoff' and 'yoff' record these offsets
            if( (x_local - xleft) > (xindex + 0.5) * dx){
                xoff = 1;
            }else{
                xoff = -1;
            }
    
            // Deal with boundary cases where we can't interpolate
            if((xindex == 0) & xoff < 0){
                xoff = 0;
            }
            if((xindex == nXSize-1) & xoff > 0){
                xoff = 0;
            }

            if( (y_local - ytop) < (yindex + 0.5) * dy){
                yoff = 1;
            }else{
                yoff = -1;
            }

            // Deal with boundary cases where we can't interpolate
            if((yindex == 0) & yoff < 0){
                yoff = 0;
            }
            if((yindex == nYSize-1) & yoff > 0){
                yoff = 0;
            }

            // Get weights for bilinear interpolation. 
            xweight = (x_local - xleft - (xindex * dx))/dx ;
            if(xoff == 1){
                xweight = xweight - 0.5;
            }else if(xoff == -1){
                xweight = 0.5 - xweight;
            }else{
                xweight = 1.0;
            }

            yweight = (y_local - ytop  - (yindex * dy))/dy ;
            if(yoff == 1){
                yweight = yweight - 0.5;
            }else if(yoff == -1){
                yweight = 0.5 - yweight;
            }else{
                yweight = 1.0;
            }

            if( (xweight > 1.0 + 1.0e-10) | (xweight < 0.0 - 1.0e-10) | (yweight > 1.0 + 1.0e-10) | (yweight < 0.0 - 1.0e-10)){
                printf("Bilinear interpolation C error. xw: %e, yw: %e\n", xweight, yweight);
            }

            // Read the line at xindex + 1, yindex
            Err = GDALRasterIO(hBand, GF_Read, xindex+xoff, yindex, 1, 1, 
                          pafScanline, 1, 1, GDT_Float64, 
                          0, 0);
            z10 = *(pafScanline);

            // Read the line at xindex, yindex + 1
            Err = GDALRasterIO(hBand, GF_Read, xindex, yindex+yoff, 1, 1, 
                          pafScanline, 1, 1, GDT_Float64, 
                          0, 0);
            z01 = *(pafScanline);
            
            // Read the line at xindex, yindex + 1
            Err = GDALRasterIO(hBand, GF_Read, xindex+xoff, yindex+yoff, 1, 1, 
                          pafScanline, 1, 1, GDT_Float64, 
                          0, 0);
            z11 = *(pafScanline);

            z_interp_x0 = z00 * (1.0-xweight) + z10 * xweight;
            z_interp_x1 = z01 * (1.0-xweight) + z11 * xweight;

            z[j] = z_interp_x0 * (1.0 - yweight) + z_interp_x1 * yweight;

            //printf("xleft: %e, ytop: %e, dx: %e, dy: %e, x: %e, y: %e, xR: %e, yR: %e, xw: %e, y: %e, xI: %d, yI: %d\n", xleft, ytop, dx, dy, x[j], y[j], xindex*dx + xleft, yindex*dy + ytop, xweight, yweight, xindex, yindex);
        }else{
            // Nearest cell interpolation
            z[j] = *(pafScanline);
        }


        if(z[j] != z[j]){
            if(verbose){
                printf("NA values found in raster \n");
            }
            exit(0);
        }

    }

    CPLFree(pafScanline);

    GDALClose(hDataset);

    if(verbose){
        printf("Read raster finished\n");
    }

    return;
}


