#' Scatter plot with bars at each point.
#'
#' Scatterplot with bars!
#' 
#' @param x, y coordinates
#' @param bar_values matrix with as many rows as x/y coordinates, and as many columns as separate values in the barplot. 
#' @param bar_width_scale overall scale of the bar width
#' @param bar_height_scale overall scale of the bar height
#' @param bar_col vector of colours for the bars, length = ncol(bar_values)
#' @param add Add to existing plot?
#' @param vertical TRUE/FALSE, use vertical or horizontal bars.
#' @return invisible(0), but make a plot
spatial_barplots<-function(x, y, bar_values, bar_width_scale=1, bar_height_scale=1, bar_col = NULL, add=FALSE, vertical=TRUE, ...){

    if(is.null(bar_col)[1]){
        bar_col = rainbow(ncol(bar_values)+1)[1:ncol(bar_values)]
    }
    
    stopifnot(length(x) == length(y))
    stopifnot(nrow(bar_values) == length(x))
    stopifnot(ncol(bar_values) == length(bar_col))

    if(!add){
        plot(x, y, col='white', 
            xlim=range(x) + c(-1,1)*bar_width_scale,
            ylim=range(y) + c(-1,1)*bar_height_scale, ...)
    }

    if(vertical){
        # Vertical bars
        rect(x, y, x+rep(bar_width_scale, nrow(bar_values)), y+bar_values[,1]*bar_height_scale, col=bar_col[1], border=bar_col[1])
        if(ncol(bar_values) > 1){
            for(i in 2:ncol(bar_values)){
                rect(xleft = x, 
                    ybottom = y+rowSums(bar_values[,1:(i-1), drop=FALSE])*bar_height_scale, 
                    xright = x+rep(bar_width_scale, nrow(bar_values)), 
                    ytop = y+rowSums(bar_values[,1:i])*bar_height_scale, 
                    col=bar_col[i], 
                    border=bar_col[i])
            }
        }
    }else{
        # Horizontal bars
        rect(x, y, x+bar_values[,1] * bar_width_scale, y + rep(bar_height_scale, nrow(bar_values)), col=bar_col[1], border=bar_col[1])
        if(ncol(bar_values) > 1){
            for(i in 2:ncol(bar_values)){
                rect(xleft = x+rowSums(bar_values[,1:(i-1), drop=FALSE])*bar_width_scale, 
                    ybottom = y, 
                    xright = x+rowSums(bar_values[,1:i])*bar_width_scale, 
                    ytop = y+rep(bar_height_scale, nrow(bar_values)), 
                    col=bar_col[i], 
                    border=bar_col[i])
            }
        }
    }

    return(invisible(0))
}

.test_spatial_barplots<-function(){

    xs = runif(10)
    ys = runif(10)
    ws = cbind(runif(10), runif(10), runif(10))
    w_scale = 1/rowSums(ws)
    for(i in 1:3) ws[,i] = ws[,i]*w_scale
    bar_col = c('red', 'green', 'blue')
    spatial_barplots(xs, ys, ws, bar_col=bar_col, add=FALSE, bar_width_scale=0.03, bar_height_scale=0.1)
}
