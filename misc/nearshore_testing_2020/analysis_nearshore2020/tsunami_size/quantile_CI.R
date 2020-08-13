# Here is a function to do most of the work
quantile_CI<-function(x, Quant=0.5, level=0.95){
    alpha=(1-level)
    n=length(x) # Size of data
    allQuants=pbinom(0:n, n, Quant) # All quantiles
    lower_index=max(which(allQuants<=alpha/2))
    upper_index=min(which(allQuants>=(1-alpha/2)))
    if(lower_index==-Inf | upper_index==n+1){
        message=paste( ' Cannot compute Quantile ' ,Quant, ' at ' , level*100, ' % confidence ' ,  ' with ' , n, ' data points ' )
        stop(message)
    }
    output=list()
    # Interval
    output$ci=sort(x)[c(lower_index,upper_index)]
    # Exact probability (>level)
    output$p = allQuants[upper_index]-allQuants[lower_index]
    return(output)
}
