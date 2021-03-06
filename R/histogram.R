`histogram` <- function( y, type="combined", grid="data", breaks=NULL, penalty="default", 
                         greedy=TRUE, right=TRUE, freq=FALSE, control=list(), verbose=TRUE, plot=TRUE,
                         ...) {

  # save y name for later (before doing anything to it)
  
  xname <- deparse( substitute(y))
  
  # check data vector

  if ( any(is.na(y)) ) {
    warning("Removing NAs from data vector")
    y <- y[!is.na(y)]
  }
    
  
  if ( length(unique(y))<2 )
    stop( "data vector must consist of at least two distinct values!" )
    
  # handle invalid penalty/type combination
  penalty = tolower( penalty )
  if ( any( penalty==c("br","nml","sc","mdl")  ) && ( tolower(type)!="regular" && tolower(type)!="r" ) ) {
    warning( "Penalty '", penalty, "' not supported for irregular histograms - creating regular histogram." ) 
    type <- "regular"
  }
  # handle invalid parameter "breaks"
  if ( length(breaks) > 1 ) {
    warning( "Breaks is a vector of length ", length(breaks), " - using first value only", call.=FALSE )
    breaks = breaks[1]
  }
  if ( ! is.null(breaks) ) {
    breaks <- floor( breaks )
    if ( breaks < 2 ) {
      warning( "Breaks must be an integer <= 2 - using breaks=2", call.=FALSE )
      breaks <- 2 
    }
  }
  
  # histogram type: regular
  if ( tolower(type)=="regular" || tolower(type)=="r" )
     out<-histogram.regular( y, penalty=penalty, breaks=breaks, control=control, right=right, freq=freq, verbose=verbose, plot=plot, yvarname=xname, ... )$H

  # histogram type: irregular
  if ( tolower(type)=="irregular" || tolower(type)=="i" )
     out<-histogram.irregular( y, grid=grid, breaks=breaks, penalty=penalty, greedy=greedy, control=control, right=right, verbose=verbose, plot=plot, yvarname=xname, ... )$H

  # histogram type: combined
  if ( tolower(type)=="combined" || tolower(type)=="c" ) {
    
    # check penalty-parameter 
    penalty = tolower( penalty )
    if ( ! any( penalty==c("default","pena","penb","penr")  ) ) {
      warning( "Penalty '", penalty, "' not supported for combined histograms - using default setting for irregular histograms", call.=FALSE )
      penalty = "default"
    }

	  if ( verbose ) 
	    message( "Choosing between regular and irregular histogram:\n\n1.", appendLF=FALSE )
	  out1 <- histogram.regular( y, penalty="br", breaks=NULL, control=control, right=right, verbose=verbose, plot=FALSE )
	  if ( verbose ) 
	    message( "2.",appendLF=FALSE )
    out2 <- histogram.irregular( y, grid=grid, breaks=NULL, penalty=penalty, greedy=greedy, control=control, right=right, verbose=verbose, plot=FALSE )

    #compare maximized likelihood or frgular and irregular histogram
    if (out1$lhvalue>=out2$lhvalue) {
		  out<-out1$H
		  if (verbose) 
		    message("\nRegular histogram chosen.\n")
		}
		else {
      out<-out2$H
      freq <- FALSE
		  if ( verbose ) 
		    message("\nIrregular histogram chosen.\n")
	  }
	         
    # Bugfix: Name of y-var gets lost above - reset it.

    out$xname <- xname
    
	  if ( plot ) 
	    plot(out, freq=freq, ...)
	}


	if ( verbose ) 
	  print( out )	
	  
  return( invisible(out) )
}
