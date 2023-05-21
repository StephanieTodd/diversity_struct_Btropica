#' Reports observed and expected heterozygosity by population or by individual. 
#' 
#' Calculates the observed and expected heterozygosities for each population (method="pop") 
#' or the observed heterozyosity for each individual (method="ind") in a genlight object.
#' 
#' Observed heterozygosity for a population takes the proportion of heterozygous
#' loci for each individual then averages over the individuals in that population. 
#' The calculations take into account missing values.
#' 
#' Expected heterozygosity for a population takes the expected proportion of
#' heterozygotes, that is, expected under Hardy-Weinberg equilibrium, for each locus, then
#' averages this across the loci for an average estimate for the population.
#' The calculations of expected heterozygosity use the unbiassed estimates of Nei, M. (1987) 
#' Molecular evolutionary genetics. New York: Columbia University Press.
#'
#' Output for method="pop" is an ordered barchart of observed heterozygosity across 
#' populations together with a table of observed and expected heterozygosity by population.
#' 
#' Observed heterozygosity for individuals is calculated as the proportion of loci that
#' are heterozygous for that individual.
#' 
#' Output for method="ind" is a histogram of heterozygosity across individuals.
#' The histogram is accompanied by a box and whisker plot presented either in standard 
#' (boxplot="standard") or adjusted for skewness (boxplot=adjusted). 
#' 
#' Refer to Tukey (1977, Exploratory Data Analysis. Addison-Wesley) for standard
#' Box and Whisker Plots and Hubert & Vandervieren (2008), An Adjusted Boxplot for Skewed
#' Distributions, Computational Statistics & Data Analysis 52:5186-5201) for adjusted
#' Box and Whisker Plots.
#' 
#' Finally, the loci that are invariant across all individuals in the dataset (that is,
#' across populations), is typically unknown. This can render estimates of heterozygosity
#' analysis specific, and so it is not valid to compare such estimates across species
#' or even across different analyses. This is a similar problem faced by microsatellites.
#' If you have an estimate of the number of invariant sequence tags (loci) in your data,
#' such as provided by gl.report.secondaries, you can specify it with the n.invariant
#' parameter to standardize your estimates of heterozygosity.
#' 
#' @param x -- a genlight object containing the SNP genotypes [Required]
#' @param method -- calculate heterozygosity by population (method='pop') or by individual (method='ind') [default 'pop']
#' @param n.invariant -- an estimate of the number of invariant sequence tags used to adjust the heterozygosity rate [default 0]
#' @param boxplot -- if 'standard', plots a standard box and whisker plot; if 'adjusted',
#' plots a boxplot adjusted for skewed distributions [default 'adjusted']
#' @param range -- specifies the range for delimiting outliers [default = 1.5 interquartile ranges]
#' @param cex.labels -- sets the size of the population labels [default 0.7]
#' @param verbose -- verbosity: 0, silent or fatal errors; 1, begin and end; 2, progress log ; 3, progress and results summary; 5, full report [default 2]
#' @return a dataframe containing population labels, heterozygosities and sample sizes
#' @export
#' @author Bernd Gruber, Arthur Georges and Renee Catullo (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom grDevices rainbow
#' @importFrom graphics par
#' @importFrom plyr join
#' @importFrom robustbase adjbox
#' @examples
#' tmp <- gl.report.heterozygosity(testset.gl,verbose=3)
#' tmp <- gl.report.heterozygosity(testset.gl,method='ind',verbose=3)

# Last amended 28-Jun-19

gl.report.heterozygosity_new <- function(x, 
                                     method="pop", 
                                     n.invariant=0,
                                     boxplot="adjusted",
                                     range=1.5,
                                     cex.labels=0.7, 
                                     verbose=2) {
  
# TIDY UP FILE SPECS

  funname <- match.call()[[1]]

# FLAG SCRIPT START

  if (verbose < 0 | verbose > 5){
    cat("  Warning: Parameter 'verbose' must be an integer between 0 [silent] and 5 [full report], set to 2\n")
    verbose <- 2
  }

  if (verbose > 0) {
    cat("Starting",funname,"\n")
  }

# STANDARD ERROR CHECKING
  
  if(class(x)!="genlight") {
    cat("  Fatal Error: genlight object required!\n"); stop("Execution terminated\n")
  }
  
  if (!(method=="pop" | method == "ind")) {
    cat("Warning: Method must either be by population or by individual, set to method='pop'\n")
    method <- "pop"   
  }
  
  if (n.invariant < 0){
    cat("Warning: Number of invariant loci must be non-negative, set to zero\n")
    n.invariant <- 0
  }
  
  if (!(boxplot=="standard" | boxplot == "adjusted")) {
    cat("Warning: Box-whisker plots must either standard or adjusted for skewness, set to boxplot='adjusted'\n")
    boxplot <- 'adjusted'   
  }


  # Set a population if none is specified (such as if the genlight object has been 
  # generated manually)
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
      if (verbose >= 2){ cat("  Population assignments not detected, 
                             individuals assigned to a single population labelled 'pop1'\n")}
      pop(x) <- array("pop1",dim = nInd(x))
      pop(x) <- as.factor(pop(x))
    }

  # Check for monomorphic loci
    tmp <- gl.filter.monomorphs(x, v=0)
    if ((nLoc(tmp) < nLoc(x)) & verbose >= 2) {
      cat("  Warning: genlight object contains monomorphic loci\n")
    }

# DO THE JOB FOR POPULATIONS
    
  if (method=="pop"){

# Split the genlight object into a list of populations
  sgl <- seppop(x)
  # sgl <- sgl[c("Spurgeon", "Davies Creek", "Emu Creek", "Tinaroo")] # reorder so in correct pop order
  
# OBSERVED HETEROZYGOSITY
  if (verbose >=2){cat("  Calculating Observed Heterozygosities, averaged across loci, for each population\n")}
    
  # Calculate heterozygosity for each population in the list
  # =======================================================================================
  # previously Ho was being calculated only for loci with 100% call rate.
  # lack of na.rm for colMeans meant loci with any missing data produced NAs
  # ======================================================================================
  Ho <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x)==1, na.rm=TRUE), na.rm=TRUE) )) # MOVED na.rm = TRUE to arg. for colMeans() instead of as.matrix()
                                                              
  # Calculate the number of loci used
  # ======================================================================================
  # as far as I cant tell, the original script for 'nl' just calculates the number of loci with 
  # any missing data. If you wanted to calculate the number of loci used above for Ho
  # you would have to use !is.na() instead. However, I'm not sure if this is really what you're after, 
  # given my correction above. Perhaps something like nLoc()?
  # =====================================================================================
  nl <- colSums(data.frame(lapply(sgl, function(x) !is.na(colMeans(as.matrix(x, na.rm=TRUE)))))) # changed is.na() to !is.na()

  nl <- data.frame(lapply(sgl, nLoc)) ## or use this instead?
  
  nl.adj <- nl + n.invariant
  Ho.adj <- Ho*nl/nl.adj
  
  # Calculate sample sizes
  #===================================================================================
  # the sample sizes this function now produces may be confusing for some users (who have missing data) -
  # but i believe they are more appropriate. 'n' is now 'average number of indiviudals with data per locus' 
  # instead of number of physical animals sampled. 
  #==================================================================================
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x)==1, na.rm=TRUE), na.rm=TRUE) )) # MOVED na.rm = TRUE to arg. for colSums() instead of as.matrix()
  n <- t(sums/Ho)
  n2 <- as.integer(n) # convert to 'whole' indiviudals
  nl <- t(nl)
  nl.adj  <- t(nl.adj)
  
  df2 <- cbind.data.frame("pop" = row.names(n), 
                          "nInd" = n2,
                          "nLoc"= nl,
                          "nLoc.adj" = nl.adj)
  
  
  # Join the sample sizes with the heteozygosities
  df1 <- data.frame(pop=names(Ho), 
                    Ho=as.numeric(Ho[1,]), 
                    Ho.adj=as.numeric(Ho.adj))

  df <- plyr::join(df1,df2, by="pop")  

# EXPECTED HETEROZYGOSITY
  if (verbose >=2){cat("  Calculating Expected Heterozygosities, averaged across loci, for each population\n")}
  
  Hexp <- array(NA,length(sgl))
  Hexp.adj <- array(NA,length(sgl))
  # For each population
  for (i in 1:length(sgl)){
    gl <- sgl[[i]]
    gl <- dartR:::utils.recalc.freqhomref(gl,v=0)
    gl <- dartR:::utils.recalc.freqhomsnp(gl,v=0)
    gl <- dartR:::utils.recalc.freqhets(gl,v=0)
    p <- gl@other$loc.metrics$FreqHomRef
    q <- gl@other$loc.metrics$FreqHomSnp
    hets <- gl@other$loc.metrics$FreqHets
    p <- (2*p + hets)/2
    q <- (2*q + hets)/2
    H <- 1 - (p*p + q*q)
    Hexp[i] <- mean(H,na.rm=T)
    Hexp.adj[i] <- Hexp[i]*nl[i]/nl.adj[i]
  }
  
  df <- data.frame(popNames(x),
                   as.numeric(table(pop(x))),
                   nl,
                   round(df$Ho,6),
                   round(df$Ho.adj,6),
                   round(Hexp,6),
                   round(Hexp.adj,6)
                   )
  names(df) <- c("pop","nInd","nLoc","Ho","Ho.adj","He","He.adj")
  
  op <- par(mfrow=c(2,1),mai=c(1.7,0.5,0.1,0),oma=c(2,2,2,0), pty="m")
  df.ordered <- df[order(df$Ho),]
  barplot(df.ordered$Ho, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0.5, border=F, col=rainbow(nrow(df.ordered)), 
          main="Observed Heterozygosity by Population")
  df.ordered <- df[order(df$He),]
  barplot(df.ordered$He, names.arg=paste(df.ordered$pop, df.ordered$nInd, sep=" | "), las=2, cex.names=cex.labels, space=0.5, border=F, col=rainbow(nrow(df.ordered)), 
          main="Expected Heterozygosity by Population")
  
# OUTPUT REPORT
  if (verbose >= 3){
    cat("Reporting Heterozygosity by Population\n")
    cat("No. of loci =", nLoc(x), "\n")
    cat("No. of individuals =", nInd(x), "\n")
    cat("No. of populations =", nPop(x), "\n")
  
    cat("  Miniumum Observed Heterozygosity: ",round(min(df$Ho,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat(" [Corrected:",round(min(df$Ho.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("  Maximum Observed Heterozygosity: ",round(max(df$Ho,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat(" [Corrected:",round(max(df$Ho.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("  Average Observed Heterozygosity: ",round(mean(df$Ho,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat(" [Corrected:",round(mean(df$Ho.adj,na.rm=TRUE),6),"]\n\n")
    } else {
      cat("\n\n")
    }  
    
    cat("  Miniumum Expected Heterozygosity: ",round(min(df$He,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat(" [Corrected:",round(min(df$He.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("  Maximum Expected Heterozygosity: ",round(max(df$He,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat(" [Corrected:",round(max(df$He.adj,na.rm=TRUE),6),"]\n")
    } else {
      cat("\n")
    }  
    cat("  Average Expected Heterozygosity: ",round(mean(df$He,na.rm=TRUE),6))
    if (n.invariant > 0) {
      cat(" [Corrected:",round(mean(df$He.adj,na.rm=TRUE),6),"]\n\n")
    } else {
      cat("\n\n")
    }  
    
    if (n.invariant > 0){
      cat("  Average correction factor for invariant loci =",mean(nl/(nl.adj),rn.na=TRUE),"\n")
    } else {
      cat("  Heterozygosity estimates not corrected for uncalled invariant loci\n")
    }
  
    if (n.invariant > 0 ) {
      print(df)
    } else {
      print(df[,c("pop","nInd","nLoc","Ho","He")])
    }
 
  }
  
  }
  
  # DO THE JOB FOR INDIVIDUALS
  
  if (method=="ind"){
    # Convert to matrix
    m <- as.matrix(x)
    
    # For each individual determine counts of hets, homs and NAs
    c.na <- array(NA, nInd(x))
    c.hets <- array(NA, nInd(x))
    c.hom0 <- array(NA, nInd(x))
    c.hom2 <- array(NA, nInd(x))
    for (i in 1:nInd(x)){
      c.na[i] <- sum(is.na(m[i,]))
      c.hets[i] <- sum(m[i,]==1,na.rm=TRUE)/(nLoc(x)-c.na[i])
      c.hom0[i] <- sum(m[i,]==0,na.rm=TRUE)/(nLoc(x)-c.na[i])
      c.hom2[i] <- sum(m[i,]==2,na.rm=TRUE)/(nLoc(x)-c.na[i])
    }
    
    # Join the sample sizes with the heteozygosities
    df <- cbind.data.frame(x@ind.names, c.hets, c.hom0, c.hom2)
    names(df)<- c("ind.name", "Ho", "f.hom.ref", "f.hom.alt")
    
    # Prepare for plotting
    # Save the prior settings for mfrow, oma, mai and pty, and reassign
      op <- par(mfrow = c(2, 1), oma=c(1,1,1,1), mai=c(0.5,0.5,0.5,0.5),pty="m")
    # Set margins for first plot
    par(mai=c(1,0.5,0.5,0.5))
    # Plot Box-Whisker plot
    if (boxplot == "standard"){
      whisker <- boxplot(df$Ho, horizontal=TRUE, col='red', range=range, main = "Heterozygosity by Individual")
      if (length(whisker$out)==0){
        cat("  Standard boxplot, no adjustment for skewness\n")
      } else {
        outliers <- data.frame(spacer="     ",ID=as.character(df$ind.name[df$Ho %in% whisker$out]),
                          Ho=whisker$out
        )
        cat("  Standard boxplot, no adjustment for skewness\n")
      }
      
    } else {
      whisker <- robustbase::adjbox(data=df, x=as.numeric(df$Ho),
           horizontal = TRUE,
           col='red',
           range=range,
           main = "Heterozygosity by Individual")
      if (length(whisker$out)==0){
        cat("  Boxplot adjusted to account for skewness\n")
      } else {
        outliers <- data.frame(ID=as.character(df$ind.name[df$Ho %in% whisker$out]),
          Ho=whisker$out
          )
        cat("  Boxplot adjusted to account for skewness\n")
      }
    }  
    # Set margins for second plot
    par(mai=c(0.5,0.5,0,0.5))  
    # Plot Histogram
      hist(c.hets, col='red', main=NULL)
      
    # OUTPUT REPORT
      if (verbose >= 3){
        cat("Reporting Heterozygosity by Individual\n")
        cat("No. of loci =", nLoc(x), "\n")
        cat("No. of individuals =", nInd(x), "\n")

        cat("  Miniumum Observed Heterozygosity: ",round(min(df$Ho),6),"\n")
        cat("  Maximum Observed Heterozygosity: ",round(max(df$Ho),6),"\n")
        cat("  Average Observed Heterozygosity: ",round(mean(df$Ho),6),"\n\n")
        cat("  Results returned as a dataframe\n\n")
        if (length(whisker$out)==0){
          cat("  No outliers detected\n")
        } else {  
          cat("  Outliers detected -- \n")
          print(outliers)
        }  
      }   
  }  
    
# FLAG SCRIPT END

  if (verbose > 0) {
    cat("Completed:",funname,"\n")
  }

  # Reset the par options    
    par(op)
  
  # Return the result
  return(df) 
}
