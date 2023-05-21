# === DEFINE FUNCTIONS====

# 1) function for creating names with object type as prefix and population name as suffix
getobjnames.popsufix <- function(objprefix, popabv_var)  { 
  obj_popnames <- c()
  for (popindex in 1:length(Population)) {
    obj_pop <- paste(objprefix, popabv_var[popindex], sep = "")
    obj_popnames <- c(obj_popnames, obj_pop)
  }
  print(obj_popnames)
}

# 2) function for creating names with object type as sufficand population name as prefix
getobjnames.popprefix <- function(objsufix, popabv_var)  { 
  obj_popnames <- c()
  for (popindex in 1:length(Population)) {
    obj_pop <- paste(popabv_var[popindex], objsufix,  sep = "")
    obj_popnames <- c(obj_popnames, obj_pop)
  }
  print(obj_popnames) # getting function to print (outside for loop) allows output to be assigned as a variable
}


# e.g.: popnames_mono <- getobjnames.popprefix(objsufix = "_mono"


# 3) function to subset mapping dataset (df or locus object) into populations
subset.mapping2pop <- function(MappingObj,PopIndexNo){
  MappingObj[MappingObj$population == Population[PopIndexNo], ]
} 

# 4) function for calculating standard error
se <- function(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))) # automatically removes NA
# length(sp_dft_new$Ho) counts length but cant remove NA  
# function for calculating both mean and SE in one go
meanandse <-  function(x) c("mean" = mean(x, na.rm = TRUE), "se" = se(x))


# 5) Edited gl.report.heterozygosity:

gl.report.heterozygosity.edited <- function (gl) 
{
  x <- gl
  sgl <- seppop(x)
  sgl <- sgl[c(3,1,2,4)]  # reorder pops to correct order!
  hs <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, 
                                                                   na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE))) # ADDED na.rm = TRUE as argfor colMeans()
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, 
                                                                    na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE)))
  n <- t(sums/hs)
  n2 <- as.integer(n)
  n <- cbind(row.names(n), n)
  d1 <- data.frame(names = names(hs), hs = as.numeric(hs[1, 
                                                         ]))
  d2 <- data.frame(n)
  d2$X2 <- n2
  names(d2) <- c("names", "freq")
  dd <- left_join(d1, d2, by = "names")
  # dd <- dd[order(dd$hs), ]
  par(mai = c(2.5, 1, 0.5, 0.2))
  barplot(dd$hs, names.arg = paste(dd$names, dd$freq, sep = " | "), 
          las = 2, cex.names = 1, space = 0, border = F, col = rainbow(nrow(dd)), 
          ylab = "Overall mean Observed Heterozygosity",
          ylim = c(0,0.2))
  return(dd)
}


# 6) Edited gl.report.heterozygosity2 - NO PLOT

gl.report.heterozygosity.edited2 <- function (gl) 
{
  x <- gl
  sgl <- seppop(x)
  hs <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, 
                                                                   na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE))) # ADDED na.rm = TRUE as argfor colMeans()
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, 
                                                                    na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE)))
  n <- t(sums/hs)
  n2 <- as.integer(n)
  n <- cbind(row.names(n), n)
  d1 <- data.frame(names = names(hs), hs = as.numeric(hs[1, 
                                                         ]))
  d2 <- data.frame(n)
  d2$X2 <- n2
  names(d2) <- c("names", "freq")
  dd <- left_join(d1, d2, by = "names")
  dd <- dd[order(dd$hs), ]
  
}


# 6) Edited gl.report.heterozygosity2 - NO PLOT

gl.report.heterozygosity.edited2 <- function (gl) 
{
  x <- gl
  sgl <- seppop(x)
  hs <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, 
                                                                   na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE))) 
  
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, 
                                                                    na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE)))
  n <- t(sums/hs)
  n2 <- as.integer(n)
  n <- cbind(row.names(n), n)
  d1 <- data.frame(names = names(hs), hs = as.numeric(hs[1, 
                                                         ]))
  d2 <- data.frame(n)
  d2$X2 <- n2
  names(d2) <- c("names", "freq")
  dd <- left_join(d1, d2, by = "names")
  dd <- dd[order(dd$hs), ]
  
}


# 7) gl.report.monomorphs.edited so it doesnt print long progress bar
# 3/8/22 changed so polymorphic if all heterozygotes
gl.report.monomorphs.edited <- function (gl) {
  x <- gl
  # cat("Identifying monomorphic loci\n")
  a <- vector(mode = "logical", length = nLoc(x))
  for (i in 1:nLoc(x)) {
    a[i] <- NA
  }
  b <- vector(mode = "logical", length = nLoc(x))
  for (i in 1:nLoc(x)) {
    b[i] <- NA
  }
  # c <- vector(mode = "logical", length = nLoc(x))
  # for (i in 1:nLoc(x)) {
  #   c[i] <- NA
  # }
  d <- vector(mode = "logical", length = nLoc(x))
  for (i in 1:nLoc(x)) {
    d[i] <- NA
  }
  
  xmat <- as.matrix(x)
  for (i in (1:nLoc(x))) {
    if (all(is.na(xmat[, i]))) {
      d[i] <- TRUE
      a[i] <- FALSE
      b[i] <- FALSE
     # c[i] <- FALSE
    }
    else {
      a[i] <- all(xmat[, i] == 0, na.rm = TRUE)
      b[i] <- all(xmat[, i] == 2, na.rm = TRUE)
     # c[i] <- all(xmat[, i] == 1, na.rm = TRUE)
      d[i] <- FALSE
    }
    
  }
  s1 <- sum(a, na.rm = TRUE) + sum(b, na.rm = TRUE) # + sum(c, na.rm = TRUE)
  s2 <- s1 - sum(d, na.rm = TRUE)
  polym <- nLoc(x) - s2
  cat("\nBreakdown of", nLoc(x), "loci\n")
  cat("  Polymorphic loci:", polym, "\n  Monomorphic loci:", 
      s1, "\n  Loci with no scores (all NA):", sum(d), 
      "\n")
  return(s2)
}

# Tukey plot with customizable labels
tuk_plot <- function (x, xlab, ylab, ylabels = NULL, ...) {
  for (i in seq_along(x)) {
    xi <- x[[i]][, -4L, drop = FALSE]
    yvals <- nrow(xi):1L
    dev.hold()
    on.exit(dev.flush())
    plot(c(xi[, "lwr"], xi[, "upr"]), rep.int(yvals, 2L), 
         type = "n", axes = FALSE, xlab = "", ylab = "", main = NULL, 
         ...)
    axis(1, ...)
    # change for custom axis labels
    if (is.null(ylabels)) ylabels <- dimnames(xi)[[1L]]
    
    axis(2, at = nrow(xi):1, labels = ylabels, 
         srt = 0, ...)
    abline(h = yvals, lty = 1, lwd = 0.5, col = "lightgray")
    abline(v = 0, lty = 2, lwd = 0.5, ...)
    segments(xi[, "lwr"], yvals, xi[, "upr"], yvals, ...)
    segments(as.vector(xi), rep.int(yvals - 0.1, 3L), as.vector(xi), 
             rep.int(yvals + 0.1, 3L), ...)
    title(main = paste0(format(100 * attr(x, "conf.level"), 
                               digits = 2L), "% family-wise confidence level\n"), 
          # change for custom axis titles
          xlab = xlab, ylab = ylab)
    
    box()
    dev.flush()
    on.exit()
  }
}


filt.maf.pop <- function(gl) {
  S = length(gl$ind.names)
  maf.popcrit = 3 / (2 * S)
  gl.filter.maf(gl, threshold =  maf.popcrit)
}


# old version much slower and works with genelight functions

# function to randomly sample number of  bettongs = to Spurgeon sample size from  gl
# rarefy.pop <- function(gl, use.Ldav = F) {
#   if(use.Ldav){
#     smpl_indnames <- sample(gl$ind.names, min(summary(mappinggl_ldav$pop)))  # randomly sample individuals
#   } else {
#     smpl_indnames <- sample(gl$ind.names, min(summary(mappinggl$pop)))  # randomly sample individuals
#   }
#   gl.keep.ind(gl, smpl_indnames, recalc = FALSE, mono.rm = FALSE)  # aha!! this is wheew
# }
# 
# # function to calculate allelic richness
# CalcAlleleRich <-function(pop_gl) {
#   total       <- pop_gl$n.loc
#   pop_poly    <- gl.filter.monomorphs(pop_gl, v = 0)
#   polymorphic <- pop_poly$n.loc
#   monomorphic <- total - polymorphic
#   ( monomorphic * 1 + polymorphic * 2 ) / total
# }  


# new version much quicker and uses matrix
rarefy.pop2 <- function(gl, use.Ldav = F) {
  if(use.Ldav) {
    mappinggl <- mappinggl_ldav
  }
  minsize <- min(summary(mappinggl$pop))   # find the pop with the smallest sample size 
  smpl_rows <- sample(nInd(gl), minsize)  # randomly sample rows (individuals)
  # subset genelight to sampled rows s youwould matrix
  return(gl[smpl_rows,])
}


CalcAlleleRich2 <-function(gl) {
  monoLoc <- colSums(as.matrix(gl))  %in% c(0, 2*nInd(gl))
  mean(1*monoLoc + 2*!monoLoc)    # one allele for monomorphic loci, two alleles for polymorphic
}  



rarefy.gls.CalcAlleleRich <- function(gl_ls) {             # gl_ls = list of genelights
  rare_gls <- sapply(gl_ls, rarefy.pop2)             # apply function created above returns list of one subset genelight per pop
  rare_gls_AR <- sapply(rare_gls,  CalcAlleleRich2) # then apply report.heterozygosity for each subset gl
  return(rare_gls_AR)
}

rarefy.gls.CalcPrivFix <- function(gl_ls){
  rare_gls <- sapply(gl_ls, rarefy.pop2)
  rare_gl  <- rbind(rare_gls[[1]], rare_gls[[2]], rare_gls[[3]],rare_gls[[4]])  # recombine
  rare_gl_pa <- gl.report.pa(rare_gl)
  return(rare_gl_pa)
}


## function to get 'lookup table' df of indivudal barcodes, index numbers in mappinggl and index numbers within each pop (with popabv prefix e.g.'sp1','sp2', 'dav1', 'dav2' etc.)
get_indrename2nr_df <- function(gl = mappinggl, gls = popsgls, popprfx = popabv) { #  currently restricted to whatever pops are pre-defined in 'popabv'
  
  indrename2nr_popsdf        <- data.frame(matrix(ncol=2,nrow=0))   # first create empty data frame with 2 cols 
  names(indrename2nr_popsdf) <- c("barcode", "popnr")
  for (pop in 1:length(gls)) {
    indnames_old <- indNames(gls[[pop]])
    indnr_pop    <- stringr::str_pad(1:length(indNames(popsgls[[pop]])), 2, pad = "0")  # zero padding
    indnames_new <- paste0(popprfx[pop], indnr_pop)                         # 
    dftemp <- cbind.data.frame("barcode" = indnames_old, "popnr" = indnames_new)
    indrename2nr_popsdf <- rbind.data.frame(indrename2nr_popsdf, dftemp)
  }
  
  indrename2nr_mappingdf <- cbind.data.frame("barcode" = indNames(gl), "nr" = 1:length(indNames(gl)))
  
  indrename2nr_df     <- merge(indrename2nr_popsdf, indrename2nr_mappingdf,       # merge by barcode
                               by = "barcode", all = T)
  return(indrename2nr_df)
  
}


# function to rename inds to shorter ind numbers, either within mapping gl or popsgls (can take list)
rename_barcode2nr <- function(gl, popprfx = popabv, use_popabv = F) {  
  if (class(gl) == "list" ) {                        # i.e. for popgls
    
    for (pop in 1:length(gl)) {
      indnames_old <- indNames(gl[[pop]])
      indnr_pop    <- str_pad(1:length(indNames(popsgls[[pop]])), 2, pad = "0")  # zero padding
      indnames_new <- paste0(popprfx[pop], indnr_pop) 
      indNames(gl[[pop]]) <-  indnames_new                    
    }
    
  } else if (use_popabv == T) {
    temp <- indrename2nr_lookup2[order(indrename2nr_lookup2$nr),]
    indNames(gl) <- temp$popnr
    
  } else {
    indNames(gl) <- 1:length(indNames(gl))
    
  }
  
  return(gl)
}


# something about knitting makes this break - 'do.write_csvs' is not found
write_csv_if <- function(df, outname = "same_as_df", flag = do.write_csvs, folder = "output/", suffix = Vsfx) {
  if (flag) {
    if (outname == "same_as_df") {
      write.csv(df, paste0(folder, deparse(substitute(df)), suffix, ".csv"))
    } else {
      write.csv(df, paste0(folder, outname, "_", suffix, ".csv"))
    }
    
  }
}

# runs and saves first time only, otheriwse loads from .Rdata file to save time
# outname = output object name & name the .Rdata file will be called (string)
# fn = name of function (as string)
# argsls = list of arguments to be parsed to fn, including input data
# folder = Rflder or other folder name, without slashes (string)

# 

run_if_first <- function(outname, fn, argsls, folder = Rflder, overwrite = F) {
  outpath <- paste0(folder, "/", outname, ".Rdata")
  if(!file.exists(file = outpath) | overwrite) {
    assign(outname, do.call(what = fn, args = argsls, quote = FALSE))         # envir = .GlobalEnv (not needed as loaded below)
    save(list = outname, file = outpath) 
  } else {
    message(paste0(outname, " loaded from existing file and was not rerun"))
  }
  load(file = outpath, envir = .GlobalEnv)  
}

# # e.g.
# run_if_first(outname = "wow", fn = "gl.recalc.metrics", argsls = list(mappinggl, "mono.rm" = T))
# run_if_first(outname = "stest", fn = "sapply", argsls = list("X" = popsgls, "FUN" = gl.recalc.metrics,  "mono.rm" = T)) # sapply - you dont have to do anything special just parse the args to sapply

run_all_opts <- function(x, fn, arg_opts) {
  outls <- list()
  for (i in 1:length(arg_opts)) {
    outls[[i]] <- fn(x, arg_opts[i])
  }
  names(outls) <- arg_opts
  return(outls)
}

# e.g. dist_all <- run_all_opts(mappinggl, fn = gl.dist.pop, arg_opts = c("euclidean", "jaccard", "bray"))
# returns list of results with the differnt options



# converts individuals to 'populations so can make use of gl.dist.pop
gl.dist.indpop <- function(gl,  method="euclidean") {
  pop(gl) <- indNames(gl)    # redefine the population information
  indist <- gl.dist.pop(gl, method = method) # calculate euclidean distance matrix
}

