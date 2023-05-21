# functions useful across projects
# use source("C:/Users/steph/OneDrive - James Cook University/Desktop (manual backup)/Rfunctions.R")
# now source("C:/Users/steph/Documents/R/Rfunctions.R")
#file.edit("C:/Users/steph/OneDrive - James Cook University/Desktop (manual backup)/Rnotes.txt")

write_csv_if <- function(df, outname = "same_as_df", flag = do.write_csvs, folder = "output/", suffix = Vsfx) {
  if (flag) {
    if (outname == "same_as_df") {
      write.csv(df, paste0(folder, deparse(substitute(df)), suffix, ".csv"))
    } else {
      write.csv(df, paste0(folder, outname, "_", suffix, ".csv"))
    }
    
  }
}

write_csv <- function(df, outname = "same_as_df", flag = T, folder = "output/", suffix = "") {
  # if (is.null(flag) & exists(do.write_csvs, envir = .GlobalEnv)) {   # cant get this to work so just changed defaults above 
  #   flag = do.write_csvs
  # } else if (is.null(flag)) {
  #   flag = T
  # }
  # if (is.null(suffix) & exists(Vsuffix, envir = .GlobalEnv)) {
  #   suffix = Vsuffix
  # } else if (is.null(Vsuffix)) {
  #   suffix = ""
  # }
  if (flag) {
    if (outname == "same_as_df") {
      write.csv(df, paste0(folder, deparse(substitute(df)), suffix, ".csv"))
    } else {
      write.csv(df, paste0(folder, outname, "_", suffix, ".csv"))
    }
    
  }
}


write__sf <- function(sf, suffix, outname = "same_as_df", folder = "output/", ...){
  for (i in 1:length(suffix)) {
    if (outname == "same_as_df") {
      write_sf(sf, paste0(folder, deparse(substitute(sf)), suffix[i]), ...)
    } else {
      write_sf(sf, paste0(folder, outname, suffix[i]), ...)
    }
  }

}



save_Rdata <- function(x, folder = Rflder, filename = NULL, suffix = NULL,  withdate = F) {
  if(is.null(filename)) {
    filename = deparse(substitute(x))
  }
  if (withdate) {
    save(x, file = paste0(folder, "/", filename, suffix, Sys.Date(), ".Rdata"))
  } else {
    save(x, file = paste0(folder, "/", filename, suffix, ".Rdata"))
  }
}


# runs and saves first time only, otheriwse loads from .Rdata file to save time
# outname = output object name & name the .Rdata file will be called (string)
# fn = name of function (as string)
# argsls = list of arguments to be parsed to fn, including input data
# folder = Rflder or other folder name, without slashes (string)


# 
# run_if_first <- function(outname, fn, argsls, folder = Rflder) {
#   outpath <- paste0(folder, "/", outname, ".Rdata")
#   if(!file.exists(file = outpath)) {
#     assign(outname, do.call(what = fn, args = argsls, quote = FALSE))         # envir = .GlobalEnv (not needed as loaded below)
#     save(list = outname, file = outpath) 
#   }
#   load(file = outpath, envir = .GlobalEnv)  
# }

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



transpose_df <- function(df,  label, varible_names) {
  dft <- as.data.frame(t(df))
  names(dft) <- varible_names
  dft <- cbind.data.frame(label = names(df), dft)
  row.names(dft) <- NULL
  return(dft)
}




transpose_df <- function(df,  label, variable_names) {
  dft <- as.data.frame(t(df))
  names(dft) <- variable_names
  dft <- cbind.data.frame(names(df), dft)
  names(dft)[1] <- label
  row.names(dft) <- NULL
  return(dft)
}



replace_nulls <- function(x, replacement = 0) {
  x[is.na(x)] <- replacement
  return(x)
}

# if replacing with more than type of value (i.e. '0' and 'FALSE'), replacement needs to be list of same length as vars
replace_nulls_multi <- function(df, vars, replacement = NULL) {
  if (is.null(replacement)) {
    replacement = 0
  }
  if (length(replacement) == 1) {
    df[,vars][is.na(df[,vars])] <- replacement
  } else {
    for (i in 1:length(replacement)) {
      df[,vars[i]][is.na(df[,vars[i]])] <- replacement[[i]]
    }
  }
  
  return(df)
}


check_twoway_mismatch <- function(x,y, by = c("names", "all")){  # returns mismatches

  if (by == "names") {
    missing_from_y <- names(x)[!names(x) %in% names(y)]
    missing_from_x <- names(y)[!names(y) %in% names(x)]
  } else {
    missing_from_y <- x[!x %in% y]
    missing_from_x <- y[!y %in% x]
  }
  message("\n in ", deparse(substitute(x)), " but missing from ", deparse(substitute(y)))
  print(unique(missing_from_y))
  message("\n in ", deparse(substitute(y)), " but missing from ", deparse(substitute(x)))
  print(unique(missing_from_x))
}


check_matching <- function(x,y, by = c("names", "all")){  # returns mismatches
  if (by == "names") {
    matching <- names(x)[names(x) %in% names(y)]
    message(" \n column names that match: ")
  } else {
    matching <- x[x %in% y]
    message(" \n values that match: ")
  }
  print(matching)
}

# check missing values for column "X" in a dataframe. 
# options to return df: "none", "missing"  or "not_missing" (default = "none")

check_missing <- function(df, x, returndf = c("none", "missing", "not_missing")) {
  if (length(returndf) >1) {
    returndf = "none"
  }
  
  for (i in 1:length(x)) {
    if (sum(is.na(df[,x[i]])) > 0) {
      message("\n There are ", sum(is.na(df[,x[i]])), " missing values of ", x[i], " in ", deparse(substitute(df)))
      missingdf <- df[is.na(df[,x[i]]),]
      
      if (returndf == "missing") {
        return(missingdf)
      } 
      
    } else  {
      message("\n There are no missing values of ", x[i], " in ", deparse(substitute(df)))
    }
  }
  
  if (returndf == "not_missing") {
    not_missingdf <- df[!(is.na(df[,x])),]
    message("\n The ", sum(is.na(df[,x])), " rows with missing values of ", x[i], " have now been removed from  ", deparse(substitute(df)))  # assumes output is assigned same name as input
    return(not_missingdf)
  }
}
  


# check_duplicates <- function(df, x, returndf = c("none", "dup_summary", "distinct")) {   # options: "dup_summary", "distinct", DEFAULT = "none"
#   if (length(returndf) >1) {
#     returndf = "none"
#   }
#   if (length(unique(df[,x])) != nrow(df)) {
#     message(" \n There are ", nrow(df) - length(unique(df[,x])), " duplicate values of ", x, " in ", deparse(substitute(df)))
#     
#     if (returndf == "dup_summary") {
#       
#       n_occur <- data.frame(table(df[,x]))
#       duplicatesdf <- n_occur[n_occur$Freq >1 ,]
#       names(duplicatesdf)[1] <- x
#       
#       out <- merge(duplicatesdf, df)  # subset to only duplicated surveys
#       return(out)
#     }
#     
#     if (returndf == "distinct") {
#       distinct_df <- distinct(df, x = (df[,x]), .keep_all = T) # subset to only duplicated surveys
#       message("\n The ",  nrow(df) - nrow(distinct_df), " rows with duplicate values of ", x, " have now been removed from  ", deparse(substitute(df)))  # assumes output is assigned same name as input
#     }
#     
#   } else {
#     message("there are no duplciates of ", x, " in ", deparse(substitute(df)))
#     distinct_df <- df
#   }
#   
#   if (returndf == "distinct") {
#     return(distinct_df)
#   }
#   
# }

# df <- wndat
# x = "SIGHTING"

check_duplicates <- function(df, x, returndf = c("none", "dup_summary", "distinct")) {   # options: "dup_summary", "distinct", DEFAULT = "none"
  if (length(returndf) >1) {
    returndf = "none"           # default
  }
  
  count_unique = nrow(data.frame(unique(df[,x]))) # data.frame allows x to be >1
  
  if (count_unique != nrow(df)) {
    message(" \n There are ", (nrow(df) - count_unique), " duplicate values of ", x, " in ", deparse(substitute(df)))
    
    if (returndf == "dup_summary") {
      
      n_occur <- data.frame(table(df[,x]))
      duplicatesdf <- n_occur[n_occur$Freq >1 ,]
      names(duplicatesdf)[1:length(x)] <- x
      
      out <- merge(df, duplicatesdf)  # subset to only duplicated surveys
      return(out)
    }
    
    if (returndf == "distinct") {
      distinct_df <- distinct(df, x = (df[,x]), .keep_all = T) # subset to only duplicated surveys
      message("\n The ",  nrow(df) - nrow(distinct_df), " rows with duplicate values of ", x, " have now been removed from  ", deparse(substitute(df)))  # assumes output is assigned same name as input
    }
    
  } else {
    message("there are no duplciates of ", x, " in ", deparse(substitute(df)))
    distinct_df <- df
  }
  
  if (returndf == "distinct") {
    return(distinct_df)
  }
  
}


check_merge <- function(df, key) {
  check_missing(df, key)
  check_duplicates(df, key, returndf = "none")
}



remove_empty_cols <- function(df) {
  empty_columns <- sapply(df, function(x) all(is.na(x) | x == ""))
  df2 <- df[, !empty_columns]
  return(df2)
}


# x needs to be ordered
# origin EHC data?
aggregate_relative_range <- function(x, range) {
  # gr <-  vector(length = length(x))
  gr <- c()
  gr[1] <- 1
  for (i in 2:length(x)) {
    gr[i] <- ifelse(x[i] <= x[i-1] + range, gr[i-1], gr[i-1] + 1)
  }
  return(gr)
}

# created for: SOKBA report conservation vols/hrs
sum_unequal_rows <- function(df, index_cols, sum_cols){   # , test_cols ??
  totals1 <- df %>%
    group_by_at(c(index_cols, all_of(sum_cols))) %>%
    summarise(rep_count = n()) %>%
    ungroup() %>%
    group_by_at(index_cols) %>%
    summarise(across(.cols = all_of(sum_cols), .fns = ~ sum(.x, na.rm = TRUE),.names = "TOTAL_{.col}"),
              unequal_rows_summed = n(),    # added 4/3/22
              total_rows =sum(rep_count))
  return(totals1)
}

sum_if_unequal <- function(a,b) {
    out <- if_else(a != b, a + b, a, missing = ifelse(!is.na(a), a, b))
    return(out)
}


  
# rename files with suffix e.g. datcode
# folder needs to be string with backslash e.g. "SambaR_output_6213_110_HWE/"
append_suffix_to_filenames <- function(suffix, folder, recursive = T){
  require(stringr)
  x          <- list.files(folder, recursive = T)
  extns      <- str_sub(x, -4, -1)
  x_no_extns <- str_sub(x, 1, -5)
  file.rename(from = paste0(folder, x), to = paste0(folder, x_no_extns, suffix, extns))
}

# 'implct2explct' - make implicit zeros explicit by mapping records onto a backbone of possibilities
# created for: GBR 'Wildnet.R'  - dealing with SubCays
# there are optionally one or two 'backbones'
# examples:
#     bb1 = list of all SubCays,  bb2 = Surveys at islands where constituent Subcays could theoretically have been surveyed, returns = identifies if both subcays were surveyed 
#     bb1 = list of all SubCays, bb2 = gbr_trig, returns = implicit zeros for potential Cays*Triggers
implct2explct <- function(dat, bb1, bb2 = NULL, explicit_val = NULL) {
  dat1 <- mutate(dat, explicit = T)  # records in data are by definition 'explicit'
  if (is.null(bb2)) {
    bb = bb1
  } else {
    bb   <- merge(bb1, bb2)  # backbone of possibilities  to map records to 
  }
 
  explicitdf <- merge(bb, dat1, all = T)  # now map actual surveys onto backbone
  
  # fill in implicit zeros (with FALSE) - implicit to explicit missing surveys
  if (!is.null(explicit_val)) {
    explicitdf <- replace_nulls_multi(explicitdf, vars = c('explicit', explicit_val) , replacement = list(F, 0))
  } else {
    explicitdf <- replace_nulls_multi(explicitdf, vars = 'explicit' , replacement = F)
  }
}


# 'export_datls2csv' takes a list of dataframes (e.g. results produced from a for loop) and exports each to .csv (and optionally to glob env.) using names provided in a lookup table
# if argument for ID is provided it will only export one df, otherwise it will export all
# assumes datls is in correct order of ID
# created for: GBR - exporting trigger data used in GLMs
export_datls2csv <- function(ID = NULL, datls, id_name_lookup, namecol, to_env = T) {
  # define function to do it once	
  runonce <- function(ID){        
    labnames <- gsub(" |-", "", id_name_lookup[,namecol])
    dat_ <- datls[[ID]]
    outname <- paste0(labnames[ID], ID)
    write_csv(dat_,  suffix = outname, folder = "output/intermediate/")
    if(to_env) {
      assign(paste0("dat_", outname), dat_, envir = .GlobalEnv)
    }
  }
    # run
    if(!is.null(ID)) 
    {
      runonce(ID)
    } 
    else 
    {
    for (i in 1:length(datls)) {
      ID = i
      runonce(ID)
    }
  }
}



# 'filter_objs_by_class' searches the global environment for objects of a specified class (or classes)
#  filtering can be inclusive or exclusive - to include everything except the specified classes, set exclude = T
#  option to specify a string pattern to filter by object name
#  Value: returns a character vector of object names matching the filter criteria
#  created for SAMBAR in the context of saving objects to file using SOAR pkg
#  useful with function (get)
filter_objs_by_class <- function(class2filt, exclude = F, namepattern = ".*") {  # woudl lie to add pattern = NULL
  out <- c()
  for (i in 1:length(class2filt)) {
    if (exclude) {
      objs <- Filter( function(x) ! class2filt[i] %in% class( get(x) ), ls(.GlobalEnv, pattern = namepattern) )  # %in% doesn't work other way around
    } else {
      objs <- Filter( function(x) class2filt[i] %in% class( get(x) ), ls(.GlobalEnv, pattern = namepattern) )  # ls(.GlobalEnv) necessary otherwise lists objs within function
    }
    out <- c(out, objs)
  }
  return(out)
}


# spatial functions
st_write_gpkg <- function(x) {
  st_write(x, paste0("output/", deparse(substitute(x)), ".gpkg"), append=FALSE)
} 

sf2df <- function(sf){
  if("geometry"  %in% names(sf)){
    sf_df <- sf %>% # save additional cols as df
      as.data.frame() %>% 
      select(-geometry)
  } else {
    sf_df <- sf %>% # save additional cols as df
      as.data.frame() %>% 
      select(-geom)
  }
  return(sf_df)
}

add_area_col <- function(sf){
  temp_area <- st_area(sf)
  out <- sf %>% 
    ungroup %>% 
    mutate(area_ha = as.numeric(temp_area)/10000)
  return(out)
}
# 
# # original
# thin_pts_by_min_dist <- function(dat, mindist) {
#   require("sf")
#   require("dplyr")
#   # checks
#   if (st_crs(dat, parameters = TRUE)$units_gdal != "metre") {
#     dat  %<>% st_transform(crs = 7845) # project
#     message("warning: points have been projected, assuming min dist units are meters")
#   }
#   seed_ptx = 1 # ptx = index of points
#   forbid_ptx <- c()
#   for (i in 1:nrow(dat)) {
#     if (i < nrow(dat)) {
#       if (!i %in% forbid_ptx) {
#         rmn_ptx <- (i+1):nrow(dat)
#         seed_ptx  <- c(seed_ptx, i)
#         seed_buf  <- st_buffer(dat[i,], mindist)
#         seed_test <- lengths(st_intersects(dat[rmn_ptx,], seed_buf))  > 0
#         forbid_ptx <- c(forbid_ptx, rmn_ptx[seed_test])
#         # ptx        <- ptx[!seed_test]
#       }
#     } else if (i == nrow(dat)) {  # special case for last point
#       if (!i %in% forbid_ptx) {
#         seed_ptx <- c(seed_ptx, i)
#       }
#     }
#   }
#   thindat <- dat[seed_ptx,]
#   return(thindat)
# }
# 
# # origignal modified
# thin_pts_by_min_dist_OGmod <- function(dat, mindist) {
#   require("sf")
#   require("dplyr")
#   # checks
#   if (st_crs(dat, parameters = TRUE)$units_gdal != "metre") {
#     dat  %<>% st_transform(crs = 7845) # project
#     message("warning: points have been projected, assuming min dist units are meters")
#   }
#   seed_ptx <- 1 # ptx = index of points
#   forbid_ptx <- c()
#   for (i in 1:nrow(dat)) {
#     if (!i %in% forbid_ptx) {
#       seed_ptx  <- c(seed_ptx, i)
#       if (i < nrow(dat)) {
#         rmn_ptx <- (i+1):nrow(dat)
#         seed_buf  <- st_buffer(dat[i,], mindist)
#         seed_test <- lengths(st_intersects(dat[rmn_ptx,], seed_buf))  > 0
#         forbid_ptx <- c(forbid_ptx, rmn_ptx[seed_test])
#       }
#     }
#   }
#   thindat <- dat[seed_ptx,]
#   return(thindat)
# }

# 
# ## opt 1 - use index of rand_pts
# thin_pts_by_min_dist1 <- function(dat, mindist) {
#   require("sf")
#   require("dplyr")
#   # checks
#   if (st_crs(dat, parameters = TRUE)$units_gdal != "metre") {
#     dat  %<>% st_transform(crs = 7845) # project
#     message("warning: points have been projected, assuming min dist units are meters")
#   }
#   
#   forbid_ptx <- c()
#   set.seed(872436)           # Set seed
#   ptx_rand <- sample(1:nrow(dat))
#   seed_ptx = ptx_rand[1] # ptx = index of points
#   last_pt = ptx_rand[length(ptx_rand)]
#   for (i in 1:length(ptx_rand)) {
#     cur_pt = ptx_rand[i]
#     if (cur_pt != last_pt) {
#       if (!cur_pt %in% forbid_ptx) {
#         rmn_ptx <- ptx_rand[(i+1):length(ptx_rand)]
#         seed_ptx  <- c(seed_ptx, cur_pt)
#         seed_buf  <- st_buffer(dat[cur_pt,], mindist)
#         seed_test <- lengths(st_intersects(dat[rmn_ptx,], seed_buf))  > 0
#         forbid_ptx <- c(forbid_ptx, rmn_ptx[seed_test])
#         # ptx        <- ptx[!seed_test]
#       }
#     } else if (cur_pt == last_pt) {  # special case for last point
#       if (!cur_pt %in% forbid_ptx) {
#         seed_ptx <- c(seed_ptx, cur_pt)
#       }
#     }
#   }
#   thindat <- dat[seed_ptx,]
#   return(thindat)
# }
# 
# test1 <- thin_pts_by_min_dist1(dat, 500)


## opt 2 - loop through only remaining non-forbidden points - random
thin_pts_by_min_dist2 <- function(dat, mindist) {
  require("sf")
  require("dplyr")
  # checks
  if (st_crs(dat, parameters = TRUE)$units_gdal != "metre") {
    dat  %<>% st_transform(crs = 7845) # project
    message("warning: points have been projected into EPSG 7845, assuming min dist units are meters")
  }
  
  forbid_ptx <- c()
  seed_ptx   <- c() # ptx = index of points
  # set.seed(872436)           # Set seed
  ptx_rand <- sample(1:nrow(dat))
  last_pt = ptx_rand[length(ptx_rand)]
  for (i in ptx_rand) {  # loops through ptx_rand directly insted of an index of an index
    if (!i %in% forbid_ptx){
      seed_ptx <- c(seed_ptx, i)
      if (i != last_pt) {
        rmn_ptx   <- ptx_rand[!ptx_rand %in% c(forbid_ptx, seed_ptx)] # subsets ptx_rand by removing seed and forbidden
        seed_buf  <- st_buffer(dat[i,], mindist)
        seed_test <- lengths(st_intersects(dat[rmn_ptx,], seed_buf))  > 0
        forbid_ptx <- c(forbid_ptx, rmn_ptx[seed_test])
      } 
    }
  }
  thindat <- dat[seed_ptx,]
  return(thindat)
}

# 
# ## opt 2 - loop through only remaining non-forbidden points 
# thin_pts_by_min_dist2_nonrand <- function(dat, mindist) {
#   require("sf")
#   require("dplyr")
#   # checks
#   if (st_crs(dat, parameters = TRUE)$units_gdal != "metre") {
#     dat  %<>% st_transform(crs = 7845) # project
#     message("warning: points have been projected into EPSG 7845, assuming min dist units are meters")
#   }
#   
#   forbid_ptx <- c()
#   seed_ptx   <- c() # ptx = index of points
#   # set.seed(872436)           # Set seed
#  # ptx_rand <- sample(1:nrow(dat))
#   ptx_rand <- 1:nrow(dat)
#   last_pt = ptx_rand[length(ptx_rand)]
#   for (i in ptx_rand) {  # loops through ptx_rand directly insted of an index of an index
#     if (!i %in% forbid_ptx){
#       seed_ptx <- c(seed_ptx, i)
#       if (i != last_pt) {
#         rmn_ptx   <- ptx_rand[!ptx_rand %in% c(forbid_ptx, seed_ptx)] # subsets ptx_rand by removing seed and forbidden
#         seed_buf  <- st_buffer(dat[i,], mindist)
#         seed_test <- lengths(st_intersects(dat[rmn_ptx,], seed_buf))  > 0
#         forbid_ptx <- c(forbid_ptx, rmn_ptx[seed_test])
#       } 
#     }
#   }
#   thindat <- dat[seed_ptx,]
#   return(thindat)
# }
# 

# from SOKBA 
# builds on addNewData function to remap (i.e. fix up) values. creates new col and option to fill NAs with unchanged original  
remap <- function(df, old_var_names, new_var_names, suffix, fill_unchanged = TRUE, fill_var = old_var_names) {
  df <- addNewData(paste0("data/remap_", suffix, ".csv"), df, allowedVars = new_var_names)
  if (fill_unchanged) { 
    for (i in 1:length(new_var_names)) {                                                              # incase multiple variables remapped at once
      df[is.na(df[,new_var_names[i]]), new_var_names[i]] <- df[is.na(df[,new_var_names[i]]), fill_var[i]]
    }
    
  }
  return(df)

}



remap2 <- function(df, new_var_names, suffix, fill_unchanged = TRUE, fill_var = NULL) {
  df     <- as.data.frame(df)  # tbl_df causes error
  outdat <- addNewData(paste0("data/remap_", suffix, ".csv"), df, allowedVars = new_var_names)
  if (fill_unchanged) {
    for (i in 1:length(new_var_names)) {                                                              # incase multiple variables remapped at once
      outdat[is.na(outdat[,new_var_names[i]]), new_var_names[i]] <- outdat[is.na(outdat[,new_var_names[i]]), fill_var[i]]
    }
    
  }
  return(outdat)
}
# df <- wwfraw_df_merged
# old_var_names = "maturity2"
# new_var_names = "maturity2"
# suffix = "field_values"
# 
# newDataFileName = paste0("data/remap_", suffix, ".csv")
# data = df
# allowedVars = new_var_names
# i = 1
# addNewData <- function(newDataFileName, data, allowedVars){
#   
#   import <- readNewData(newDataFileName, allowedVars)
#   
#   if( !is.null(import)){    
#     for(i in seq_len(nrow(import))){  #Make replacements
#       col.to <- import$newVariable[i] 
#       col.from <- import$lookupVariable[i]
#       if(is.na(col.from)){ # apply to whole column
#         data[col.to] <- import$newValue[i]
#       } else { # apply to subset
#         rows <- data[[col.from]] == import$lookupValue[i]
#         rows <- ifelse(is.na(rows), FALSE, rows)
#         data[rows,col.to] <- import$newValue[i]
#       }
#     }   
#   }      
#   data
# }
# test <- data[col.from][484]
# 


# function for matching headers - from cage trap data wrangling
check_headers_match <- function(data1, data2) {
  # check for data1 in data2: FALSE is extra cols in data1 not found in data2. result is index of data1
  not_found_in_data2 <- !(colnames(data1) %in% colnames(data2)) 
  not_found_in_data1 <- !(colnames(data2) %in% colnames(data1))
  
  if (sum(not_found_in_data2) == 0) {
    message("all data1 columns found in data2")
  } else { 
    message("There are data1 columns missing or different from data2 \n") 
    print(colnames(data1)[not_found_in_data2])
  } 
  
  if (sum(not_found_in_data1) == 0) {
    message("all data2 columns found in data1")
  } else { 
    message("There are data2 columns missing or different from data1 \n") 
    print(colnames(data2)[not_found_in_data1])
  } 
  
  
  return(list("data1_names_not_found_in_data2" = colnames(data1)[not_found_in_data2], 
              "data1_index_not_found_in_data2" = which(not_found_in_data2),
              "data2_names_not_found_in_data1" = colnames(data2)[not_found_in_data1],
              "data2_index_not_found_in_data1" = which(not_found_in_data1)))
}


# get list of all unique, cleaned column names across a list of dataframes - from cage trap data wrangling
get_all_colnames <- function(dflist) {
  allnames <- data.frame()
  for (i in 1:length(dflist)) {
    dflist[[i]] <- as.data.frame(dflist[[i]])
    names(dflist[[i]])  %<>%  
      clean_strings(.) %>%    # make lower case, remove special characters etc
      gsub(" ", "_", .)   # make snake case
    curnames <- cbind.data.frame("colname" = names(dflist[[i]]), "example_datname" = names(dflist[i]))
    allnames <- rbind(allnames, curnames)
  }
  return(distinct(allnames, colname, .keep_all = T))
}

# replace column headers with chosen synonym, based on a lookup table called 'remap_field_names.csv' - from cage trap data wrangling
replace_header_synonyms <- function(dflist) {
  if (class(dflist) == "list") {
    for (i in 1:length(dflist)) {
      namesdf <- data.frame("oldname" = names(dflist[[i]]), "newname" = NA)
      namesdf %<>%  
        remap(., old_var_names = "oldname",
              new_var_names = "newname", 
              suffix = "field_names")
      names(dflist[[i]]) <- namesdf$newname
    }
    return(dflist)
    
  } else {
    namesdf <- data.frame("oldname" = names(dflist), "newname" = NA)
    namesdf %<>%  
      remap(., old_var_names = "oldname",
            new_var_names = "newname", 
            suffix = "field_names")
    names(dflist) <- namesdf$newname
  }
  
}

# function that keeps non-null rows within each group of a dataframe (different row kept for different fields)
# ouput if one row per group of the most complete data
# id is a column that should be unique in output, that identifies duplicates from different datsets
# - from cage trap data wrangling

# subdat <- awc_df_merged2 %>%
#   filter(capture_id == "TDC1_2020-10-22")
# subdat <- as.data.frame(subdat)
# j = which(names(subdat) == "date")
# i = 28 
# dat = awc_df_merged2
# concat_col = "filename"
# id = "capture_id"
# unique_ids[i]
# # use function that keeps non-null rows within each group (different row kept for different fields)
# awc_df_merged_distinct <- pool_across_duplicates(awc_df_merged2, "capture_id", "filename")

pool_across_duplicates <- function(dat, id, concat_col = NULL) {
  require(dttr2)
  dat <- as.data.frame(dat)
  unique_ids <- unique(dat[,id])
  newdat <- data.frame()
  for (i in 1:length(unique_ids)) {
    # first filter to current group
    subdat <- dat[dat[,id] == unique_ids[i],]
    # within each group:
    if (nrow(subdat) > 1) {
      for (j in 1:ncol(dat)) {
        val2keep <- unique(subdat[,j][!is.na(subdat[,j])]) # keep non-null values
        if (class(val2keep) == "character") {
          val2keep <- unique(val2keep[val2keep != ""]) # keep non-blank values
        }
        
        if (length(val2keep) == 0) {
          is.na(val2keep) <- 1 # sets element 1 to NA of whatever class val2keep is
        }
        
        if (length(val2keep) > 1 & !is.null(concat_col)) {
          val2keep = ifelse(names(dat)[j] == concat_col, paste(val2keep, collapse = ' and '), val2keep[1])
        } else {
          val2keep <- val2keep[1] # doesnt matter if length(val2keep) already == 1
        }
        
        if (class(dat[,j]) == "Date" & class(val2keep) == "numeric") { # fixes 'Error in as.Date.numeric(value) : 'origin' must be supplied'
          val2keep <- as.Date(val2keep, origin = "1970-01-01")
        }
        newdat[i,j]   <- val2keep  # each group now only has one row
        
      } 
    } else {
      newdat[i,1:ncol(subdat)]   <- subdat # or only one row originally present
    }
  }
  
  names(newdat) <- names(dat)
  return(newdat)
  
}

# https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


# convert stupd excel decimal time format to readable character  
dec_hrs2hm <- function(x) {
  paste(floor(x), round((x-floor(x))*60), sep=":")
} 
