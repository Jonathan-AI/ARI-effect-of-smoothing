# Importing the packages


library(aws)
library(ARIbrain)
library(scatterplot3d)
library(TCIU)
library(RNifti)
library(ggplot2)
library(ggthemes)

# Function

## Version that is optimal in terms of largest area found by clusters 

opt_fwhm = function(zstats, min_fwhm = 2, max_fwhm = 12, tdp, max_size = T){
  
  ### Calculating voxelsize??????
  cluster_calc = function(fwhm, tdp, zstats){
    
    sm_arr = attr(kernsm(zstats, h = fwhm, unit = "FWHM")[], "yhat")
    new_mat = sm_arr
    mask_mat = ifelse(zstats == 0, F, T)
    p_mat = (1 - pnorm(abs(new_mat)))*2
    
    fx = function(p_mat, mask_mat, min_tdp){
      
      tdpsmcl = TDPQuery(ARIBrainCluster(Pmap = p_mat, mask = mask_mat), tdp)
      sm = as.data.frame(summary(tdpsmcl))
      return(sum(sm$Size) / length(p_mat[mask_mat]))
    }
    
    suppressErrorsWithTry <- function(f, p_mat, mask_mat, min_tdp) {
  result <- try(f(p_mat, mask_mat, min_tdp), silent = TRUE)
  if (inherits(result, "try-error")) {
    return(0)
  } else {
    return(result)
  }
  }
  
  suppressErrorsWithTry(fx, p_mat, mask_mat, min_tdp)
  }
  
  size_calc = function(fwhm, tdp, zstats){
    
    sm_arr = attr(kernsm(zstats, h = fwhm, unit = "FWHM")[], "yhat")
    new_mat = sm_arr
    mask_mat = ifelse(zstats == 0, F, T)
    p_mat = (1 - pnorm(abs(new_mat)))*2

    fx = function(p_mat, mask_mat, min_tdp){
      
      tdpsmcl = TDPQuery(ARIBrainCluster(Pmap = p_mat, mask = mask_mat), tdp)
      sm = as.data.frame(summary(tdpsmcl))
      return(nrow(sm))
    }
    
    suppressErrorsWithTry <- function(f, p_mat, mask_mat, min_tdp) {
  result <- try(f(p_mat, mask_mat, min_tdp), silent = TRUE)
  if (inherits(result, "try-error")) {
    return(0)
  } else {
    return(result)
  }
  }
  
  suppressErrorsWithTry(fx, p_mat, mask_mat, min_tdp)
  }
    
  
  
 # tdp_set = tdp
  #zstats_set = zstats
  fwhm_options = c(min_fwhm:max_fwhm)
  if (max_size == T){
  return(fwhm_options[which.max(sapply(fwhm_options, function(fwhm) cluster_calc(fwhm, 0.7, zstats)))])
  }
  else{
    return(fwhm_options[which.max(sapply(fwhm_options, function(fwhm) size_calc(fwhm, 0.7, zstats)))])
  }
}
