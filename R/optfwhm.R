#' @title Optimal smoothing parameter for All-resolutions inference cluster thresholding in neuroimaging (opt_fwhm)
#' @description The function is made to guide researchers in defining their optimal smoothing parameter
#' @param zstat Z-stats 3D array of activation values, or a nifti file name.
#' @param min_fwhm Minimum level of spatial smoothing a researcher wishes to set, default is 2.
#' @param max_fwhm Maximum level of spatial smoothing a researcher wishes to set, default is 2.
#' @param tdp Minimum TDP threshold set to apply All-Resolutions Inference framework
#' @param max_size Optimization based on maximal size (True) or maximal number of clusters (False)
#' @examples 
#' 
#' optfwhm(zstat, min_fwhm = 2, max_fwhm = 12, tdp = 0.7, max_size = T)
#' 
#' 
#' 
#' 
#' 
#' @return Returns a value that can be used to define the optimal smoothing parameter.
#' @author Jonathan Ornstein, Wouter Weeda.
#'
#' @export


#library(devtools)
#devtools::install_github("wdweeda/ARIBrain")
#library(ARIBrain)

optfwhm = function(zstats, min_fwhm = 2, max_fwhm = 12, tdp, max_size = T){
 
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
