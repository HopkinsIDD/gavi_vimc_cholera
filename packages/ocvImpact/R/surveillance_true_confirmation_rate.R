#' @name surveillance_true_confirmation_rate
#' @title surveillance_true_confirmation_rate
#' @description surveillance_true_confirmation_rate
#' @param country country code
#' @param admin_level admin level
#' @return 
#' @export
#' @include
surveillance_true_confirmation_rate <- function(country, admin_level){
  confirm_rate_fn <- paste0("intermediate_raster/", country, "_trueconfirmrate_admin", admin_level, ".tif")
  confirm_rate_raster <- terra::rast(confirm_rate_fn)

  return(confirm_rate_raster)
}
