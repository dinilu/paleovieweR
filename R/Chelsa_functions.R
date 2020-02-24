#' Title
#'
#' @param chelsa_var
#' @param months
#'
#' @return
#' @export
#'
#' @examples
download_chelsa <- function(chelsa_var, months, directory){
  output <- paste(directory, chelsa_var, sep="/")
  .check_and_create_directory(output)
  # if(!dir.exists(chelsa_var)){ dir.create(chelsa_var) }
  months <- sprintf("%02d", months)
  for(month in months){
    if(chelsa_var %in% c("tmax", "tmin")){
      urldir <- paste0("https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/", chelsa_var, "/")
      file <- paste0("CHELSA_", chelsa_var, "10_", month, "_1979-2013_V1.2_land.tif")
    }else{
      urldir <- paste0("https://www.wsl.ch/lud/chelsa/data/climatologies/", chelsa_var, "/")
      file <- paste0("CHELSA_", chelsa_var, "_", month, "_V1.2_land.tif")
    }
    download.file(paste0(urldir, file), paste0(output, "/", file), "curl")
  }
}

#' Title
#'
#' @param var
#' @param batim
#' @param w
#' @param elev_thrs
#' @param ncell_thrs
#' @param ext
#'
#' @return
#' @export
#'
#' @examples
extrapolate_climate <- function(var, batim, w = c(25, 25),
                                elev_thrs = 20, ncell_thrs = 50,
                                ext = FALSE) {

  lgm_seashore <- batim
  values(lgm_seashore) <- ifelse(values(lgm_seashore) < -120, NA, values(lgm_seashore))

  mask <- var
  values(mask) <- ifelse(is.na(values(mask)), NA, 1)

  mask_mod <- batim
  values(mask_mod) <- ifelse(values(mask_mod) <= elev_thrs, NA, 1)

  dem <- batim
  ll_matrix <- raster::xyFromCell(batim, 1:ncell(batim))
  lat <- batim
  values(lat) <- ll_matrix[, 2]
  lon <- batim
  values(lon) <- ll_matrix[, 1]
  preds <- stack(dem, lat, lon)
  preds <- preds * mask
  names(preds) <- c("dem", "lat", "lon")

  intercept <- preds[[1]]
  intercept[] <- NA
  dem_est <- intercept
  lat_est <- intercept
  lon_est <- intercept

  pb <- txtProgressBar(min = 0, max = nrow(preds), initial = 1, style=3)
  for (rl in 1:nrow(preds)) {
    setTxtProgressBar(pb, rl)
    x <- getValuesFocal(preds * mask_mod, row = rl, nrows = 1, ngb = w, array = FALSE)
    x_int <- rep(NA, nrow(x[[1]]))
    x1 <- rep(NA, nrow(x[[1]]))
    x2 <- rep(NA, nrow(x[[1]]))
    x3 <- rep(NA, nrow(x[[1]]))
    y <- getValuesFocal(var * mask_mod, row = rl, nrows = 1, ngb = w, array = FALSE)

    for (i in 1:nrow(x[[1]])) {
      xy <- na.omit(data.frame(x1 = x[[1]][i, ],
                               x2 = x[[2]][i, ],
                               x3 = x[[3]][i, ],
                               y = y[i, ]))

      if (nrow(xy) > ncell_thrs & nrow(xy) <= 624) {
        # if (nrow(xy) > ncell_thrs) {
        coefs <- coefficients(lm(as.numeric(xy$y) ~ as.numeric(xy$x1) +
                                   as.numeric(xy$x2) + as.numeric(xy$x3)))

        x_int[i] <- coefs[1]
        x1[i] <- coefs[2]
        x2[i] <- coefs[3]
        x3[i] <- coefs[4]
      } else {
        x_int[i] <- NA
        x1[i] <- NA
        x2[i] <- NA
        x3[i] <- NA
      }
    }

    intercept[rl, ] <- x_int
    dem_est[rl, ] <- x1
    lat_est[rl, ] <- x2
    lon_est[rl, ] <- x3

  }
  close(pb)

  coeffs <- stack(intercept, dem_est, lat_est, lon_est,
                  (intercept + dem * dem_est +
                     lat * lat_est + lon * lon_est), var)
  names(coeffs) <- c("intercept", "dem_est", "lat_est", "lon_est", "fitted", "observed")

  if(ext == TRUE){
    var_pred <- coeffs$fitted
    values(var_pred) <- ifelse(is.na(values(var)) & !is.na(values(lgm_seashore)), values(var_pred), values(var))
    return(var_pred)
  }else{
    return(coeffs)
  }
}
