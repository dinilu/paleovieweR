.check_existing_file <- function(pathfile, overwrite) {
  if(file.exists(pathfile)){
    if(overwrite == FALSE){
      stop("File already exists and overwrite option set to FALSE, if you want to overwrite turn the option to TRUE.")
    }
  }
}

.check_and_create_directory <- function(out.path) {
  if(!dir.exists(out.path)){
    dir.create(out.path)
  }
}

.open_or_create_ncfile <- function(nc.source, vars, suppress_dimvals = TRUE) {
  if(file.exists(nc.source)){
    cat("Opening existing NetCDF file.", "\n")
    nc.trg <- ncdf4::nc_open(nc.source, write = TRUE, suppress_dimvals = supress_dimvals)
  }else{
    cat("Creating new NetCDF file.", "\n")
    nc.trg <- ncdf4::nc_create(nc.source, vars)
  }
  return(nc.trg)
}

.flush_ncfile <- function(nc, suppress_dimvals = TRUE){
  nc.file <- nc$filename
  ncdf4::nc_close(nc)
  nc <- ncdf4::nc_open(nc.file, write=T, suppress_dimvals = supress_dimvals)
  return(nc)
}

#' Crop PaleoView NetCDF files
#'
#' @param file
#' @param ext
#' @param out.path
#' @param overwrite
#'
#' @return
#' @export
#'
#' @examples
crop_paleoview <- function(nc.source, ext, out.path, overwrite = FALSE, suppress_dimvals = TRUE){
  # nc.source <- "precipitation-5000BP-1989AD.nc"
  # ext <- study_area
  # out.path <- "croped"
  # overwrite <- FALSE
 
  in.filename <- basename(nc.source)
  
  outpathfile <- paste(out.path, in.filename, sep="/")

  .check_existing_file(outpathfile, overwrite)
  
  nc <- ncdf4::nc_open(nc.source, suppress_dimvals = supress_dimvals)
  lon <- ncdf4::ncvar_get(nc, "longitudes")
  lon_i <- which(lon >= ext[1] & lon <= ext[2])
  lat <- ncdf4::ncvar_get(nc,"latitudes")
  lat_i <- which(lat >= ext[3] & lat < ext[4])
  mon <- ncdf4::ncvar_get(nc, "months")
  mon_i <- which(mon >= 1 & mon <= 12)

  londim <- ncdf4::ncdim_def("lon", "degrees_east", lon[lon_i], longname="longitude")
  latdim <- ncdf4::ncdim_def("lat", "degrees_north", lat[lat_i], longname="latitude")
  timedim <- ncdf4::ncdim_def("month", "", as.integer(1:12), longname="months")

  fillvalue <- 1e30

  varnames <- sapply(1:nc$nvars, function(i, x){x$var[[i]]$name}, nc)

  vars <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 7:length(varnames)){
    utils::setTxtProgressBar(pb, i)
    varname <- nc$var[[i]]$name
    varunits <- nc$var[[i]]$units
    vardim <- list(latdim, londim, timedim)
    varlongname <- nc$var[[i]]$longname
    varprec <- nc$var[[i]]$prec
    varshuffle <- nc$var[[i]]$suffle
    if(is.null(varshuffle)){ varshuffle <- FALSE }
    varcompression <- nc$var[[i]]$compression
    varchunksizes <- c(length(lat_i), length(lon_i), length(mon_i))
    vars[[i-6]] <- ncdf4::ncvar_def(varname, varunits, vardim, fillvalue, varlongname, varprec, varshuffle, varcompression, varchunksizes)
  }
  close(pb)
  
  .check_and_create_directory(out.path)

  nc.target <- .open_or_create_ncfile(outpathfile, vars, suppress_dimvals = supress_dimvals)

  pb <- utils::txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 7:length(varnames)){
    utils::setTxtProgressBar(pb, i)
    varvals <- ncdf4::ncvar_get(nc, varnames[[i]])
    varvals <- varvals[lon_i, lat_i, mon_i]
    varvals <- aperm(varvals, c(2, 1, 3))
    ncdf4::ncvar_put(nc.target, varnames[[i]], varvals)
  }
  close(pb)
  
  ncdf4::nc_close(nc.target)
  ncdf4::nc_close(nc)
}

#' Title
#'
#' @param nc.source
#' @param nc.baseline
#' @param baseline
#' @param out.path
#'
#' @return
#' @export
#'
#' @examples
calculate_anomalies <- function(nc.source, nc.baseline, baseline, out.path, overwrite = FALSE, suppress_dimvals = TRUE){
  # nc.source <- "precipitation-22000BP-15000BP.nc"
  # nc.baseline <- "precipitation-5000BP-1989AD.nc"
  # baseline <- "1951AD-1989AD/1989AD"
  # out.path <- "anomalies"

  in.filename <- basename(nc.source)
  
  outpathfile <- paste(out.path, in.filename, sep="/")

  .check_existing_file(outpathfile, overwrite)
  
  # Extract values of baseline conditions
  nc.bl <- ncdf4::nc_open(nc.baseline, suppress_dimvals = supress_dimvals)
  var.bl <- ncdf4::ncvar_get(nc.bl, varid=baseline)

  nc.src <- ncdf4::nc_open(nc.source, suppress_dimvals = supress_dimvals)
  varnames <- names(nc.src$var)
  vars <- list()
  
  pb <- utils::txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 1:length(varnames)){
    utils::setTxtProgressBar(pb, i)
    varname <- nc.src$var[[i]]$name
    varunits <- nc.src$var[[i]]$units
    vardim <- nc.src$var[[i]]$dim
    varfillvalue <- nc.src$var[[i]]$missval
    varlongname <- nc.src$var[[i]]$longname
    varprec <- nc.src$var[[i]]$prec
    varshuffle <- nc.src$var[[i]]$suffle
    if(is.null(varshuffle)){ varshuffle <- FALSE }
    varcompression <- nc.src$var[[i]]$compression
    varchunksizes <- nc.src$var[[i]]$chunksizes
    vars[[i]] <- ncdf4::ncvar_def(varname, varunits, vardim, varfillvalue, varlongname, varprec, varshuffle, varcompression, varchunksizes)
  }
  close(pb)

  .check_and_create_directory(out.path)
  
  nc.trg <- .open_or_create_ncfile(outpathfile, vars, suppress_dimvals = supress_dimvals)
 
  # Extract values of all time periods in the file and compute anomalies
  pb <- utils::txtProgressBar(min = 0, max = length(nc.src$var), initial = 1, style=3)
  for(i in 1:length(nc.src$var)){
    utils::setTxtProgressBar(pb, i)
    var <- ncdf4::ncvar_get(nc.src, varid=nc.src$var[[i]]$name)
    var.an <- var - var.bl
    var.an <- round(var.an, 3)
    ncdf4::ncvar_put(nc.trg, varid=nc.src$var[[i]]$name, var.an)
  }
  close(pb)

  ncdf4::nc_close(nc.bl)
  ncdf4::nc_close(nc.src)
  ncdf4::nc_close(nc.trg)
}



#' Title
#'
#' @param nc.source
#' @param out.path
#' @param downscale.factor
#' @param res.src
#'
#' @return
#' @export
#'
#' @examples
interpolate_paleoview <- function(nc.source, out.path, res.src = 2.5, downscale.factor = 1/5, flush.seq = seq(1000, 10000, by=1000), overwrite = FALSE, suppress_dimvals = TRUE){
  # nc.source <- "../Data/PaleoView/NetCDF/croped/anomalies/precipitation-22000BP-15000BP.nc"
  # out.path <- "../Data/PaleoView/NetCDF/croped/anomalies/interpolated"
  # res.src <- 2.5
  # downscale.factor <- 1/60
  # flush.seq <- seq(1000, 10000, by=1000)
  # overwrite <- FALSE

  in.filename <- basename(nc.source)
  
  outpathfile <- paste(out.path, in.filename, sep="/")

  .check_existing_file(outpathfile, overwrite)
  
  nc.src <- ncdf4::nc_open(nc.source, suppress_dimvals = supress_dimvals)
  varnames <- names(nc.src$var)
  vars <- list()

  res.trg <- res.src * downscale.factor
  half.res.src <- res.src / 2
  half.res.trg <- res.trg / 2

  lon <- ncdf4::ncvar_get(nc.src, "lon")
  lon.bb <- c(min(lon) - half.res.src, max(lon) + half.res.src)
  lon.trg <- c(lon.bb[1] + half.res.trg, lon.bb[2] - half.res.trg)
  lon.trg <- seq(lon.trg[1], lon.trg[2], by=res.trg)

  lat <- ncdf4::ncvar_get(nc.src,"lat")
  lat.bb <- c(min(lat) - half.res.src, max(lat) + half.res.src)
  lat.trg <- c(lat.bb[1] + half.res.trg, lat.bb[2] - half.res.trg)
  lat.trg <- seq(lat.trg[1], lat.trg[2], by=res.trg)

  mon <- ncdf4::ncvar_get(nc.src, "month")

  londim <- ncdf4::ncdim_def("lon", "degrees_east", lon.trg, longname="longitude")
  latdim <- ncdf4::ncdim_def("lat", "degrees_north", lat.trg, longname="latitude")
  timedim <- ncdf4::ncdim_def("month", "", as.integer(1:12), longname="months")

  cat("Creating variables for the new data file:", "\n")
  pb <- utils::txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 1:length(varnames)){
    utils::setTxtProgressBar(pb, i)
    varname <- nc.src$var[[i]]$name
    varunits <- nc.src$var[[i]]$units
    vardim <- list(latdim, londim, timedim)
    varfillvalue <- nc.src$var[[i]]$missval
    varlongname <- nc.src$var[[i]]$longname
    varprec <- nc.src$var[[i]]$prec
    varshuffle <- nc.src$var[[i]]$suffle
    if(is.null(varshuffle)){ varshuffle <- FALSE }
    varcompression <- nc.src$var[[i]]$compression
    varchunksizes <- c(length(lat.trg), length(lon.trg), length(mon))
    vars[[i]] <- ncdf4::ncvar_def(varname, varunits, vardim, varfillvalue, varlongname, varprec, varshuffle, varcompression, varchunksizes)
  }
  close(pb)

  .check_and_create_directory(out.path)
  
  nc.trg <- .open_or_create_ncfile(outpathfile, vars, suppress_dimvals = supress_dimvals)

  raster.trg <- raster::brick(nrows=length(lat.trg), ncols=length(lon.trg), xmn=lon.bb[1], xmx=lon.bb[2], ymn=lat.bb[1], ymx=lat.bb[2], nl=length(mon))

  cat("Interpolating variables to the new resolution:", "\n")
  pb <- utils::txtProgressBar(min = 0, max = length(nc.src$var), initial = 1, style=3)
  for(i in 1:length(nc.src$var)){
    utils::setTxtProgressBar(pb, i)
    var <- ncdf4::ncvar_get(nc.src, varid=nc.src$var[[i]]$name)
    var <- raster::brick(var)
    var <- raster::setExtent(var, c(lon.bb[1], lon.bb[2], lat.bb[1], lat.bb[2]))
    var.trg <- raster::resample(var, raster.trg)
    var.trg <- as.vector(raster::t(var.trg))
    var.trg <- round(var.trg, 3)
    ncdf4::ncvar_put(nc.trg, varid=nc.trg$var[[i]]$name, var.trg)
    if(i %in% flush.seq){
      nc.trg <- .flush_ncfile(nc.trg, suppress_dimvals = supress_dimvals)
      # nc_sync(nc.trg) # Not use this. It doesn't free memory.
      # ncdf4::nc_close(nc.trg)
      # nc.trg <- ncdf4::nc_open(outpathfile, write=T, suppress_dimvals = supress_dimvals)
    }
  }
  close(pb)

  ncdf4::nc_close(nc.src)
  ncdf4::nc_close(nc.trg)
}

#' Title
#'
#' @param file
#' @param var
#' @param ext
#'
#' @return
#' @export
#'
#' @examples
rasterize_paleoview <- function(file, var = names(nc$var)[[1]], ext, suppress_dimvals = TRUE){
  cat(file, "\n")
  nc <- ncdf4::nc_open(file, suppress_dimvals = supress_dimvals)
  if(is.numeric(var) == T){
    var <- names(nc$var)[[var]]
  }
  if(!is.character(var)){
    stop("'var' has to be numeric or character.")
  }
  cat(var, "\n")
  data <- ncdf4::ncvar_get(nc, var)
  v <- raster::brick(data, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4])
  ncdf4::nc_close(nc)
  return(v)
}
