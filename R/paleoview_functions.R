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

.open_or_create_ncfile <- function(nc.source, vars) {
  if(file.exists(nc.source)){
    cat("Opening existing NetCDF file.", "\n")
    nc.trg <- nc_open(nc.source, write = TRUE)
  }else{
    cat("Creating new NetCDF file.", "\n")
    nc.trg <- nc_create(nc.source, vars)
  }
  return(nc.trg)
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
crop_paleoview <- function(nc.source, ext, out.path, overwrite = FALSE){
  # nc.source <- "precipitation-5000BP-1989AD.nc"
  # ext <- study_area
  # out.path <- "croped"
  # overwrite <- FALSE
 
  in.path <- dirname(nc.source)
  nc.source <- basename(nc.source)
  
  outpathfile <- paste(out.path, nc.source, sep="/")
  inpathfile <- paste(in.path, nc.source, sep="/")
  
  .check_existing_file(outpathfile, overwrite)
  
  nc <- nc_open(inpathfile)
  lon <- ncvar_get(nc, "longitudes")
  lon_i <- which(lon >= ext[1] & lon <= ext[2])
  lat <- ncvar_get(nc,"latitudes")
  lat_i <- which(lat >= ext[3] & lat < ext[4])
  mon <- ncvar_get(nc, "months")
  mon_i <- which (mon >= 1 & mon <= 12)

  londim <- ncdim_def("lon", "degrees_east", lon[lon_i], longname="longitude")
  latdim <- ncdim_def("lat", "degrees_north", lat[lat_i], longname="latitude")
  timedim <- ncdim_def("month", "", as.integer(1:12), longname="months")

  fillvalue <- 1e30

  varnames <- sapply(1:nc$nvars, function(i, x){x$var[[i]]$name}, nc)

  vars <- list()
  pb <- txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 7:length(varnames)){
    setTxtProgressBar(pb, i)
    varname <- nc$var[[i]]$name
    varunits <- nc$var[[i]]$units
    vardim <- list(latdim, londim, timedim)
    varlongname <- nc$var[[i]]$longname
    varprec <- nc$var[[i]]$prec
    varshuffle <- nc$var[[i]]$suffle
    if(is.null(varshuffle)){ varshuffle <- FALSE }
    varcompression <- nc$var[[i]]$compression
    varchunksizes <- c(length(lat_i), length(lon_i), length(mon_i))
    vars[[i-6]] <- ncvar_def(varname, varunits, vardim, fillvalue, varlongname, varprec, varshuffle, varcompression, varchunksizes)
  }
  close(pb)
  
  .check_and_create_directory(out.path)

  nc.target <- .open_or_create_ncfile(outpathfile, vars)

  pb <- txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 7:length(varnames)){
    setTxtProgressBar(pb, i)
    varvals <- ncvar_get(nc, varnames[[i]])
    varvals <- varvals[lon_i, lat_i, mon_i]
    varvals <- aperm(varvals, c(2, 1, 3))
    ncvar_put(nc.target, varnames[[i]], varvals)
  }
  close(pb)
  
  nc_close(nc.target)
  nc_close(nc)
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
calculate_anomalies <- function(nc.source, nc.baseline, baseline, out.path, overwrite = FALSE){
  # nc.source <- "precipitation-22000BP-15000BP.nc"
  # nc.baseline <- "precipitation-5000BP-1989AD.nc"
  # baseline <- "1951AD-1989AD/1989AD"
  # out.path <- "anomalies"

  in.path <- dirname(nc.source)
  nc.source <- basename(nc.source)
  
  outpathfile <- paste(out.path, nc.source, sep="/")
  inpathfile <- paste(in.path, nc.source, sep="/")
  basepathfile <- paste(in.path, nc.baseline, sep="/")
  
  .check_existing_file(outpathfile, overwrite)
  
  # Extract values of baseline conditions
  nc.bl <- nc_open(basepathfile)
  var.bl <- ncvar_get(nc.bl, varid=baseline)

  nc.src <- nc_open(inpathfile)
  varnames <- names(nc.src$var)
  vars <- list()
  
  pb <- txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 1:length(varnames)){
    setTxtProgressBar(pb, i)
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
    vars[[i]] <- ncvar_def(varname, varunits, vardim, varfillvalue, varlongname, varprec, varshuffle, varcompression, varchunksizes)
  }
  close(pb)

  .check_and_create_directory(out.path)
  
  nc.trg <- .open_or_create_ncfile(outpathfile, vars)
 
  # Extract values of all time periods in the file and compute anomalies
  pb <- txtProgressBar(min = 0, max = length(nc.src$var), initial = 1, style=3)
  for(i in 1:length(nc.src$var)){
    setTxtProgressBar(pb, i)
    var <- ncvar_get(nc.src, varid=nc.src$var[[i]]$name)
    if(nc.src$var[[i]]$units == "degrees C"){
      var.an <- var - var.bl
    }
    if(nc.src$var[[i]]$units == "mm/day"){
      var.an <- var / (var.bl + 1)
    }
    var.an <- round(var.an, 3)
    ncvar_put(nc.trg, varid=nc.src$var[[i]]$name, var.an)
  }
  close(pb)

  nc_close(nc.bl)
  nc_close(nc.src)
  nc_close(nc.trg)
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
interpolate_paleoview <- function(nc.source, out.path, res.src = 2.5, downscale.factor = 1/5, flush.seq = seq(1000, 10000, by=1000), overwrite = FALSE){
  # nc.source <- "precipitation-22000BP-15000BP.nc"
  # out.path <- "interpolated"
  # downscale.factor <- 0.2
  # res.src <- 2.5

  in.path <- dirname(nc.source)
  nc.source <- basename(nc.source)
  
  outpathfile <- paste(out.path, nc.source, sep="/")
  inpathfile <- paste(in.path, nc.source, sep="/")

  .check_existing_file(outpathfile, overwrite)
  
  nc.src <- nc_open(inpathfile)
  varnames <- names(nc.src$var)
  vars <- list()

  res.trg <- res.src * downscale.factor
  half.res.src <- res.src / 2
  half.res.trg <- res.trg / 2

  lon <- ncvar_get(nc.src, "lon")
  lon.bb <- c(min(lon) - half.res.src, max(lon) + half.res.src)
  lon.trg <- c(lon.bb[1] + half.res.trg, lon.bb[2] - half.res.trg)
  lon.trg <- seq(lon.trg[1], lon.trg[2], by=res.trg)

  lat <- ncvar_get(nc.src,"lat")
  lat.bb <- c(min(lat) - half.res.src, max(lat) + half.res.src)
  lat.trg <- c(lat.bb[1] + half.res.trg, lat.bb[2] - half.res.trg)
  lat.trg <- seq(lat.trg[1], lat.trg[2], by=res.trg)

  mon <- ncvar_get(nc.src, "month")

  londim <- ncdim_def("lon", "degrees_east", lon.trg, longname="longitude")
  latdim <- ncdim_def("lat", "degrees_north", lat.trg, longname="latitude")
  timedim <- ncdim_def("month", "", as.integer(1:12), longname="months")

  cat("Creating variables for the new data file:", "\n")
  pb <- txtProgressBar(min = 0, max = length(varnames), initial = 1, style=3)
  for(i in 1:length(varnames)){
    setTxtProgressBar(pb, i)
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
    vars[[i]] <- ncvar_def(varname, varunits, vardim, varfillvalue, varlongname, varprec, varshuffle, varcompression, varchunksizes)
  }
  close(pb)

  .check_and_create_directory(out.path)
  
  nc.trg <- .open_or_create_ncfile(outpathfile, vars)

  raster.trg <- brick(nrows=length(lat.trg), ncols=length(lon.trg), xmn=lon.bb[1], xmx=lon.bb[2], ymn=lat.bb[1], ymx=lat.bb[2], nl=length(mon))

  cat("Interpolating variables to the new resolution:", "\n")
  pb <- txtProgressBar(min = 0, max = length(nc.src$var), initial = 1, style=3)
  for(i in 1:length(nc.src$var)){
    setTxtProgressBar(pb, i)
    var <- ncvar_get(nc.src, varid=nc.src$var[[i]]$name)
    var <- brick(var)
    extent(var) <- c(lon.bb[1], lon.bb[2], lat.bb[1], lat.bb[2])
    var.trg <- resample(var, raster.trg)
    var.trg <- as.vector(t(var.trg))
    var.trg <- round(var.trg, 3)
    ncvar_put(nc.trg, varid=nc.trg$var[[i]]$name, var.trg)
    if(i %in% flush.seq){
      # nc_sync(nc.trg) # Not use this. It doesn't free memory.
      nc_close(nc.trg)
      nc.trg <- nc_open(nc.source, write=T)
    }
  }
  close(pb)

  nc_close(nc.src)
  nc_close(nc.trg)
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
rasterize_paleoview <- function(file, var = names(nc$var)[[1]], ext){
  cat(file, "\n")
  nc <- nc_open(file)
  if(is.numeric(var) == T){
    var <- names(nc$var)[[var]]
  }
  if(!is.character(var)){
    stop("'var' has to be numeric or character.")
  }
  cat(var, "\n")
  data <- ncvar_get(nc, var)
  v <- brick(data, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4])
  nc_close(nc)
  return(v)
}
