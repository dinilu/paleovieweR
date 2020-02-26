#' Title
#'
#' @param nc.source
#' @param baseline
#' @param out.path
#'
#' @return
#' @export
#'
#' @examples
# combine_anomalies_and_baseline <- function(nc.source, baseline, out.path, init.i=NULL, end.i=NULL){
combine_anomalies_and_baseline <- function(nc.source, baseline, out.path, flush.seq = seq(1000, 10000, by=1000), overwrite = FALSE, suppress_dimvals = TRUE){
  # nc.source <- "precipitation-22000BP-15000BP.nc"
  # baseline <- prec
  # out.path <- "/home/dinilu/Med-Refugia/Data/PaleoView_Downscaled"

  nc.src <- ncdf4::nc_open(nc.source, suppress_dimvals = suppress_dimvals)
  varnames <- names(nc.src$var)

  pathfile <- paste(out.path, nc.source, sep="/")
  if(file.exists(pathfile)){
    if(overwrite == FALSE){
      stop("File already exists and overwrite option set to FALSE (default), if you want to overwrite turn the option to TRUE.")
    }
  }

  # if(is.null(init.i)){ init.i <- 1}
  # if(is.null(end.i)){ end.i <- length(nc.src$var)}

  lon <- ncdf4::ncvar_get(nc.src, "lon")
  half.res <- ((lon[[2]] - lon[[1]])/2)
  lon.bb <- c(min(lon) - half.res, max(lon) + half.res)
  lat <- ncdf4::ncvar_get(nc.src, "lat")
  half.res <- ((lat[[2]] - lat[[1]])/2)
  lat.bb <- c(min(lat) - half.res, max(lat) + half.res)

  cat("Creating variables for the new data file:", "\n")
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

  if(!dir.exists(out.path)){
    cat("Creating new output directory.", "\n")
    dir.create(out.path)
  }

  wd.init <- getwd()
  setwd(out.path)

  # if(init.i != 1){
  #   cat("Opening NetCDF file.", "\n")
  #   nc.trg <- nc_open(nc.source, write=TRUE, suppress_dimvals = suppress_dimvals)
  # }else{
  #   cat("Creating new file.", "\n")
  #   nc.trg <- nc_create(nc.source, vars)
  # }

  if(file.exists(nc.source)){
    cat("Opening existing NetCDF file.", "\n")
    nc.trg <- ncdf4::nc_open(nc.source, write = TRUE, suppress_dimvals = suppress_dimvals)
  }else{
    cat("Creating new NetCDF file.", "\n")
    nc.trg <- ncdf4::nc_create(nc.source, vars)
  }

  cat("Adding anomalies to baseline:", "\n")
  pb <- utils::txtProgressBar(min = 0, max = length(nc.src$var), initial = 1, style=3)
  for(i in 1:length(nc.src$var)){
    # pb <- txtProgressBar(min = init.i, max = end.i, initial = 1, style=3)
    # for(i in init.i:end.i){
    utils::setTxtProgressBar(pb, i)
    var <- ncdf4::ncvar_get(nc.src, varid=nc.src$var[[i]]$name)
    var <- raster::brick(var)
    var <- raster::setExtent(var, c(lon.bb[[1]], lon.bb[[2]], lat.bb[[1]], lat.bb[[2]]))

    if(raster::extent(var) != raster::extent(baseline) || raster::res(var)[1] != raster::res(baseline)[1] || raster::res(var)[2] != raster::res(baseline)[2]){
      stop("NetCDF files and baseline do not match spatially (different extent or resolution).")
    }

    if(nc.src$var[[i]]$units == "degrees C"){
      var <- baseline + var
      var <- as.vector(t(var))
    }
    if(nc.src$var[[i]]$units == "mm/day"){
      var <- var * (baseline + 1)
      var <- as.vector(t(var))
    }
    var <- round(var, 3)
    ncdf4::ncvar_put(nc.trg, varid=nc.trg$var[[i]]$name, var)
    if(i %in% flush.seq){
      nc.trg <- .flush_ncfile(nc.trg, suppress_dimvals = suppress_dimvals)
      # ncdf4::nc_close(nc.trg)
      # nc.trg <- ncdf4::nc_open(nc.source, write=T, suppress_dimvals = suppress_dimvals)
    }
  }
  close(pb)

  cat("Closing files.", "\n")
  ncdf4::nc_close(nc.src)
  ncdf4::nc_close(nc.trg)
  setwd(wd.init)
}

