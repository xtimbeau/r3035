#' Raster vers data.table
#'
#' Transforme un raster::ratser en data.table
#' en ajoutant un idINS. Fonction inverse de dt2r
#'
#' @param raster le raster \code{RatserLayer}, \code{RatserBrick},  \code{RatserStack}
#' @param resolution résolution en mètre pour la transformation (200m par défaut)
#' @param fun une fonction pour l'agrégation (mean par défaut)
#' @import data.table
#' @return un data.table avec un idINS
#' @export
#'
r2dt <- function(raster, resolution=NULL, fun=mean, shortidINS = FALSE)
{
  base_res <- max(raster::res(raster))
  vars <- names(raster)
  dt <- raster::as.data.frame(raster, xy=TRUE, centroids=TRUE)
  data.table::setDT(dt)
  dt <- na.omit(melt(dt, measure.vars=vars), "value")
  dt <- data.table::dcast(dt, x+y~variable, value.var="value")
  dt[, idINS := idINS3035(x, y, resolution=base_res)]
  id <- stringr::str_c("idINS", base_res)
  data.table::setnames(dt, "idINS",id)
  navars <- setdiff(vars, names(dt))
  rvars <- setdiff(vars, navars)
  if(!is.null(resolution))
  {
    if(shortidINS)
      ids <- sidINS3035(x,y,resolution) else
        ids <- idINS3035(x,y,resolution)

    id <- stringr::str_c("idINS",resolution)
    dt[, (stringr::str_c("idINS",resolution)):=ids]
    dt <- dt[, lapply(.SD, function(x) fun(x, na.rm=TRUE)), by=c(id), .SDcols=rvars]
  }
  if (length(navars)>0)
    dt[, (navars):=rep(list(rep(NA, nrow(dt))), length(navars))]
  data.table::setkeyv(dt, cols=id)
  dt
}

# vérifie qu'il y a un idINS et détermine la résolution, interne
getresINS <- function(dt, idINS="idINS") {
  purrr::map(
    purrr::keep(names(dt), ~stringr::str_detect(.x,idINS)),
    ~{
      r<-stringr::str_extract(dt[[.x]], "(?<=r)[0-9]+")  |>
        as.numeric()
      ur <- unique(r)
      if (length(ur)==0)
        list(idINS=.x, res=NA_integer_)
      else
        list(idINS=.x, res=ur)
    })
}

# vérifie qu'il y a un idINS et détermine la résolution, interne
getINSres <- function(dt, resolution, idINS="idINS") {
  rr <- getresINS(dt, idINS)
  ncol <- names(dt)
  if(length(rr)==0)
    return(FALSE)
  if (resolution==Inf)
    resolution <- min(unlist(purrr::transpose(rr)$res))
  isresin <- purrr::map_lgl(rr, ~.x[["res"]]==resolution)
  if(any(isresin))
    return(list(idINS = rr[which(isresin)][[1]]$idINS, res = resolution))
  else
    return(FALSE)
}

#' Crée un raster à partir d'un data.table avec idINS.
#'
#' @param dt data.table avec idINS.
#' @param resolution résolution du raster.
#' @param idINS nom de la variable idINS, par défaut "idINS".
#'
#' @import data.table
#'
#' @export
dt2r <- function (dt, resolution = NULL, idINS = "idINS", shortidINS = FALSE)
{
  dt <- setDT(dt)
  if(!shortidiNS) rr <- r3035:::getresINS(dt, idINS)
  if(shortidiNS) {
    rr <- NULL
    resolution <- 200}
  ncol <- names(dt)
  if (length(rr) == 0)
    idINSin <- FALSE
  else {
    if (!is.null(resolution))
      isresin <- purrr::map_lgl(rr, ~.x[["res"]] == resolution)
    else isresin <- c(TRUE, rep(FALSE, length(rr) - 1))
    if (any(isresin)) {
      res <- rr[which(isresin)][[1]]$res
      idINS <- rr[which(isresin)][[1]]$idINS
      idINSin <- TRUE
    }
    else idINSin <- FALSE
  }
  if (!idINSin) {
    stopifnot(!is.null(res))
    stopifnot("x" %in% ncol & "y" %in% ncol)
    if(shortidINS)
      dt[, `:=`(idINS, sidINS3035(x, y, resolution = resolution))]
    else
      dt[, `:=`(idINS, idINS3035(x, y, resolution = resolution))]
    idINS <- "idINS"
    res <- resolution
  }
  if(shortidINS)
    xy <- sidINS2point(dt[[idINS]], resolution = res)
  if(!shortidINS)
    xy <- idINS2point(dt[[idINS]], resolution = res)

  dt[, `:=`(x = xy[, 1], y = xy[, 2])]
  rref <- raster_ref(dt, resolution = res, crs = 3035)
  cells <- raster::cellFromXY(rref, xy)
  layers <- purrr::keep(ncol, ~(!.x %in% c("x", "y", idINS)))
  brickette <- raster::brick(purrr::map(layers, ~{
    r <- raster::raster(rref)
    r[cells] <- dt[[.x]]
    r
  }))
  names(brickette) <- layers
  raster::crs(brickette) <- 3035
  brickette
}

# fabrique un raster de référence, vide
#
#' Title
#'
#' @param data un sf ou n'importe quoi qui répond à st_bbox
#' @param resolution la résolution
#' @param crs le systeme de coordonnées
#'
#' @return un raster
raster_ref <- function(data, resolution=200, crs=3035)
{
  alignres <- max(resolution, 200)
  if(checkmate::testMultiClass(data,c("sf", "sfc")))
  {
    b <- sf::st_bbox(data)
    crss <- sf::st_crs(data)$proj4string
  }
  else
  {
    stopifnot("x"%in%names(data)&"y"%in%names(data))
    b <- list(xmin=min(data$x, na.rm=TRUE),
              xmax=max(data$x, na.rm=TRUE),
              ymin=min(data$y, na.rm=TRUE),
              ymax=max(data$y, na.rm=TRUE))
    crss <- sp::CRS(st_crs(crs)$proj4string)
  }
  ext <- raster::extent(floor(b$xmin / alignres )*alignres,
                        ceiling(b$xmax/alignres)*alignres,
                        floor(b$ymin/alignres)*alignres,
                        ceiling(b$ymax/alignres)*alignres)
  raster::raster(ext, crs=crss,  resolution=resolution)
}


#' Transform a data.table with idINS to a star object.
#' Assumes crs 3035
#'
#' @param data a data.table with an idINS column
#' @param idINS name of the idINS column. Default to idINS.
#'
#' @import data.table
#'
#' @export
dt2stars <- function(data, idINS = "idINS")
{
  if(!is.data.table(data)) setDT(data)
  xy <- idINS2coord(data[[idINS]]) |> as.data.table()

  data <- cbind(data[, .SD, .SDcol = !idINS], xy)

  #sf::st_as_sf(data, coords = c("x", "y"), crs = st_crs(3035))
  data <- stars::st_as_stars(data, dims = c("x", "y"))
  st_crs(data) <- sf::st_crs(3035)

  data
}
