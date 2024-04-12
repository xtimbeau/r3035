#' convertit les idINS longs en courts
#'
#' @param idINS les id longs
#'
#' @return des id courts (100000*x/r + y/r)
#' @export
#'
contract_idINS <- function(idINS) {
  cr_pos <- stringr::str_locate(idINS[[1]], "r(?=[0-9])")[,"start"]+1
  cy_pos <- stringr::str_locate(idINS[[1]], "N(?=[0-9])")[,"start"]+1
  cx_pos <- stringr::str_locate(idINS[[1]], "E(?=[0-9])")[,"start"]+1
  lcoord <- cx_pos-cy_pos-1
  r <- as.integer(stringr::str_sub(idINS[[1]],cr_pos,cy_pos-cr_pos))
  y <- as.integer(stringr::str_sub(idINS,cy_pos,cy_pos+lcoord)) %/% r
  x <- as.integer(stringr::str_sub(idINS,cx_pos,cx_pos+lcoord)) %/% r
  as.integer(x * 100000 + y)
}

#' convertit les idINS courts en longs
#'
#' @param ids les id courts
#' @param resolution la résolution (defaut 200m)
#'
#' @return un vecteur d'idINS longs (ex. r200N2051600E4258800)
#' @export
#'
expand_idINS <- function(ids, resolution=200) {
  x <- ids%/%100000L
  y <- as.integer(resolution)*(ids-x*100000L)
  x <- as.integer(resolution)*x
  stringr::str_c("r", as.integer(resolution), "N", as.integer(y), "E", as.integer(x))
}

#' Convertit l'idINS en polygone
#'
#' @param ids vecteur d'idINS.
#' @param resolution resolution, par défaut celle attachée à idINS.
#'
#' @export
#'
sidINS2square <- function(ids, resolution=200)
{
  ids <- as.integer(ids)
  resolution <- as.integer(resolution)
  x <- ids%/%100000L
  y <- resolution*(ids-x*100000L)
  x <- resolution*x
  purrr::pmap(list(x,y), ~sf::st_polygon(
    list(matrix(
      c(..1, ..1+resolution,
        ..1+resolution,..1,
        ..1, ..2,
        ..2, ..2+resolution,
        ..2+resolution,..2),
      nrow=5, ncol=2))))  |>
    sf::st_sfc(crs=3035)
}

#' Calculate euclidean distance between to short idINS, in meter
#'
#' @param fromidINS character vector of starting idINS
#' @param toidINS character vector of ending idINS
#' @param resolution default to NULL. Set if no resolution is provided in idINS
#'
#' @export
sidINS2dist <- function(fromidINS, toidINS, resolution=200) {
  stopifnot(length(fromidINS)==length(toidINS))
  if(length(fromidINS)==0)
    return(numeric())

  fromx <- as.integer(fromidINS)%/%100000
  fromy <- as.integer(fromidINS)-fromx*100000

  tox <- as.integer(toidINS)%/%100000
  toy <- as.integer(toidINS)-tox*100000

  return(resolution*sqrt((tox-fromx)^2 + (toy-fromy)^2))
}

#' retourne un sf à partir d'un data.frame avec un champ idINS
#'
#' @param data le data frame
#' @param resolution 200m par défaut
#' @param idINS le champ qui contient l'idINS court (défaut "idINS")
#'
#' @return un sf
#' @export
#'
sidINS2sf <- function(data, resolution = 200, idINS = "idINS") {
  data$geometry <- data[[idINS]] |>
    sidINS2square()
  sf::st_as_sf(data)
}

#' Crée idINS à partir de coordonnées.
#'
#' @param x Ou un vecteur de la coordonnée x, ou un df avec les colonnes x et y.
#' @param y vecteur de la coordonnée y, si nécessaire. NULL par défaut.
#' @param resolution resolution, par défaut 200.
#' @param resinstr faut-il inscrire la résolution dans idINS ? Par défaut TRUE.
#'
#' @export
sidINS3035 <- function(x, y=NULL, resolution=200)
{
  if(is.null(y))
  {
    y <- x[,2]
    x <- x[,1]
  }

  x <- as.integer(floor(x / resolution ))
  y <- as.integer(floor(y / resolution ))
  resultat <- x*100000L + y

  return(resultat)
}

#' Récupère les coordonnées X et Y de idINS court.
#'
#' @param ids vecteur d'idINS.
#' @param resolution resolution, par défaut celle attachée à idINS.
#'
#' @export
sidINS2point <- function(ids, resolution=200)
{
  x <- round(ids/100000)
  y <- round(resolution*(ids-x*100000))
  x <- round(resolution*x)
  r <- resolution
  m <- matrix(c(as.integer(x+r/2),as.integer(y+r/2)), ncol=2)
  colnames(m) <- c("X", "Y")
  m
}

#' Raster vers data.table avec idINS short
#'
#' Transforme un raster::ratser en data.table
#' en ajoutant un idINS. Fonction inverse de sdt2r
#'
#' @param raster le raster \code{RatserLayer}, \code{RatserBrick},  \code{RatserStack}
#' @param resolution résolution en mètre pour la transformation (200m par défaut)
#' @param fun une fonction pour l'agrégation (mean par défaut)
#' @import data.table
#' @return un data.table avec un idINS
#' @export
#'
r2sdt <- function(raster, resolution=200, fun=mean)
{
  base_res <- resolution
  vars <- names(raster)
  dt <- raster::as.data.frame(raster, xy=TRUE, centroids=TRUE)
  data.table::setDT(dt)
  dt <- na.omit(melt(dt, measure.vars=vars), "value")
  dt <- data.table::dcast(dt, x+y~variable, value.var="value")
  dt[, idINS := sidINS3035(x, y, resolution=base_res)]
  navars <- setdiff(vars, names(dt))
  rvars <- setdiff(vars, navars)
  if(!is.null(resolution))
  {
    id <- "idINS"
    dt[, idINS:=sidINS3035(x,y,resolution)]
    dt <- dt[, lapply(.SD, function(x) fun(x, na.rm=TRUE)), by=c(id), .SDcols=rvars]
  }
  if (length(navars)>0)
    dt[, (navars):=rep(list(rep(NA, nrow(dt))), length(navars))]
  data.table::setkeyv(dt, cols=id)
  dt
}

#' Crée un raster à partir d'un data.table avec idINS short
#'
#' @param dt data.table avec idINS.
#' @param resolution résolution du raster.
#' @param idINS nom de la variable idINS, par défaut "idINS".
#'
#' @import data.table
#'
#' @export
sdt2r <- function (dt, resolution = 200, idINS = "idINS")
{
  dt <- setDT(dt)
  ncol <- names(dt)
  stopifnot(!is.null(res))
  res <- resolution
  xy <- sidINS2point(dt[[idINS]], resolution = res)
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

#' Retrieves (lon, lat) coordinates from long idINS.
#'
#' @param idINS character vector of idINS
#' @param resolution default to NULL. Set to convert id with different resolutions
#'
#' @export
sidINS2lonlat <- function(idINS, resolution=200) {
  if(length(idINS)==0)
    return(tibble::tibble(lon=numeric(), lat=numeric()))
  r <- as.integer(resolution)
  x <- idINS%/%100000L
  y <- r*(idINS-x*100000L)
  x <- r*x
  x <- x+r/2
  y <- y+r/2
  lonlat <- sf_project(
    from=st_crs(3035),
    to=st_crs(4326),
    pts = matrix(c(x,y), nrow=length(x), ncol=2))
  return(
    tibble(lon = lonlat[, 1], lat = lonlat[, 2]))
}

#' Retrieves (lon, lat) coordinates from long or short idINS.
#'
#' @param idINS character or integer vector of idINS
#' @param resolution default to NULL. Set to convert id with different resolutions
#'
#' @export
idINS2lonlat <- function(idINS, resolution=NULL) {
  if(is.character(idINS))
    return(lidINS2lonlat(idINS, resolution))
  if(is.null(resolution))
    resolution <- 200
  if(is.integer(idINS))
    return(sidINS2lonlat(idINS, resolution))
  stop("idINS must be character or integer")
}

#' Récupère les coordonnées X et Y de idINS court ou long.
#'
#' @param ids vecteur d'idINS.
#' @param idINS le champ qui contient l'idINS court (défaut "idINS")
#' @param resolution resolution, par défaut celle attachée à idINS.
#'
#' @export
idINS2point <- function(ids, resolution=NULL)
{
  if(is.character(ids))
    return(lidINS2point(ids, resolution=resolution))
  if(is.null(resolution))
    resolution <- 200
  if(is.integer(ids))
    return(sidINS2point(ids, resolution=resolution))
  stop("idINS must be character or integer")
}

#' Récupère les coordonnées X et Y de idINS court ou long.
#'
#' @param ids vecteur d'idINS.
#' @param resolution resolution, par défaut celle attachée à idINS.
#'
#' @export
idINS2square <- function(ids, resolution=NULL)
{
  if(is.character(idINS))
    return(lidINS2square(idINS, resolution))
  if(is.null(resolution))
    resolution <- 200
  if(is.integer(idINS))
    return(sidINS2square(idINS, resolution))
  stop("idINS must be character or integer")
}

#' Calculate euclidean distance between to long or short idINS, in meter
#'
#' @param fromidINS character vector of starting idINS
#' @param toidINS character vector of ending idINS
#' @param resolution default to NULL. Set if no resolution is provided in idINS
#'
#' @export
#'

idINS2dist <- function(fromidINS, toidINS, resolution=NULL) {
  if(is.character(fromidINS)&is.character(toidINS))
    return(lidINS2dist(fromidINS, toidINS, resolution))
  if(is.null(resolution))
    resolution <- 200
  if(is.integer(fromidINS)&is.integer(toidINS))
    return(sidINS2dist(fromidINS, toidINS, resolution))
  stop("idINS must be character or integer")
}

#' retourne un sf à partir d'un data.frame avec un champ idINS
#'
#' @param data le data frame
#' @param resolution 200m par défaut
#' @param idINS le champ qui contient l'idINS court (défaut "idINS")
#'
#' @return un sf
#' @export
#'
idINS2sf <- function(data, resolution = 200, idINS = "idINS") {
  if(is.character(idINS))
    return(lidINS2sf(data, resolution, idINS))
  if(is.null(resolution))
    resolution <- 200
  if(is.integer(idINS))
    return(sidINS2sf(data, resolution, idINS))
  stop("idINS must be character or integer")
  }
