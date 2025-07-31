#' Convertit l'idINS en polygone
#'
#' @param ids vecteur d'idINS.
#' @param resolution resolution, par défaut celle attachée à idINS.
#'
#' @export
#'
lidINS2square <- function(ids, resolution=NULL)
{
  cr_pos <- stringr::str_locate(ids[[1]], "r(?=[0-9])")[,"start"]+1
  cy_pos <- stringr::str_locate(ids[[1]], "N(?=[0-9])")[,"start"]+1
  cx_pos <- stringr::str_locate(ids[[1]], "E(?=[0-9])")[,"start"]+1
  lcoord <- cx_pos-cy_pos-2
  y <- as.numeric(stringr::str_sub(ids,cy_pos,cy_pos+lcoord))
  x <- as.numeric(stringr::str_sub(ids,cx_pos,cx_pos+lcoord))
  r <- if(is.null(resolution))
    as.numeric(stringr::str_sub(ids,cr_pos,cy_pos-cr_pos))
  else
    rep(resolution, length(x))
  purrr::pmap(list(x,y,r), ~sf::st_polygon(
    list(matrix(
      c(..1, ..1+..3,..1+..3,..1, ..1,
        ..2, ..2, ..2+..3,..2+..3,..2),
      nrow=5, ncol=2))))  |>
    sf::st_sfc(crs=3035)
}

#' Récupère les coordonnées X et Y de idINS long
#'
#' @param ids vecteur d'idINS.
#' @param resolution resolution, par défaut celle attachée à idINS.
#'
#' @export
lidINS2point <- function(ids, resolution=NULL)
{
  cr_pos <- stringr::str_locate(ids[[1]], "r(?=[0-9])")[,"start"]+1
  cy_pos <- stringr::str_locate(ids[[1]], "N(?=[0-9])")[,"start"]+1
  cx_pos <- stringr::str_locate(ids[[1]], "E(?=[0-9])")[,"start"]+1
  lcoord <- cx_pos-cy_pos-2
  y <- as.numeric(stringr::str_sub(ids,cy_pos,cy_pos+lcoord))
  x <- as.numeric(stringr::str_sub(ids,cx_pos,cx_pos+lcoord))
  r <- if(is.null(resolution))
    as.numeric(stringr::str_sub(ids,cr_pos,cy_pos-cr_pos))
  else
    rep(resolution, length(x))
  m <- matrix(c(x+r/2,y+r/2), ncol=2)
  colnames(m) <- c("X", "Y")
  m
}

#' Crée idINS à partir de coordonnées.
#'
#' @param x Ou un vecteur de la coordonnée x, ou un df avec les colonnes x et y.
#' @param y vecteur de la coordonnée y, si nécessaire. NULL par défaut.
#' @param resolution résolution, par défaut 200.
#' @param resinstr faut-il inscrire la résolution dans idINS ? Par défaut TRUE.
#'
#' @export
idINS3035 <- function(x, y=NULL, resolution=200, resinstr=TRUE)
{
  if(is.null(y))
  {
    y <- x[,2]
    x <- x[,1]
  }
  resolution <- round(resolution)
  x <- formatC(floor(x / resolution )*resolution, format="d")
  y <- formatC(floor(y / resolution )*resolution, format="d")
  resultat <- if(resinstr)
    paste0("r", resolution, "N", y, "E", x)
  else
    paste0("N", y, "E", x)

  nas <- which(is.na(y)|is.na(x))
  if (length(nas)>0) resultat[nas] <- NA

  resultat
}

#' Produce a stars object with zero value on a fixed resolution
#'
#' @param data spatial feature of the relevant locations
#' @param resolution cell size
#' @param values values to populate the stars values with, default to 0
#'
#' @export
stars_ref <- function(data, resolution, values = 0L) {

  stars::st_as_stars(sf::st_bbox(data), dx = resolution, dy = resolution, values = values)
}

#' Creates identifiers for (x,y) coordinates.
#' Used to work with crs 3035
#'
#' @param x matrix of coordinates (first column = x, second = y). Or a vector of x coordinates.
#' @param y default to NULL (if x is a matrix). Or a vector of y coordinates.
#' @param resolution default to 200. Size of raster cells which are identified.
#'
#' @export
coord2idINS <- function(x, y = NULL, resolution = 200) {

  if (!is.null(y)) x <- cbind(x, y)
  resolution <- round(resolution)

  res <- apply(x, MARGIN = 1:2, FUN = function(z) formatC(floor(z / resolution) * resolution, format = "d"))
  res <- paste0("r", resolution, "N", res[, 2], "E", res[, 1])

  nas <- apply(x, MARGIN = 1, FUN = anyNA)
  res[nas] <- NA

  res
}

#' Retrieves (x, y) coordinates from idINS.
#' Used to work with crs 3035
#'
#' @param ids character vector of idINS
#' @param stop_if_res_not_cst default to TRUE. Refuse to convert id with different resolutions
#'
#' @export
idINS2coord <- function(ids, stop_if_res_not_cst = TRUE, resolution = 200) {

  cr_pos <- stringr::str_locate(ids[[1]], "(?<=r)[0-9]+|(?<=RES)[0-9]+")
  cy_pos <- stringr::str_locate(ids[[1]], "(?<=N)[0-9]+")
  cx_pos <- stringr::str_locate(ids[[1]], "(?<=E)[0-9]+")
  y <- as.numeric(stringr::str_sub(ids,cy_pos[[1]],cy_pos[[2]]))
  x <- as.numeric(stringr::str_sub(ids,cx_pos[[1]],cx_pos[[2]]))
  r <- if(is.null(resolution))
    as.numeric(stringr::str_sub(ids[[1]],cr_pos[[1]], cr_pos[[2]]))
  else
    resolution
  res <- matrix(c(x + resolution/2, y + resolution/2), ncol = 2)
  colnames(res) <- c("x", "y")

  res
}

#' Retrieves (lon, lat) coordinates from idINS.
#'
#' @param idINS character vector of idINS
#' @param resolution default to NULL. Set to convert id with different resolutions
#'
#' @export
idINS2lonlat <- function(idINS, resolution=NULL) {
  if(length(idINS)==0)
    return(tibble::tibble(lon=numeric(), lat=numeric()))
  cr_pos <- stringr::str_locate(idINS[[1]], "(?<=r)[0-9]+|(?<=RES)[0-9]+")
  cy_pos <- stringr::str_locate(idINS[[1]], "(?<=N)[0-9]+")
  cx_pos <- stringr::str_locate(idINS[[1]], "(?<=E)[0-9]+")
  y <- as.numeric(stringr::str_sub(idINS,cy_pos[[1]],cy_pos[[2]]))
  x <- as.numeric(stringr::str_sub(idINS,cx_pos[[1]],cx_pos[[2]]))
  r <- if(is.null(resolution))
    as.numeric(stringr::str_sub(idINS[[1]],cr_pos[[1]], cr_pos[[2]]))
  else
    resolution
  x <- x+r/2
  y <- y+r/2
  lonlat <- sf_project(
    from=st_crs(3035),
    to=st_crs(4326),
    pts = matrix(c(x,y), nrow=length(x), ncol=2))
  return(
    tibble(lon = lonlat[, 1], lat = lonlat[, 2]))
}

#' Retrieves (lon, lat) coordinates from long idINS.
#'
#' @param idINS character vector of idINS
#' @param resolution default to NULL. Set to convert id with different resolutions
#'
#' @export
lidINS2lonlat <- function(idINS, resolution=NULL) {
  if(length(idINS)==0)
    return(tibble::tibble(lon=numeric(), lat=numeric()))
  cr_pos <- stringr::str_locate(idINS[[1]], "(?<=r)[0-9]+|(?<=RES)[0-9]+")
  cy_pos <- stringr::str_locate(idINS[[1]], "(?<=N)[0-9]+")
  cx_pos <- stringr::str_locate(idINS[[1]], "(?<=E)[0-9]+")
  y <- as.numeric(stringr::str_sub(idINS,cy_pos[[1]],cy_pos[[2]]))
  x <- as.numeric(stringr::str_sub(idINS,cx_pos[[1]],cx_pos[[2]]))
  r <- if(is.null(resolution))
    as.numeric(stringr::str_sub(idINS[[1]],cr_pos[[1]], cr_pos[[2]]))
  else
    resolution
  x <- x+r/2
  y <- y+r/2
  lonlat <- sf_project(
    from=st_crs(3035),
    to=st_crs(4326),
    pts = matrix(c(x,y), nrow=length(x), ncol=2))
  return(
    tibble(lon = lonlat[, 1], lat = lonlat[, 2]))
}


#' Calculate euclidean distance between to long idINS, in meter
#'
#' @param fromidINS character vector of starting idINS
#' @param toidINS character vector of ending idINS
#' @param resolution default to NULL. Set if no resolution is provided in idINS
#'
#' @export
#'

lidINS2dist <- function(fromidINS, toidINS, resolution=NULL) {
  stopifnot(length(fromidINS)==length(toidINS))
  if(length(fromidINS)==0)
    return(numeric())

  cr_pos <- stringr::str_locate(fromidINS[[1]], "(?<=r)[0-9]+|(?<=RES)[0-9]+")
  cy_pos <- stringr::str_locate(fromidINS[[1]], "(?<=N)[0-9]+")
  cx_pos <- stringr::str_locate(fromidINS[[1]], "(?<=E)[0-9]+")
  r <- if(is.null(resolution))
    as.numeric(stringr::str_sub(fromidINS[[1]],cr_pos[[1]], cr_pos[[2]]))
  else
    resolution
  fromy <- as.numeric(stringr::str_sub(fromidINS,cy_pos[[1]],cy_pos[[2]]))+r/2
  fromx <- as.numeric(stringr::str_sub(fromidINS,cx_pos[[1]],cx_pos[[2]]))+r/2

  cr_pos <- stringr::str_locate(toidINS[[1]], "(?<=r)[0-9]+|(?<=RES)[0-9]+")
  cy_pos <- stringr::str_locate(toidINS[[1]], "(?<=N)[0-9]+")
  cx_pos <- stringr::str_locate(toidINS[[1]], "(?<=E)[0-9]+")
  r <- if(is.null(resolution))
    as.numeric(stringr::str_sub(toidINS[[1]],cr_pos[[1]], cr_pos[[2]]))
  else
    resolution
  toy <- as.numeric(stringr::str_sub(toidINS,cy_pos[[1]],cy_pos[[2]]))+r/2
  tox <- as.numeric(stringr::str_sub(toidINS,cx_pos[[1]],cx_pos[[2]]))+r/2

  return(sqrt((tox-fromx)^2 + (toy-fromy)^2))
}

#' retourne un sf à partir d'un data.frame avec un champ idINS long
#'
#' @param data le data frame
#' @param resolution 200m par défaut
#' @param idINS le champ qui contient l'idINS court (défaut "idINS")
#'
#' @return un sf
#' @export
#'
lidINS2sf <- function(data, resolution = 200, idINS = "idINS") {
  data$geometry <- data[[idINS]] |>
    lidINS2square()
  sf::st_as_sf(data)
}
