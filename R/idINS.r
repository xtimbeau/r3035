idINS2square <- function(ids, resolution=NULL)
{
  cr_pos <- stringr::str_locate(ids[[1]], "r(?=[0-9])")[,"start"]+1
  cy_pos <- stringr::str_locate(ids[[1]], "N(?=[0-9])")[,"start"]+1
  cx_pos <- stringr::str_locate(ids[[1]], "E(?=[0-9])")[,"start"]+1
  lcoord <- cx_pos-cy_pos-1
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

#' Récupère les coordonnées X et Y de idINS.
#'
#' @param ids vecteur d'idINS.
#' @param resolution resolution, par défaut celle attachée à idINS.
#'
#' @export
idINS2point <- function(ids, resolution=NULL)
{
  cr_pos <- stringr::str_locate(ids[[1]], "r(?=[0-9])")[,"start"]+1
  cy_pos <- stringr::str_locate(ids[[1]], "N(?=[0-9])")[,"start"]+1
  cx_pos <- stringr::str_locate(ids[[1]], "E(?=[0-9])")[,"start"]+1
  lcoord <- cx_pos-cy_pos-1
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
#' @param resolution resolution, par défaut 200.
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
idINS2coord <- function(ids, stop_if_res_not_cst = TRUE) {

  ids <- strsplit(as.character(ids), "[rNE]")
  resolution <- lapply(ids, FUN = function(liste) liste[2]) |> unlist() |> as.numeric()
  if (stop_if_res_not_cst) stopifnot(length(unique(resolution)) == 1L)

  x <- lapply(ids, FUN = function(liste) liste[4]) |> unlist() |> as.numeric()
  y <- lapply(ids, FUN = function(liste) liste[3]) |> unlist() |> as.numeric()

  res <- matrix(c(x + resolution/2, y + resolution/2), ncol = 2)
  colnames(res) <- c("x", "y")

  res
}
