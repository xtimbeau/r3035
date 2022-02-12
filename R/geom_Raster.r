#' ajoute une couche raster à un ggplot
#'
#' geom_raster demande de préciser le système de coordonnées (aes(x=lon, y=lat)). \code{geom_Ratser} (avec un R majuscule)
#' permet de s'en passer. Attention le système de coordonnées est uniquement 3035.
#'
#' @param raster le raster (package raster) au format rasterlayer, rasterbrick ou rasterstack à afficher
#' @param mapping un mapping supplémentaire (définit la colonne employé pour la couleur)
#' @param ... autres paramètres passés à geom_tile
#' @param long permet de passer des raster au format long si on veut pour les différentes couches du raster
#' @param style kmeans
#' @param k nombre de groupes (kmeans)
#'
#'
#' @import data.table
#' @return une couche ggplot2
#' @export
#'
geom_Raster <- function(raster, mapping=NULL, ..., long=FALSE, style="cont", k = 5)
{
  checkmate::assert(checkmate::checkMultiClass(raster, c("RasterLayer", "RasterBrick", "RasterStack")))
  coords <- raster::coordinates(raster)
  data <- setDT(as.data.frame(raster))
  vv <- names(raster)
  data[, `:=`(x=coords[,1],y=coords[,2])]
  nas <- reduce(vv, function(x,v) x & is.na(data[[v]]), .init=TRUE)
  data <- data[!nas]
  if(long)
  {
    data <- melt(data, measure.vars=vv, na.rm=TRUE)
    vv <- "value"
  }
  if(style=="kmeans")
    for(v in vv)
    {
      clus <- kbins(data[[v]], k=k, bins="factor")
      data[, (v):= clus]
    }
  list(ggplot2::geom_tile(data=data, mapping=modifyList(aes(x=x, y=y), mapping), ...), ggplot2::coord_sf(crs=st_crs(3035)))
}


# utilisée dans geom_Ratser, interne
#
kbins <- function(x, k=5, bins="factor")
{
  kmeans <- kmeans(x, centers=k, nstart=1L, iter.max = 1000, algorithm = "Lloyd")
  cc <- tibble::tibble(centers = as.vector(kmeans$centers), i = 1:nrow(kmeans$centers))
  cc <- cc  |>
    dplyr::arrange(cc, centers)  |>
    dplyr::mutate(ii=row_number())  |>  arrange(i)  |>  pull(ii)
  res <- tibble::tibble( x = x, cluster=cc[kmeans$cluster])
  minmax <- res  |>
    dplyr::group_by(cluster)  |>
    dplyr::summarize(max=max(x), min=min(x), mean=mean(x), median=median(x))
  if(bins=="factor")
    return(factor(res$cluster, labels=str_c("[", uf2si2(minmax$min), ", ", uf2si2(minmax$max), "]")))
  if(bins=="mean")
    return(minmax$mean[res$cluster])
  if(bins=="median")
    return(minmax$median[res$cluster])
  if(bins=="cluster")
    return(res$cluster)
  if(bins=="fullcluster")
    return(list(x=res$cluster, labels=str_c("[", uf2si2(minmax$min), ", ", uf2si2(minmax$max), "]")))
  return(rep(NA, length(x)))
}
