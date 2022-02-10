#' Affiche une carte comme un raster
#'
#' A partir d'un tableau (tibble, data.frame, data.table) équipé d'un champ idINS ou un sf,
#' on fabrique un raster pour la variable d'intérêt et on affiche une carte
#'
#' @param data données (tibble, data.frame) avec un idINS ou un sf (l'idINS sera fabriqué par projection sur une grille 3035)
#' @param var l'expression à afficher sur la carte, évaluée selon le schéma tidyverse
#' @param label le nom de la variable construite affichée sur la carte
#' @param fun une fonction pour aggréger si nécessaire (mean par défaut)
#' @param dropth élemine de l'échantillon les valeurs extrêmes
#' @param n nombre de groupe pour la carte
#' @param rev palette à l'endroit ou à l'envers (FALSE par défaut)
#' @param palette une palette de couleur (terrain 2 par défaut)
#' @param style passé à tmap, style de la carte
#' @param resolution résolution du raster (200m par défaut)
#' @param decor sous couche de carte à ajouter (doit contenir un $fdc en dessous, et un $hdc au dessus)
#' @param bbox boite englobant la carte
#' @param ... paramètres supplémentaires passés à tmap
#'
#' @return
#' @export
#'
rastermap <-
  function(data, var,
           label=NULL, # pour le graphe
           fun=mean, # opérateur d'aggrégation
           dropth = 0, # drop 1% des valeurs extrêmes
           n=5,
           rev=FALSE,
           palette=sequential_hcl(n=n, palette="Terrain 2", rev=rev),
           style="kmeans",
           resolution=50,
           decor=NULL,
           bbox=NULL, ...) {
    library("tmap", quietly=TRUE)

    quo_var <- rlang::enquo(var)

    raster.temp <- rastervar(data=data, var={{var}}, fun=fun, dropth=dropth, resolution=resolution)

    text.temp <-ifelse(!is.null(label), label, rlang::quo_name(quo_var))
    decor$fdc+
      tm_shape(raster.temp, bbox=bbox) +
      tm_raster(title = str_c("Dens. ", text.temp), style=style, palette = palette, n=n, ...)+
      decor$hdc
  }

# fabrique la variable pour le raster, interne
rastervar <-
  function(data, ...,
           fun=mean, # opérateur d'aggrégation
           dropth = 0, # drop 1% des valeurs extrêmes
           resolution=50, idINS="idINS") {

    library("data.table", quietly=TRUE)

    quo_var <- rlang::enquos(...)
    idinspire <- getINSres(data,resolution=resolution,idINS=idINS)
    if (any(idinspire==FALSE))
      idinspire <- idINS3035(st_coordinates(st_centroid(data %>% st_as_sf)), resolution=resolution)
    else
      idinspire <- data[[idinspire]]

    data.temp <- purrr::map_dfc(quo_var, ~{data %>%
        dplyr::as_tibble() %>%
        dplyr::transmute(!!rlang::quo_name(.x) := !!.x)})

    data.table::setDT(data.temp)
    vars <- rlang::set_names(names(data.temp))
    isnum <- purrr::map_lgl(data.temp, is.numeric)
    data.temp <- data.temp[, lapply(.SD, factor2num)]
    obs_na <- purrr::map(data.temp, ~is.na(.x))
    obs_drop <- purrr::reduce(obs_na, `&`)
    data.temp[, idINS := idinspire]
    data.temp <- data.temp[!obs_drop]

    if (dropth>0)
      for(v in vars(isnum))
        data.temp[, (v):=ifelse(selxth(get(v), dropth), get(v), NA_real_)]

    data.temp <- data.temp[, lapply(.SD, function(x) fun(x, na.rm=TRUE)), by=idINS]

    dt2r(data.temp, resolution=resolution, idINS="idINS")
  }
