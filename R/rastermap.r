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
#' @param resolution résolution du raster (par défaut cherche la meilleure résolution disponible, sinon, caclule res_def sur le sf)
#' @param res_def résolution par défaut si aucun idINS n'est fournit
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
           palette=colorspace::sequential_hcl(n=n, palette="Terrain 2", rev=rev),
           style="kmeans",
           resolution=Inf,
           res_def = 200,
           decor=NULL,
           bbox=NULL, ...) {

    quo_var <- rlang::enquo(var)

    raster.temp <- rastervar(data=data, var={{var}}, fun=fun, dropth=dropth, resolution=resolution, res_def=res_def)

    text.temp <-ifelse(!is.null(label), label, rlang::quo_name(quo_var))
    decor$fdc+
      tmap::tm_shape(raster.temp, bbox=bbox) +
      tmap::tm_raster(title = stringr::str_c("Dens. ", text.temp), style=style, palette = palette, n=n, ...)+
      decor$hdc
  }

#' ratervar
#'
#' @param data données (tibble, data.frame) avec un idINS ou un sf (l'idINS sera fabriqué par projection sur une grille 3035)
#' @param ... expressions en quasiquotation
#' @param fun fonction d'agrégagtion
#' @param dropth Elimine les extrêmes de la distribtion (un double, inférieur à 1)
#' @param resolution résolution par défaut, cherche la plus fine, sinon, calcule à res_def
#' @param idINS nom de la variable idINS
#' @param res_def résolution par défaut si pas de idINS fourni
#'
#' @return
#' @import data.table
#' @export
rastervar <-
  function(data, ...,
           fun=mean, # opérateur d'aggrégation
           dropth = 0, # drop 1% des valeurs extrêmes
           resolution=Inf, idINS="idINS", res_def = 200) {
    quo_var <- rlang::enquos(...)
    findres <- getINSres(data,resolution=resolution,idINS=idINS)
    idinspire <- findres$idINS
    resolution <- findres$res
    if (any(idinspire==FALSE))
      idinspire <- idINS3035(
        sf::st_coordinates(
          sf::st_centroid(data %>% sf::st_as_sf())),
        resolution=resolution)
    else
      idinspire <- data[[idinspire]]

    data.temp <- purrr::map_dfc(quo_var, ~{data %>%
        dplyr::as_tibble() %>%
        dplyr::transmute(!!rlang::quo_name(.x) := !!.x)})

    data.table::setDT(data.temp)

    vars <- rlang::set_names(names(data.temp))
    isnum <- purrr::map_lgl(data.temp, is.numeric)
    data.temp <- data.temp[, lapply(.SD, factor2num)]
    # for (j in names(data.temp)) set(data.temp, j = j, value = factor2num(data.temp[[j]]))
    obs_na <- purrr::map(data.temp, ~is.na(.x))
    obs_drop <- purrr::reduce(obs_na, `&`)
    set(data.temp, j="idINS", value=idinspire)
    data.temp <- data.temp[!obs_drop,]

    if (dropth>0)
      for(v in names(isnum))
      {
        data.temp[, (v):=ifelse(selxth(get(v), dropth), get(v), NA_real_)]
      }
    data.temp <- data.temp[, lapply(.SD, function(x) fun(x, na.rm=TRUE)), by=idINS]

    dt2r(data.temp, resolution=resolution, idINS="idINS")
  }
