% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rastermap.r
\name{rastermap}
\alias{rastermap}
\title{Affiche une carte comme un raster}
\usage{
rastermap(
  data,
  var,
  label = NULL,
  fun = mean,
  dropth = 0,
  n = 5,
  rev = FALSE,
  palette = colorspace::sequential_hcl(n = n, palette = "Terrain 2", rev = rev),
  style = "kmeans",
  title = NULL,
  resolution = Inf,
  res_def = 200,
  decor = NULL,
  bbox = NULL,
  ...
)
}
\arguments{
\item{data}{données (tibble, data.frame) avec un idINS ou un sf (l'idINS sera fabriqué par projection sur une grille 3035)}

\item{var}{l'expression à afficher sur la carte, évaluée selon le schéma tidyverse}

\item{label}{le nom de la variable construite affichée sur la carte}

\item{fun}{une fonction pour aggréger si nécessaire (mean par défaut)}

\item{dropth}{élemine de l'échantillon les valeurs extrêmes}

\item{n}{nombre de groupe pour la carte}

\item{rev}{palette à l'endroit ou à l'envers (FALSE par défaut)}

\item{palette}{une palette de couleur (terrain 2 par défaut)}

\item{style}{passé à tmap, style de la carte}

\item{title}{titre de la légende, si NULL construit à partir de la variable}

\item{resolution}{résolution du raster (par défaut cherche la meilleure résolution disponible, sinon, caclule res_def sur le sf)}

\item{res_def}{résolution par défaut si aucun idINS n'est fournit}

\item{decor}{sous couche de carte à ajouter (doit contenir un $fdc en dessous, et un $hdc au dessus)}

\item{bbox}{boite englobant la carte}

\item{...}{paramètres supplémentaires passés à tmap}
}
\value{
une carte
}
\description{
A partir d'un tableau (tibble, data.frame, data.table) équipé d'un champ idINS ou un sf,
on fabrique un raster pour la variable d'intérêt et on affiche une carte
}
