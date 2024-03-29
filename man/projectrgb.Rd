% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{projectrgb}
\alias{projectrgb}
\title{Projette un raster RGB}
\usage{
projectrgb(rgb, crs = "3035", res = 200)
}
\arguments{
\item{rgb}{le raster en entrée}

\item{crs}{le système de coordonnées en sortie}
}
\value{
un raster RGB projetté dans
}
\description{
En procédant à une interpolation, la projection simple d'un raster
produit un raster avec des valeurs non entières même si l'input
est composé de composantes RGB et entière (0-255)
La fonction projette donc dans un système de coordonnées
et rebase et arrondi les différentes couches.
}
