\name{vhica_0_1-package}
\alias{vhica_0_1-package}
\alias{vhica_0_1}
\docType{package}
\encoding{utf8}


\title{
Vertical and Horizontal Inheritance Consistence Analysis
}

\description{
The package implements the VHICA method described in Wallau et al. (in prep). The purpose of the method is to detect horizontal transfers of transposable elements, by contrasting the divergence of transposable element sequences with that of regular genes. Two files should be provided, for both a set of reference genes and transposable element sequences: (i) pairwise divergence across species (e.g., dS), (ii) codon usage bias for all genes and elements in all species.
}

\details{
\tabular{ll}{
Package: \tab vhica_0_1\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-11\cr
License: \tab GPL-v2\cr
}
This package contains three main functions.
\itemize{
	\item {\code{\link{read.vhica}}: } {reads two data files (divergence and codon usage) and generates an object of class \code{vhica} that will be used for further analysis.}
	\item {\code{\link{plot.vhica}}: } {plots the VHICA regression between two species, and displays how far transposable elements (or any kind of other sequences) are from the reference genes.}
	\item {\code{\link{image.vhica}}: } {plots the consistency of a specific element across all species, which makes it possible to build evolutionary scenarios.}
}
}

\author{
Implementation: Arnaud Le Rouzic <lerouzic@legs.cnrs-gif.fr> \cr
Scientists who designed the method: Gabriel Wallau, Aurélie Hua-Van, Arnaud Le~Rouzic.

Maintainer: Arnaud Le Rouzic <lerouzic@legs.cnrs-gif.fr>
}
\references{
Gabriel Luz Wallau, Arnaud Le Rouzic, Pierre Capy, Elgion Loreto, Aurélie Hua-Van. VHICA: A new method to discriminate between vertical and horizontal transposon transfer: application to the mariner family within Drosophila. in prep. 
}

\keyword{ package }
% \seealso{ }
\examples{
file.cb <- system.file("extdata", "mini-cbias.txt", package="vhica")
file.div <- system.file("extdata", "mini-div.txt", package="vhica")
file.tree <- if(require("ape")) system.file("extdata", "phylo.nwk", package="vhica") else NULL
vc <- read.vhica(cb.filename=file.cb, div.filename=file.div)
plot(vc, "dere", "dana")
image(vc, "mellifera:6", treefile=file.tree, skip.void=TRUE)
}
