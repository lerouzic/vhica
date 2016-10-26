\name{summary.vhicaimage}
\alias{summary.vhicaimage}
\encoding{utf8}

\title{
	Provides a data.frame that nicely displays the information returned by \code{\link{image.vhica}}. 
}

\description{
	The \code{\link{image.vhica}} routine displays visually the statistical support for horizontal transfers, and can return an object of class \code{vhicaimage}. The current \code{summary} method reorganizes this object into a dataframe that can be displayed or reused in further analysis. 
}

\usage{
\method{summary}{vhicaimage}(object, divrate=NA, p.thresh=1, ...)
}
\arguments{
  \item{object}{
	An object of class \code{vhicaimage}, created by \code{\link{image.vhica}}. 
}
  \item{divrate}{
	Optional divergence rate (in neutral substitutions per Myr). 
}
  \item{p.thresh}{
	Optional p-value threshold. By default, all data is returned. 
	}
  \item{\dots}{
	Additional options for \code{\link{summary}} (unused). 
}
}


\value{	The resulting \code{data.frame} has 4 or 5 columns. The two first columns are \code{sp1} and \code{sp2}, the two species between which the horizontal transfer is tested, in an arbitrary order. The column \code{p.value} contains the p-value calculated as in \code{\link{image.vhica}} (including the possible correction for multiple testing). The \code{dS} column is a copy of the corresponding divergence from the original data. The last, optional column \code{ Time(Mya)} is a molecular clock estimate of the time of divergence between the two sequences, based on the divergence rate (when provided). }

\references{
Gabriel Luz Wallau, Arnaud Le Rouzic, Pierre Capy, Elgion Loreto, Aurélie Hua-Van. VHICA: A new method to discriminate between vertical and horizontal transposon transfer: application to the mariner family within Drosophila. Molecular biology and evolution 33 (4), 1094-1109.
}

\author{
Implementation: Arnaud Le Rouzic <lerouzic@legs.cnrs-gif.fr> and Gabriel Wallau\cr
Scientists who designed the method: Gabriel Wallau, Aurélie Hua-Van, Arnaud Le~Rouzic.
}

%% \note{}

\seealso{
\code{\link{read.vhica}}, \code{\link{image.vhica}}
}

\examples{
file.cb <- system.file("extdata", "mini-cbias.txt", package="vhica")
file.div <- system.file("extdata", "mini-div.txt", package="vhica")
file.tree <- if(require("ape")) system.file("extdata", "phylo.nwk", package="vhica") else NULL
vc <- read.vhica(cb.filename=file.cb, div.filename=file.div)
plot(vc, "dere", "dana")
im <- image(vc, "mellifera:6", treefile=file.tree, skip.void=TRUE)
summary(im)
}


\keyword{print}
