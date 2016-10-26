\name{plot.vhica}
\alias{plot.vhica}
\encoding{utf8}

\title{
	Plots a VHICA regression between two species. 
}

\description{
	The VHICA method is based on a contrast between gene divergence and codon usage bias. A regression between divergence and codon usage provides a reference, and sequences of interest (typically, transposable elements) will be compared to the reference genes. 
}

\usage{
\method{plot}{vhica}(x, sp1 = NULL, sp2 = NULL, ...)
}
\arguments{
  \item{x}{
	An object of class \code{vhica}, created by \code{\link{read.vhica}}. 
}
  \item{sp1}{
	Name of the first species, as in the data files. 
}
  \item{sp2}{
	Name of the second species, as in the data files. 
	}
  \item{\dots}{
	Additional options for \code{\link{plot}}. 
}
}

\details{
	The resulting figure displays genes as circles, and transposable elements as symbols. 
}

%% \value{}

\references{
Gabriel Luz Wallau, Arnaud Le Rouzic, Pierre Capy, Elgion Loreto, Aurélie Hua-Van. VHICA: A new method to discriminate between vertical and horizontal transposon transfer: application to the mariner family within Drosophila. Molecular biology and evolution 33 (4), 1094-1109. 
}

\author{
Implementation: Arnaud Le Rouzic <lerouzic@legs.cnrs-gif.fr>\cr
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
image(vc, "mellifera:6", treefile=file.tree, skip.void=TRUE)
}


\keyword{regression}
\keyword{hplot}
