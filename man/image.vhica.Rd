\name{image.vhica}
\alias{image.vhica}
\encoding{utf8}

\title{
Consistency matrix for a transposable element in the VHICA analysis.  
}
\description{
This function plots a composite figure summarizing the evolutionary properties of a transposable element in a group of related species. Discrepancies may indicate horizontal transfers. 
}
\usage{
\method{image}{vhica}(x, element = "", H1.test = "bilat", treefile = NULL, 
skip.void = FALSE, species = NULL, p.threshold = 0.05, 
p.adjust.method = "bonferroni",  ncolors = 1024, 
main = element, threshcol=0.1, colsqueeze=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	An object of class \code{vhica}, created by the function \code{\link{read.vhica}}.
}
  \item{element}{
	The name of the transposable element, as specified in the data files. If the element is not present in the data, the program halts. 
}
  \item{H1.test}{
	A value among \code{"bilat"}, \code{"lower"}, or \code{"greater"}. 
}
  \item{treefile}{
	A Newick file containing a phylogenetic tree. Species names in the tree need to match the data. If absent, the figure will not display the phylogenetic relationship (which makes the interpretation impossible). 
}
  \item{skip.void}{
	Whether or not the figure should show species that do not contain the transposable element. 
}
  \item{species}{
	A named character vector to display pretty species names. The names of the vector are the real species names (as they will appear in the figure), the content of the vector are species codes as in the data files. 
}
  \item{p.threshold}{
	Threshold for the p-value (above which the color gradient increases). 
}
  \item{p.adjust.method}{
	As documented in \code{\link{p.adjust}}. 
}
  \item{ncolors}{
	Number of colors in the gradient. 
}
  \item{main}{
	Main title of the figure (default: the name of the transposable element). 
}
  \item{threshcol}{
	Part of the color spectrum devoted to non-significant values.
}
 \item{colsqueeze}{
	Values larger than 1 shrink the color gradient around the threshold.
}
  \item{\dots}{
	Further arguments to the generic function \code{\link{image}}. 
}
}
\details{
	The figure displays in blue TE copies that are more divergent than expected between species, and in red copies that are less divergent than expected. If several lineages of copies are present in a species, the table will be split to display both lineages. Keys for the interpretation of the pattern and the reconstruction of an evolutionary scenario are provided in the original publication. 
}
\value{
	The function returns (invisibly) a list of elements (object of class \code{vhicaimage}) which can be used for further analysis: \code{tree} contains the phylogenetic tree (object of class \code{"phylo"}), \code{species} is the vector of species, \code{stats} is a matrix of log10(P-values) (positive elements are minus log10(P) corresponding to positive residuals), code{dS} is a matrix reminding the divergence rates from the data. Calling the method \code{\link{summary.vhicaimage}} on this object returns a nicely formatted data frame. 
}
\references{
Gabriel Luz Wallau, Arnaud Le Rouzic, Pierre Capy, Elgion Loreto, Aurélie Hua-Van. VHICA: A new method to discriminate between vertical and horizontal transposon transfer: application to the mariner family within Drosophila. Molecular biology and evolution 33 (4), 1094-1109.
}
\author{
Implementation: Arnaud Le Rouzic <lerouzic@legs.cnrs-gif.fr> \cr
Scientists who designed the method: Gabriel Wallau, Aurélie Hua-Van, Arnaud Le~Rouzic.
}
%% \note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.vhica}}, \code{\link{plot.vhica}}, \code{\link{summary.vhicaimage}}.
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


\keyword{hplot}
\keyword{regression}
\keyword{tree}
