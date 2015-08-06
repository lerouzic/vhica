\name{read.vhica}
\alias{read.vhica}
\encoding{utf8}

\title{
	Reads divergence and codon usage data files for the VHICA method. 
}

\description{
	The VHICA method relies on two sources of information: (i) the divergence between sequences, and (ii) the codon usage bias. This function reads two data files and creates an object of class \code{vhica} that can be further explored by \code{\link{plot.vhica}} and \code{\link{image.vhica}}.
}

\usage{
read.vhica(cb.filename, div.filename, reference = "Gene", 
divergence = "dS", ...)
}

\arguments{
  \item{cb.filename}{
	File name for the codon usage bias data. 
}
  \item{div.filename}{
	File name for the divergence data. 
}
  \item{reference}{
	Name of the reference type in the codon usage file. Default is "Gene". 
}
  \item{divergence}{
	Name of the divergence column in the divergence file. Default is "dS". 
}
  \item{\dots}{
	Further parameters for the internal function \code{.reference.regression}. 
}
}

\details{
	\itemize{
		\item{Codon usage bias}{
			Example of data file:
			
\preformatted{
        Type    sp1     sp2     sp3
CG4231  Gene    42.3    51.1    47.2
CG2214  Gene    47.2    44.9    53.2
Pelem1  TE      36.2    47.0    44.4
...}

			\itemize{
				\item {Row names (or first column)}{sequence index}
				\item {Type} {whether the sequence is a reference (default: Gene) or a focal sequence (transposable element, ...)}
				\item {Following columns} {a measurement of codon bias (ENC, CBI...) for every species}
			}
		}
		\item{Divergence}{
			Example of data file:
			
\preformatted{
seq     dS      sp1     sp2
CG4231  0.84    Dmel    Dsim
CG4231  0.46    Dmel    Dana
CG4231  0.58    Dsim    Dana
CG2214  0.10    Dmel    Dsim
...}

			\itemize{
				\item {First column (or row names):} {sequence index}
				\item {Second column:} {divergence measurement}
				\item {Columns 3 and 4:} {the pair of species on which the divergence is calculated}
				\item {}{Row names and Col names are allowed but disregarded}
			}
		}
	}
}

\value{
	The function returns an object of class \code{vhica}, a list containing:
	\itemize{
		\item{cbias:} {A codon bias array}
		\item{div:} {The divergence matrix}
		\item{reg:} {The result of all pairwise regressions}
		\item{reference:} {The \code{reference} option}
		\item{target:} {The sequence type that is not the reference}
		\item{divergence:} {The \code{divergence} option}
	}
}

\references{
Gabriel Luz Wallau, Arnaud Le Rouzic, Pierre Capy, Elgion Loreto, Aurélie Hua-Van. VHICA: A new method to discriminate between vertical and horizontal transposon transfer: application to the mariner family within Drosophila. in prep. 
}

\author{
Implementation: Arnaud Le Rouzic <lerouzic@legs.cnrs-gif.fr> \cr
Scientists who designed the method: Gabriel Wallau, Aurélie Hua-Van, Arnaud Le~Rouzic.
}

%% \note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.vhica}}, \code{\link{image.vhica}}
}
\examples{
file.cb <- system.file("extdata", "mini-cbias.txt", package="vhica")
file.div <- system.file("extdata", "mini-div.txt", package="vhica")
file.tree <- if(require("ape")) system.file("extdata", "phylo.nwk", package="vhica") else NULL
vc <- read.vhica(cb.filename=file.cb, div.filename=file.div)
plot(vc, "dere", "dana")
image(vc, "mellifera:6", treefile=file.tree, skip.void=TRUE)
}

\keyword{ data }
