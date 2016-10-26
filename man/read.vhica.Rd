\name{read.vhica}
\alias{read.vhica}
\encoding{utf8}

\title{
	Reads divergence and codon usage data files for the VHICA method. 
}

\description{
	The VHICA method relies on two sources of information: (i) the divergence between sequences, and (ii) the codon usage bias. This function reads two data files and creates an object of class \code{vhica} that can be further explored by \code{\link{plot.vhica}} and \code{\link{image.vhica}}. Input can be either (1) two vectors of fasta file names (one for the genes, one for the putatively transfered genes), or (2) already processed files containing codon usage bias and divergence data (see Details). 
}

\usage{
read.vhica(gene.fasta=NULL, target.fasta=NULL, 
	cb.filename=NULL, div.filename=NULL, 
	reference = "Gene", divergence = "dS", 
	CUB.method="ENC", div.method="LWL85", div.pairwise=TRUE, 
	div.max.lim=3, species.sep="_", gene.sep=".", family.sep=".", ...)
}

\arguments{
  \item{gene.fasta}{
    Sequence files (FASTA format) containing the aligned sequences 
    (respecting the translation phase) for all species of the reference genes.}
  \item{target.fasta}{
    Sequence files (FASTA format) containing the aligned sequence 
    of the putatively transfered genes.}
  \item{cb.filename}{
	File name for the codon usage bias data. If FASTA files are provided, this file will be created.}
  \item{div.filename}{
	File name for the divergence data. If FASTA files are provided, this file will be created.}
  \item{reference}{
	Name of the reference type in the codon usage file. Default is "Gene".}
  \item{divergence}{
	Name of the divergence column in the divergence file. Default is "dS".}
  \item{CUB.method}{
    Method to be used for Codon Usage Bias calculation (see \code{\link{CUB}}).}
  \item{div.method}{
    Method to be used for divergence calculation (see \code{\link{div}}).}
  \item{div.pairwise}{
	Whether divergence should be calculated from the whole alignment of between pairs of sequences
	(see \code{\link{div}}).}
  \item{div.max.lim}{Maximum divergence score. Estimated divergence much larger than 100\% are likely 
	to be problematic and should not be considered.}
  \item{species.sep}{
	Separator for species (or equivalent) labels in sequence names. Any character string following 
	this separator will be disregarded -- be careful about potential duplicates.}
  \item{gene.sep}{
    Separator for gene names from gene sequence files.}
  \item{family.sep}{
	Separator for target sequence sub-families.}
  \item{\dots}{
	Further parameters for the internal function \code{.reference.regression}.}
}

\details{
	Details about CUB and divergence calculations can be found in \code{\link{CUB}} and \code{\link{div}}. If CUB and/or divergence need to be calculated by an external program, it is possible to provide them in the following format: 
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
		\item{family.sep:} {The character used to indicate TE sub-families}
	}
}

\references{
Gabriel Luz Wallau, Arnaud Le Rouzic, Pierre Capy, Elgion Loreto, Aurelie Hua-Van. VHICA: A new method to discriminate between vertical and horizontal transposon transfer: application to the mariner family within Drosophila. Molecular biology and evolution 33 (4), 1094-1109. 
}

\author{
Implementation: Arnaud Le Rouzic \cr
Scientists who designed the method: Gabriel Wallau, Aurelie Hua-Van, Arnaud Le Rouzic.
}

%% \note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.vhica}}, \code{\link{image.vhica}}, \code{\link{CUB}}, \code{\link{div}}
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
