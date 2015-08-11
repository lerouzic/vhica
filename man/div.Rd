\name{div}
\alias{div}
\title{
	Computation of the synonymous divergence between sequences
}
\description{
	The divergence between DNA sequences can be synonymous (neutral) or non-synonymous. Synonymous differences are generally considered as a better proxy for evolutionary divergence, as it is not affected by selection. This function computes the synonymous divergence between sequences. 
}
\usage{
div(file = NULL, sequence = NULL, sqs = NULL, method = "LWL85", 
	pairwise = TRUE, max.lim = 3)
}
\arguments{
  \item{file}{FASTA file in which aligned sequences are stored.
}
  \item{sequence}{Alternatively, the result of seqinr::read.fasta.
}
  \item{sqs}{Vector of sequence names to be compared. If not provided, all pairwise comparisons will be performed.
}
  \item{method}{Method used to compute the divergence. So far, only the LWL85 method (from Li et al. 1985).
}
  \item{pairwise}{Boolean: should the divergence be calculated for each pair of sequences or on the whole dataset? This is of particular importance when indels (gaps) are present in sequences, as codons with gaps are generally discarded by most methods. Setting this option to \code{TRUE} is thus more likely to give accurate results with multiple-gap sequences, but the calculation will also be slower.
}
  \item{max.lim}{Maximum value for divergence. Depending on the algorithm, various corrections can bring the divergence value above 100\%. Values larger than \code{max.lim} will be replaced by \code{NA}s, as they can be problematic for further statistical tests.}
}
\details{
	The \code{LWL85} method is a wrapper around the \code{kaks} function from the \code{seqinr} package. 
}
\value{
	A 3-column data frame with the following fields:
	\itemize{
		\item{div:} {The divergence score}	
		\item{sq1:} {The first sequence in the comparison}
		\item{sq2:} {The second sequence in the comparison}
	}
}
\references{
Li, W. H., Wu, C. I., & Luo, C. C. (1985). A new method for estimating synonymous and nonsynonymous rates of nucleotide substitution considering the relative likelihood of nucleotide and codon changes. Molecular biology and evolution, 2(2), 150-174.
}
\author{Aurelie Hua-Van and Arnaud Le Rouzic
}
\seealso{
	\code{\link{CUB}}
}
\examples{
	seq.file <- system.file("extdata/Genes", "Amd.fas", package="vhica")
	div(seq.file)
}

