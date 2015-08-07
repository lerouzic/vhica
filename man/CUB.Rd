\name{CUB}
\alias{CUB}
\title{
	Computes the Codon Usage Bias of DNA sequences
}
\description{
	The function reads aligned sequences in a fasta file and estimates the codon usage bias for each sequence. Several methods exist to estimate CUB; so far, only the "Effective Number of Codons" (ENC) calculation is available. 
}
\usage{
CUB(file = NULL, sequence = NULL, method = "ENC")
}
\arguments{
  \item{file}{FASTA file in which aligned sequences are stored.
}
  \item{sequence}{Alternatively, the result of seqinr::read.fasta.
}
  \item{method}{The method used to compute CUB. "ENC": Effective Number of Codons, as described in Wright (1990). 
}
}
\value{
A named vector of CUB scores. Names correspond to sequence names in the dataset. 
}
\references{
Wright, F. (1990). The 'effective number of codons' used in a gene. Gene, 87(1), 23-29.
}
\author{
Aurelie Hua-Van and Arnaud Le Rouzic.
}
\seealso{
	\code{\link{div}}
}
\examples{
	seq.file <- system.file("extdata/Genes", "Amd.fas", package="vhica")
	CUB(seq.file)
}
