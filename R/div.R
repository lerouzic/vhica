div <- function
(file=NULL, sequence=NULL, sqs=NULL, method="LWL85", pairwise=FALSE)
{
	stopifnot(
		!(is.null(file) && is.null(sequence)),
		method[1] %in% c("LWL85"),
		requireNamespace("gtools", quietly=TRUE))

	if (!is.null(file)) {
		if (!requireNamespace("seqinr", quietly=TRUE)) {
			stop("Reading FASTA files require package seqinr")
		}		
		sequence <- seqinr::read.fasta(file)
	}
	sequence <- .checkseq(sequence)
	if (is.null(sqs)) {
		sqs <- names(sequence)
	}
	combn <- gtools::combinations(n=length(sqs), r=2, v=sqs)		
	if (method[1]=="LWL85") {
		return(data.frame(sq1=combn[,1], sq2=combn[,2], div=.LWL85(sequence, combn[,1], combn[,2], pairwise=pairwise)))
	}	
}
