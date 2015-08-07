CUB <-
function(file=NULL, sequence=NULL, method="ENC")
{
	stopifnot(
		!(is.null(file) && is.null(sequence)),
		method[1] %in% c("ENC"))

	if (!is.null(file)) {
		if (!requireNamespace("seqinr", quietly=TRUE)) {
			stop("Reading FASTA files require package seqinr")
		}		
		sequence <- seqinr::read.fasta(file)
	}
	sequence <- .checkseq(sequence)
	if (method[1]=="ENC") {
		return(sapply(sequence, .ENC))
	}
}
