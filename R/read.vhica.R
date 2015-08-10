read.vhica <-
function (gene.fasta=NULL, target.fasta=NULL, cb.filename=NULL, div.filename=NULL, 
	reference = "Gene", divergence = "dS", CUB.method="ENC", div.method="LWL85", div.pairwise=FALSE, species.sep="_", gene.sep=".", ...) 
{
	stopifnot( 
		!(is.null(gene.fasta) && is.null(target.fasta)) || 
		!(is.null(cb.filename) && is.null(div.filename)))
    vhica.obj <- list()
    if (!is.null(gene.fasta)) {
		vhica.obj$cbias <- 
			.seq.codon.bias(gene.fasta=gene.fasta, target.fasta=target.fasta, method=CUB.method, species.sep=species.sep)
		vhica.obj$div <- 
			.seq.divergence(sequence.fasta=c(gene.fasta, target.fasta), method=div.method, pairwise=div.pairwise, species.sep=species.sep)
		if (!is.null(cb.filename))
			write.table(vhica.obj$cbias, file=cb.filename, sep="\t", quote=FALSE)
		if (!is.null(div.filename)) 
			write.table(vhica.obj$div, file=div.filename, sep="\t", quote=FALSE)
    } else {
		vhica.obj$cbias <- .read.codon.bias(file = cb.filename, reference = reference)
		vhica.obj$div <- .read.divergence(file = div.filename, divergence = divergence)    
    }
    vhica.obj$reg <- .reference.regression(vhica.obj$cbias, vhica.obj$div, 
        reference = reference, divergence = divergence, ...)
    vhica.obj$reference <- reference
    tmp.target <- levels(vhica.obj$cbias[, "Type"])
    vhica.obj$target <- tmp.target[tmp.target != reference][1]
    vhica.obj$divergence <- divergence
    class(vhica.obj) <- c("vhica", class(vhica.obj))
    return(vhica.obj)
}
