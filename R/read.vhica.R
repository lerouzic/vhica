read.vhica <-
function (cb.filename, div.filename, reference = "Gene", divergence = "dS", 
    ...) 
{
    vhica.obj <- list()
    vhica.obj$cbias <- .read.codon.bias(file = cb.filename, reference = reference)
    vhica.obj$div <- .read.divergence(file = div.filename, divergence = divergence)
    vhica.obj$reg <- .reference.regression(vhica.obj$cbias, vhica.obj$div, 
        reference = reference, divergence = divergence, ...)
    vhica.obj$reference <- reference
    tmp.target <- levels(vhica.obj$cbias[, "Type"])
    vhica.obj$target <- tmp.target[tmp.target != reference][1]
    vhica.obj$divergence <- divergence
    class(vhica.obj) <- c("vhica", class(vhica.obj))
    return(vhica.obj)
}
