summary.vhicaimage <- function(object, divrate=NA, p.thresh=1, ...)
{
    st <- object$stats
    ds <- object$dS
    st[lower.tri(st)] <- NA
    ds[lower.tri(ds)] <- NA
    tokeep <- c(!is.na(ds)) # should be the same as !is.na(st)
    ans <- data.frame(
        expand.grid(colnames(object$stats), rownames(object$stats))[tokeep,], 
        'p.value'= c(10^object$stats)[tokeep], 
        dS=     c(object$dS)[tokeep])
    colnames(ans)[1:2] <- c("sp1","sp2")
    rownames(ans) <- NULL
    if (!is.na(divrate)) 
        ans[,"Time(Mya)"] <- ans$dS/(2*divrate)
    return(ans[ans$'p.value' <= p.thresh,])
}
