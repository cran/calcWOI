buildperiodic <-
function(x) {
    x1 <- x [, dim(x)[2]:1, drop = FALSE]
    x2 <- cbind(x, x1)
    x3 <- x2 [dim(x)[2]:1, , drop = FALSE]
    outarray <- rbind(x2, x3)
    return(outarray)
}
