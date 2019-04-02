blowup <-
function(x, M, number = 0) {
    y <- array(number, dim = c(M, M))
    nr <- nrow(x)
    nc <- ncol(x)
    dr <- (M - nr) %/% 2
    dc <- (M - nc) %/% 2
    y [(dr+1):(dr+nr), (dc+1):(dc+nc)] <- x
    return(y)
}
