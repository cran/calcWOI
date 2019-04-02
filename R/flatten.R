flatten <-
function(x, filter) {
    k <- length(filter)
    nr <- nrow(x)
    nc <- ncol(x)
    x [(1):(k), (1):(nc)] <- sweep(x [(1):(k), (1):(nc)], 1, filter, '*')
    x [(nr-k+1):(nr), (1):(nc)] <- sweep(x [(nr-k+1):(nr), (1):(nc)], 1, rev(filter), '*')
    x [(1):(nr), (1):(k)] <- sweep(x [(1):(nr), (1):(k)], 2, filter, '*')
    x [(1):(nr), (nc-k+1):(nc)] <- sweep(x [(1):(nr), (nc-k+1):(nc)], 2, rev(filter), '*')
    return(x)
}
