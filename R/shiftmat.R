shiftmat <-
function(x, dx = 0, dy = 0) {
    m = dim (x)[1]
    n = dim (x)[2]
    if (dy > 0) {
        x <- cbind(x [, (n-dy+1):n], x [, 1:(n-dy)])
    }
    if (dx > 0) {
        x <- rbind(x [(m-dx+1):m, ], x [1:(m-dx) ,]) 
    }
    if (dy < 0) {
        x <- cbind(x [, (abs(dy)+1):n], x [, 1:abs(dy)]) 
    }
    if (dx < 0) {
        x <- rbind (x [(abs(dx)+1):m, ], x [1:abs(dx), ]) 
    }
    return(x)
}
