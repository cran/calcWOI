if(getRversion() >= "2.15.1")  utils::globalVariables(c("AICEN"))
wavtra <-
function (x) {
    # Most of this code is copied from the cddews function within the LS2W package
    filter.number <- 4
    family <- "DaubExPhase"
    data.wd <- imwd(x, filter.number = filter.number, family = family, type = "station")
    RawPer <- getdata(data.wd, switch = "direction")
    first.last.d <- data.wd$fl.dbase$first.last.d
    first.last.c <- data.wd$fl.dbase$first.last.c
    firstD <- first.last.d[data.wd$nlevels, 1]
    lastD <- first.last.d[data.wd$nlevels, 2]
    LengthD <- lastD - firstD + 1
    LEVELS <- data.wd$nlevels
    AI <- AICEN [[paste("AI", 2 ^ LEVELS, sep = "")]]
    centers <- AICEN [[paste("CEN", 2 ^ LEVELS, sep = "")]]
    TMP <- matrix(aperm(RawPer), nrow = 3 * LEVELS, ncol = LengthD ^ 2, byrow = TRUE)
    TMP2 <- AI %*% TMP
    data2 <- array(0, dim(RawPer))
    for (i in (1:(3 * LEVELS))) {
        data2[i, , ] <- matrix(TMP2[i, ], nrow = LengthD, ncol = LengthD, byrow = TRUE)
        data2[i, , ] <- shiftmat (x = data2[i, , ], dx = as.integer ((2 ^ LEVELS)/2 - centers$x [i]), dy = as.integer ((2 ^ LEVELS)/2 - centers$y [i]))
    }
    data2 <- aperm(data2, c(2, 3, 1))
    return(data2)
}
