WOI <-
function(x = x, s = c(1,3), l = c(4,7), thres = 0.1, flat = 25, verbose = FALSE, periodic = FALSE){
    
    ## pre-processing array
    tstart <- Sys.time()
    if (max(s) > max(l)) {
        stop("max(s) is larger than max(l)!")
    }
    if (min(s) > min(l)) {
        stop("min(s) is larger than min(l)!")
    }
    if (max(s) >= min(l)) {
        stop("max(s) is larger than or equal min(l)!")
    }
    nx <- dim(x)[1]
    ny <- dim(x)[2]
    if (max(nx, ny) < 8) {
        stop("Sorry, but dim(x) is too small! dim(x) should be >= 8.")
    }
    if (max(nx, ny) > 4096) {
        stop("Sorry, but dim(x) is too large! dim(x) should be < 4096.")
    }
    quad = nx == ny && log2(nx) %% 1 == 0
    if (quad) {
        if (verbose && periodic) print(paste("Object size is ", dim(x)[1], "x", dim(x)[2], " and is periodic. -> WOI calculation needs no mirroring or flatten. Start: ", date(), sep = ""))
        flat = 0
        if (! periodic) {
        if (verbose) print(paste("Object size is ", dim(x)[1], "x", dim(x)[2], ". -> WOI calculation needs mirroring. Start: ", date(), sep = ""))
        x <- buildperiodic(x = x)
        }
    } else {
        if (min(nx, ny) <= 2 * flat) {
        stop(paste("Sorry, but dim(x) is too small! dim(x) should be >", 2 * flat, "."))
        }
        if (verbose) print(paste("Object size is ", dim(x)[1], "x", dim(x)[2], ". -> WOI calculation needs flatten and blowup. Start: ", date(), sep = ""))
        x <- flatten(x = x, filter = seq(0, 1, , flat))
        x <- blowup(x = x, M = 2 ^ ceiling(log2(max(nx, ny))))
    }
    r <- sum(x > thres)

    ## wavelet transform
    if (verbose) print(paste("Wavelet transform of ", dim(x)[1], "x",  dim(x)[2], " object.", sep = ""))
    wav = wavtra(x = x)
    if (quad) {
    wav <- wav [1:nx, 1:ny, ]
    x <- x [1:nx, 1:ny]
    } else {
    mask <- array(1, dim = c(nx, ny))
    mask <- blowup(x = mask, M = 2 ^ ceiling(log2(max(nx, ny))), number = NA)
    ind <- apply (which (!is.na (mask), arr.ind = T), 2, range, na.rm = T)
    wav <- wav [(ind [1,1] + flat):(ind [2,1] - flat), (ind [1,2] + flat):(ind [2,2] - flat), ]
    x <- x [(ind [1,1] + flat):(ind [2,1] - flat), (ind [1,2] + flat):(ind [2,2] - flat)]
    }

    ## calculate WOI via Fortran
    if (verbose) print(paste("Calculate WOI and LWOI.", sep = ""))
    RR <- x
    RR [RR < thres] <- -9999.
    r <- sum(RR > 0)
    WOI = .Fortran ("woifortran", nx = as.integer (dim (wav)[1]), ny = as.integer (dim (wav)[2]), nz = as.integer (dim (wav)[3]), r = as.integer (r), s1 = as.integer (s [1]), s2 = as.integer (s [2]), l1 = as.integer (l [1]), l2 = as.integer (l [2]), wav = wav, WOI = array (-9999, dim = c (as.integer (6))), PACKAGE = "calcWOI")$WOI
    LWOI = .Fortran ("lwoifortran", nx = as.integer (dim (wav)[1]), ny = as.integer (dim (wav)[2]), nz = as.integer (dim (wav)[3]), RR = RR, s1 = as.integer (s [1]), s2 = as.integer (s [2]), l1 = as.integer (l [1]), l2 = as.integer (l [2]), wav = wav, WOI = array (-9999., dim = c (as.integer (dim (wav)[1]), as.integer (dim (wav)[2]), as.integer (3))), PACKAGE = "calcWOI")$WOI

    ## post-preocessing WOI
    if (verbose) print(paste("Post-processing WOI and LWOI.", sep = ""))
    WOI [WOI == -9999] <- NA
    LWOI [LWOI == -9999] <- NA
    RR [RR == -9999] <- NA
    
    ## finished
    if (verbose) print(paste("Finished WOI and LWOI. End:", date()))
    l <- list(WOI1orig = WOI [1], WOI2orig = WOI [2], WOI3orig = WOI [3], WOIorig = log (WOI [1] * WOI [2] * WOI [3]), WOI1 = WOI [4], WOI2 = WOI [5], WOI3 = WOI [6], WOI = 1/3 * (WOI [4] + WOI [5] * WOI [6]), LWOI1 = LWOI [, , 1], LWOI2 = LWOI [, , 2], LWOI3 = LWOI [, , 3], LWOI = 1/3 * (LWOI [, , 1] + LWOI [, , 2] + LWOI [, , 3]), s = s, l = l, flat = flat, quad = quad, thres = thres, RR = RR, ts = round(difftime(Sys.time(), as.POSIXct(tstart), units = "secs"), 3))
    if (verbose) print(paste("WOI calculation took", l$ts, "seconds."))
    return(l)
}
