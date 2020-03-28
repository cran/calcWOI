LW <- function(x, thres = 0.1, Nx = 2^ceiling(log2(max(dim(x)))), Ny = 2^ceiling(log2(max(dim(x)))), boundaries = "pad", verbose = FALSE) {

  # array checks
  if (any(is.na(x))) stop("x contains NA or NaN values!")
  if (Nx < 4) stop("Nx is too small! Use Nx >= 4.")
  if (Ny < 4) stop("Ny is too small! Use Ny >= 4.")
  if (Nx > 1024) stop("dim(x)[1] is too large! Use dimensions <= 1024.")
  if (Ny > 1024) stop("dim(x)[2] is too large! Use dimensions <= 1024.")
  
  # start calculation
  tstart <- Sys.time()
  if (verbose) print(paste(date(), "Set mask array, dimensions:", dim(x)[1], "x",  dim(x)[2]))
  # mask
  mask <- x
  mask [x < thres] <- NA
  mask <- mask / mask
  # wavelet transform
  if (verbose) print(paste(date(), "Start wavelet transform of size:", Nx, "x", Ny))
  E <- fld2dt(x, Nx = Nx, Ny = Ny, correct = TRUE, boundaries = boundaries)
  # central 
  if (verbose) print(paste(date(), "Calculate LWsc and LWai."))
  cen <- dt2cen(E)
  # direction components
  if (verbose) print(paste(date(), "Calculate LWuu and LWvv."))
  uv <- cen2uv(cen)
  # LWOI components
  if (verbose) print(paste(date(), "Calculate LWsc."))
  LWsc <- (cen [, , 3]-1) / (ceiling(log2(max(dim(x))))-3-1)
  LWin <- 1 - exp(-apply(E, 2:3, mean))
  LWai <- cen [, , 1]
  LWuu <- uv [, , 1]
  LWvv <- uv [, , 2]
  angle <- cen [, , 2]

  if (verbose) print(paste(date(), "LW calculation took", round(difftime(Sys.time(), as.POSIXct(tstart), units = "secs"), 3), "seconds."))
  l <- list(LWsc = LWsc*mask, LWin = LWin*mask, LWai = LWai*mask, LWuu = LWuu*mask, LWvv = LWvv*mask, angle = angle*mask, thres = thres, mask = mask, x = x, ts = round(difftime(Sys.time(), as.POSIXct(tstart), units = "secs"), 3))
  return(l)
}
