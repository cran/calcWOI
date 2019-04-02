! compile with
! gfortran -fpic -g -O2 -c fortran.f95 -o fortran.o && gfortran -shared -L/usr/lib64/R/lib -L/usr/local/lib64 -o fortran.so fortran.o -L/usr/lib64/R/lib -lR

SUBROUTINE woifortran (nx, ny, nz, r, s1, s2, l1, l2, wav, WOI)
!
! calculates WOI and its components for given wavelet coefficients
!
implicit none
!
! declare some variables
!
integer (kind = 4), intent (in) :: nx, &  ! x-dimension of the wavelet field
  ny, &  ! y-dimension of the wavelet field
  nz, &  ! z-dimension of the wavelet field
  r, &  ! number of grid points with rain 
  s1, &  ! smallest small scale
  s2, &  ! largest small scale
  l1, &  ! smallest large scale
  l2  ! largest large scale
  
real (kind = 8), intent (in) :: wav (nx, ny, nz)  ! wavelet coefficients
  
real (kind = 8), intent (out) :: WOI (6)  ! field of WOI output
  
real (kind = 8) :: Es, &  ! mean small scale energy
  El, &  ! mean large scale energy
  Es1, &  ! small scale energy in East West direction
  Es2, &  ! small scale energy in North South direction
  Es3, &  ! small scale energy in diagonal direction
  El1, &  ! large scale energy in East West direction
  El2, &  ! large scale energy in North South direction
  El3, &  ! large scale energy in diagonal direction
  Ej (nz/3),  &  ! direction averaged energy
  E, &  ! total energy 
  W (nz), &  ! averaged energy 
  wavgt0 (nx,ny,nz)  ! wav gt 0
  
integer (kind = 4) :: Nj, i, j, s, l, z, k  ! loop

! number of scales Nj
Nj = nz/3

! initialize arrays
Es1 = 0.
Es2 = 0.
Es3 = 0.
El1 = 0.
El2 = 0.
El3 = 0.
W = 0.

! calculate domain average
DO z = 1, nz
  W (z) = sum (wav (:,:,z)) / (nx*ny)
END DO

! small scales
DO s = s1, s2
  Es1 = Es1 + W (s)
  Es2 = Es2 + W (s+Nj)
  Es3 = Es3 + W (s+2*Nj)
END DO

! large scales
DO l = l1, l2      
  El1 = El1 + W (l)
  El2 = El2 + W (l+Nj)
  El3 = El3 + W (l+2*Nj)
END DO
      
! calculate average
Es1 = Es1 / (s2-s1+1.)
Es2 = Es2 / (s2-s1+1.)
Es3 = Es3 / (s2-s1+1.)
El1 = El1 / (l2-l1+1.)
El2 = El2 / (l2-l1+1.)
El3 = El3 / (l2-l1+1.)

! calculate means over directions
Es = (Es1 + Es2 + Es3) / 3.
El = (El1 + El2 + El3) / 3.

! total energy averaged over directions
E = sum (W) / 3.

! direction averaged energy for each scale
DO z = 1, nz/3
  Ej (z) = (W (z) + W (z+Nj) + W (z+2*Nj)) / 3.
END DO

! WOI1orig 
WOI (1) = El / (Es + El)

! WOI2orig
WOI (2) = (Es + El) * nx*ny / r

! WOI3orig
WOI (3) = 1./3. * sqrt (((Es1 - Es)/Es)**2 &
        + ((El1 - El)/El)**2 &
        + ((Es2 - Es)/Es)**2 &
        + ((El2 - El)/El)**2 &
        + ((Es3 - Es)/Es)**2 &
        + ((El3 - El)/El)**2)


        
! same calculations without negative energy
wavgt0 = wav
DO i = 1, nx
  DO j = 1, ny
    DO k = 1, nz
      if (wavgt0 (i,j,k) .lt. 0.) wavgt0 (i,j,k) = 0.
    END DO
  END DO
END DO

! initialize arrays
Es1 = 0.
Es2 = 0.
Es3 = 0.
El1 = 0.
El2 = 0.
El3 = 0.
W = 0.

! calculate domain average
DO z = 1, nz
  W (z) = sum (wavgt0 (:,:,z)) / (nx*ny)
END DO

! small scales
DO s = s1, s2
  Es1 = Es1 + W (s)
  Es2 = Es2 + W (s+Nj)
  Es3 = Es3 + W (s+2*Nj)
END DO

! large scales
DO l = l1, l2      
  El1 = El1 + W (l)
  El2 = El2 + W (l+Nj)
  El3 = El3 + W (l+2*Nj)
END DO
      
! calculate average
Es1 = Es1 / (s2-s1+1.)
Es2 = Es2 / (s2-s1+1.)
Es3 = Es3 / (s2-s1+1.)
El1 = El1 / (l2-l1+1.)
El2 = El2 / (l2-l1+1.)
El3 = El3 / (l2-l1+1.)

! calculate means over directions
Es = (Es1 + Es2 + Es3) / 3.
El = (El1 + El2 + El3) / 3.

! total energy averaged over directions
E = sum (W) / 3.

! direction averaged energy for each scale
DO z = 1, nz/3
  Ej (z) = (W (z) + W (z+Nj) + W (z+2*Nj)) / 3.
END DO


        
! WOI1 
WOI (4) = El / (Es + El)

! WOI2
WOI (5) = 1. - exp (- ((sum (wav) / (nx*ny*nz)) / (1.*r/(nx*ny))))

! WOI3
WOI (6) = 1./ (2*sqrt(3.)) * sqrt (((Es1 - Es)/Es)**2 &
        + ((El1 - El)/El)**2 &
        + ((Es2 - Es)/Es)**2 &
        + ((El2 - El)/El)**2 &
        + ((Es3 - Es)/Es)**2 &
        + ((El3 - El)/El)**2)
        
if (WOI (4) .lt. 0. .or. WOI (4) .gt. 1.) WOI (4) = -9999.
if (WOI (5) .lt. 0. .or. WOI (5) .gt. 1.) WOI (5) = -9999.
if (WOI (6) .lt. 0. .or. WOI (6) .gt. 1.) WOI (6) = -9999.

END SUBROUTINE woifortran





SUBROUTINE lwoifortran (nx, ny, nz, RR, s1, s2, l1, l2, wav, WOI)
!
! calculates WOI and its components for given wavelet coefficients
!
implicit none
!
! declare some variables
!
integer (kind = 4), intent (in) :: nx, &  ! x-dimension of the wavelet field
  ny, &  ! y-dimension of the wavelet field
  nz, &  ! z-dimension of the wavelet field
  s1, &  ! smallest small scale
  s2, &  ! largest small scale
  l1, &  ! smallest large scale
  l2  ! largest large scale
  
real (kind = 8), intent (in) :: wav (nx, ny, nz), &  ! wavelet coefficients
  RR (nx, ny)  ! array of mask rain
  
real (kind = 8), intent (out) :: WOI (nx, ny, 3)  ! field of WOI output
  
real (kind = 8) :: Es (nx,ny), &  ! mean small scale energy
  El (nx,ny), &  ! mean large scale energy
  Es1 (nx,ny), &  ! small scale energy in East West direction
  Es2 (nx,ny), &  ! small scale energy in North South direction
  Es3 (nx,ny), &  ! small scale energy in diagonal direction
  El1 (nx,ny), &  ! large scale energy in East West direction
  El2 (nx,ny), &  ! large scale energy in North South direction
  El3 (nx,ny), &  ! large scale energy in diagonal direction
  Ej (nx,ny,nz/3),  &  ! direction averaged energy
  E (nx,ny), &  ! total energy
  wavgt0 (nx,ny,nz)  ! wav gt 0
  
integer (kind = 4) :: Nj, i, j, s, l, z, k  ! loop

integer, dimension (8) :: time

! number of scales Nj
Nj = nz/3

! remove negative energy
wavgt0 = wav*0.
DO i = 1, nx
  DO j = 1, ny
    IF (RR (i,j) .gt. -1.) THEN
      DO k = 1, nz
        if (wav (i,j,k) .gt. 0.) wavgt0 (i,j,k) = wav (i,j,k)
      END DO
    END IF
  END DO
END DO

! initialize arrays
Es1 (:,:) = 0.
Es2 (:,:) = 0.
Es3 (:,:) = 0.
El1 (:,:) = 0.
El2 (:,:) = 0.
El3 (:,:) = 0.

! small scales
DO i = 1, nx
  DO j = 1, ny
    IF (RR (i,j) .gt. -1.) THEN
      DO s = s1, s2
        Es1 (i,j) = Es1 (i,j) + wavgt0 (i,j,s)
        Es2 (i,j) = Es2 (i,j) + wavgt0 (i,j,s+Nj)
        Es3 (i,j) = Es3 (i,j) + wavgt0 (i,j,s+2*Nj)
      END DO
! large scales
      DO l = l1, l2
        El1 (i,j) = El1 (i,j) + wavgt0 (i,j,l)
        El2 (i,j) = El2 (i,j) + wavgt0 (i,j,l+Nj)
        El3 (i,j) = El3 (i,j) + wavgt0 (i,j,l+2*Nj)
      END DO
    END IF
  END DO
END DO

! calculate average
Es1 (:,:) = Es1 (:,:) / (s2-s1+1.)
Es2 (:,:) = Es2 (:,:) / (s2-s1+1.)
Es3 (:,:) = Es3 (:,:) / (s2-s1+1.)
El1 (:,:) = El1 (:,:) / (l2-l1+1.)
El2 (:,:) = El2 (:,:) / (l2-l1+1.)
El3 (:,:) = El3 (:,:) / (l2-l1+1.)

! calculate means over directions
Es (:,:) = (Es1 (:,:) + Es2 (:,:) + Es3 (:,:)) / 3.
El (:,:) = (El1 (:,:) + El2 (:,:) + El3 (:,:)) / 3.

! total energy averaged over directions
DO i = 1, nx
  DO j = 1, ny
    IF (RR (i,j) .gt. -1.) E (i,j) = sum (wavgt0 (i,j,:)) / 3.
  END DO
END DO

! direction averaged energy for each scale
DO z = 1, nz/3
  DO i = 1, nx
    DO j = 1, ny
      IF (RR (i,j) .gt. -1.) THEN
        Ej (i,j,z) = (wavgt0 (i,j,z) + wavgt0 (i,j,z+Nj) + wavgt0 (i,j,z+2*Nj)) / 3.
      END IF
    END DO
  END DO
END DO

! LWOI1
DO i = 1, nx
  DO j = 1, ny
    IF (RR (i,j) .gt. -1.) WOI (i,j,1) = El (i,j) / (Es (i,j) + El (i,j))
  END DO
END DO

! LWOI2
DO i = 1, nx
  DO j = 1, ny
    IF (RR (i,j) .gt. -1.) WOI (i,j,2) = 1. - exp (- (sum (wavgt0 (i,j,:)) / (nz)))
  END DO
END DO

! LWOI3
DO i = 1, nx
  DO j = 1, ny
    IF (RR (i,j) .gt. -1.) THEN
      WOI (i,j,3) = 1. / (2.*sqrt (3.)) * sqrt (((Es1 (i,j) - Es (i,j))/Es (i,j))**2 &
              + ((El1 (i,j) - El (i,j))/El (i,j))**2 &
              + ((Es2 (i,j) - Es (i,j))/Es (i,j))**2 &
              + ((El2 (i,j) - El (i,j))/El (i,j))**2 &
              + ((Es3 (i,j) - Es (i,j))/Es (i,j))**2 &
              + ((El3 (i,j) - El (i,j))/El (i,j))**2)
    END IF
  END DO
END DO



! missing values
DO i = 1, nx
  DO j = 1, ny
    if (WOI (i,j,1) .lt. 0. .or. WOI (i,j,1) .gt. 1.) WOI (i,j,1) = -9999.
    if (WOI (i,j,2) .lt. 0. .or. WOI (i,j,2) .gt. 1.) WOI (i,j,2) = -9999.
    if (WOI (i,j,3) .lt. 0. .or. WOI (i,j,3) .gt. 1.) WOI (i,j,3) = -9999.
  END DO
END DO

END SUBROUTINE lwoifortran



