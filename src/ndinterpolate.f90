!Copyright 2024 Paolo Lampitella
!This code is licensed under the terms of the MIT license

FUNCTION ndinterpolate(ndim,x,ng,grid,values)
!Multilinear interpolation in ndim dimensions
USE sps_vars
IMPLICIT NONE
INTEGER, INTENT(IN) :: ndim     !Number of dimensions
INTEGER, INTENT(IN) :: ng(ndim) !Number of points along each dimension (i.e., [nx, ny, ...])
REAL(SP),     INTENT(IN) :: x(ndim)  !The interpolation point
REAL(SP),     INTENT(IN) :: grid(SUM(ng)) !Coordinates, one dimension after the other (i.e., [x(:), y(:), ...])
REAL(SP),     INTENT(IN) :: values(PRODUCT(ng)) !Tabulated values, linearized in column major format
REAL(SP)     :: ndinterpolate, eta(2,ndim), wei, xb
INTEGER :: ind(ndim), i, j0, j1, j2, j
INTEGER :: kp, ki, ks, delta, middle

!Loop to find the cell (hypercube) of the table containing the interplation point x
!and compute the resulting interpolation factors eta
j0 = 0
DO i = 1, ndim
   j1     = j0 + 1
   j2     = j0 + ng(i)
   xb     = MIN(MAX(x(i),grid(j1)),grid(j2)) !Clipping the interpolation point to the grid
   !Binary search to find the cell containing xb
   DO
      delta  = j2-j1
      IF (delta<=1) EXIT
      middle = j1+delta/2
      IF (xb>grid(middle)) THEN
         j1  = middle
      ELSE
         j2 = middle
      ENDIF
   ENDDO
   ind(i)   = j1 - j0   !Index of the cell containing xb
   eta(2,i) = (xb - grid(j1))/(grid(j1+1)-grid(j1)) !Using x(i) here (instead of xb) performs extrapolation
   eta(1,i) = 1.0_SP-eta(2,i)
   j0       = j0 + ng(i)
ENDDO

ndinterpolate = 0.0_SP
!Loop over all the 2^ndim vertices of the cell containing the point
DO j  = 1, 2**ndim
   wei = 1.0_SP
   kp  = 1
   ki  = 1
   !Loop over the dimensions to retrieve the index and weight of the given node
   !We use the j bit pattern (which is made of ndim bits) to pick the lower/upper node in
   !each dimension of the given cell
   DO i = 1, ndim
      !This magic number will give 1 if the (i-1)-th bit of j-1 is set, 0  otherwise
      !It is derived from a more general formula for computing permutations with repetitions
      ks = MOD((2*j-1)/2**i,2)

      !Retrieve and use the weight for the i-th dimension of the j-th vertex
      wei = wei*eta(1+ks,i)

      !This is just the i-th step of sub2ind applied to sub(i)=ind(i)+ks which, in the end, returns
      !the linear index ki corresponding to subscripts sub(i) for a ndim-dimensional array
      !stored in column-major order (as it is values). It is defined as:
      !ki - 1 = SUM_i((sub(i)-1)*PROD_j(ng(j),j=1...i-1),i=1...ndim)
      ki = ki + (ind(i)+ks-1)*kp
      kp = kp*ng(i)
   ENDDO
   !Summing up this vertex contribution to the final interpolation
   ndinterpolate = ndinterpolate + wei*values(ki)
ENDDO
ENDFUNCTION ndinterpolate