module hipblasmodule
  use iso_c_binding

  implicit none

  type(c_ptr) :: hipblas_handle

  enum, bind(c)
    enumerator :: HIPBLAS_STATUS_SUCCESS = 0 ! Function succeeds
    enumerator :: HIPBLAS_STATUS_NOT_INITIALIZED = 1 ! HIPBLAS library not initialized
    enumerator :: HIPBLAS_STATUS_ALLOC_FAILED = 2 ! resource allocation failed
    enumerator :: HIPBLAS_STATUS_INVALID_VALUE = 3 ! unsupported numerical value was passed to function
    enumerator :: HIPBLAS_STATUS_MAPPING_ERROR = 4 ! access to GPU memory space failed
    enumerator :: HIPBLAS_STATUS_EXECUTION_FAILED = 5 ! GPU program failed to execute
    enumerator :: HIPBLAS_STATUS_INTERNAL_ERROR = 6 ! an internal HIPBLAS operation failed
    enumerator :: HIPBLAS_STATUS_NOT_SUPPORTED = 7 ! function not implemented
    enumerator :: HIPBLAS_STATUS_ARCH_MISMATCH = 8 !
    enumerator :: HIPBLAS_STATUS_HANDLE_IS_NULLPTR = 9 !hipBLAS handle is null pointer
 end enum

 enum, bind(c) !:: hipblasFillMode_t
    enumerator :: HIPBLAS_FILL_MODE_LOWER = 121
    enumerator :: HIPBLAS_FILL_MODE_UPPER = 122
    enumerator :: HIPBLAS_FILL_MODE_FULL =  123
 end enum !hipblasFillMode_t

 enum, bind(c) !:: hipblasDiag type_t
    enumerator :: HIPBLAS_DIAG_NON_UNIT = 131
    enumerator :: HIPBLAS_DIAG_UNIT = 132
 end enum !hipblasDiag    type_t

 enum, bind(c) !:: hipblasSideMode_t
    enumerator :: HIPBLAS_SIDE_LEFT = 141
    enumerator :: HIPBLAS_SIDE_RIGHT = 142
    enumerator :: HIPBLAS_SIDE_BOTH = 143
 end enum !hipblasSideMode_t

 enum, bind(c)
    enumerator :: HIPBLAS_OP_N = 111
    enumerator :: HIPBLAS_OP_T = 112
    enumerator :: HIPBLAS_OP_C = 113
 end enum

  interface

    integer(c_int) function &
        & hipblasCreate(handle) &
        & bind(c, name="hipblasCreate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: handle
    end function hipblasCreate

    integer(c_int) function &
        & hipblasDestroy(handle) &
        & bind(c, name="hipblasDestroy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
    end function hipblasDestroy

    integer(c_int) function &
        & hipblasGetStream(handle, stream) &
        & bind(c, name="hipblasGetStream")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      type(c_ptr) :: stream
    end function hipblasGetStream

    integer(c_int) function &
        & hipblasSetStream(handle, stream) &
        & bind(c, name="hipblasSetStream")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      type(c_ptr), value :: stream
    end function hipblasSetStream

    integer(c_int) function &
        & hipblasGetVector(n, elemSize, dx_src, incx, hy_dst, incy) &
        & bind(c, name="hipblasGetVector")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: hy_dst
      integer(c_int), value :: incy
    end function hipblasGetVector

    integer(c_int) function &
        & hipblasGetVectorAsync(n, elemSize, dx_src, incx, hy_dst, incy, stream) &
        & bind(c, name="hipblasGetVectorAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: hy_dst
      integer(c_int), value :: incy
      type(c_ptr), value :: stream
    end function hipblasGetVectorAsync

    integer(c_int) function &
        & hipblasSetVector(n, elemSize, hx_src, incx, dy_dst, incy) &
        & bind(c, name="hipblasSetVector")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: dy_dst
      integer(c_int), value :: incy
    end function hipblasSetVector

    integer(c_int) function &
        & hipblasSetVectorAsync(n, elemSize, hx_src, incx, dy_dst, incy, stream) &
        & bind(c, name="hipblasSetVectorAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: n
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hx_src
      integer(c_int), value :: incx
      type(c_ptr), value :: dy_dst
      integer(c_int), value :: incy
      type(c_ptr), value :: stream
    end function hipblasSetVectorAsync

    integer(c_int) function &
        & hipblasSetMatrix(rows, cols, elemSize, hA_src, lda, dB_dst, lddb) &
        & bind(c, name="hipblasSetMatrix")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hA_src
      integer(c_int), value :: lda
      type(c_ptr), value :: dB_dst
      integer(c_int), value :: lddb
    end function hipblasSetMatrix

    integer(c_int) function &
        & hipblasSetMatrixAsync(rows, cols, elemSize, hA_src, lda, dB_dst, lddb, stream) &
        & bind(c, name="hipblasSetMatrixAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: hA_src
      integer(c_int), value :: lda
      type(c_ptr), value :: dB_dst
      integer(c_int), value :: lddb
      type(c_ptr), value :: stream
    end function hipblasSetMatrixAsync

    integer(c_int) function &
        & hipblasGetMatrix(rows, cols, elemSize, dA_src, ldda, hB_dst, ldb) &
        & bind(c, name="hipblasGetMatrix")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dA_src
      integer(c_int), value :: ldda
      type(c_ptr), value :: hB_dst
      integer(c_int), value :: ldb
    end function hipblasGetMatrix

    integer(c_int) function &
        & hipblasGetMatrixAsync(rows, cols, elemSize, dA_src, ldda, hB_dst, ldb, stream) &
        & bind(c, name="hipblasGetMatrixAsync")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: rows
      integer(c_int), value :: cols
      integer(c_size_t), value :: elemSize
      type(c_ptr), value :: dA_src
      integer(c_int), value :: ldda
      type(c_ptr), value :: hB_dst
      integer(c_int), value :: ldb
      type(c_ptr), value :: stream
    end function hipblasGetMatrixAsync

    integer(c_int) function &
        & hipblasDnrm2(handle, n, dx, incx, xnorm) &
        & bind(c, name="hipblasDnrm2")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      real(c_double) :: xnorm
    end function hipblasDnrm2

    integer(c_int) function &
        & hipblasDaxpy(handle, n, alpha, dx, incx, dy, incy) &
        & bind(c, name="hipblasDaxpy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      type(c_ptr), value :: dy
      integer(c_int), value :: incy
    end function hipblasDaxpy

    integer(c_int) function &
        & hipblasDgemv(handle, trans, m, n, alpha, dA, ldda, dx, incx, beta, dy, incy) &
        & bind(c, name="hipblasDgemv")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: trans
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      real(c_double) :: beta
      type(c_ptr), value :: dy
      integer(c_int), value :: incy
    end function hipblasDgemv

    integer(c_int) function &
        & hipblasDgetrfBatched(handle, n, dA, ldda, dP, dInfo, nbatch) &
        & bind(c, name="hipblasDgetrfBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dP
      type(c_ptr), value :: dInfo
      integer(c_int), value :: nbatch
    end function hipblasDgetrfBatched

    integer(c_int) function &
        & hipblasDgemmBatched(handle, transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc, nbatch) &
        & bind(c, name="hipblasDgemmBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
      integer(c_int), value :: nbatch
    end function hipblasDgemmBatched

    integer(c_int) function &
        & hipblasDgemmStridedBatched(handle, transa, transb, m, n, k, alpha, &
        & dA, ldda, strideA, dB, lddb, strideB, beta, dC, lddc, strideC, nbatch) &
        & bind(c, name="hipblasDgemmStridedBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      integer(c_int), value :: strideA
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      integer(c_int), value :: strideB
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: strideC
      integer(c_int), value :: lddc
      integer(c_int), value :: nbatch
    end function hipblasDgemmStridedBatched

    integer(c_int) function &
        & hipblasDtrsv(uplo, trans, diag, n, dA, ldda, dx, incx) &
        & bind(c, name="hipblasDtrsv")
      use, intrinsic :: iso_c_binding
      character(c_char), value :: uplo
      character(c_char), value :: trans
      character(c_char), value :: diag
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
    end function hipblasDtrsv

    integer(c_int) function &
        & hipblasDtrsv(handle, uplo, trans, diag, n, dA, ldda, dx, incx) &
        & bind(c, name="hipblasDtrsv")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: uplo
      integer(c_int), value :: trans
      integer(c_int), value :: diag
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
    end function hipblasDtrsv

    integer(c_int) function &
        & hipblasDtrsm(side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb) &
        & bind(c, name="hipblasDtrsm")
      use, intrinsic :: iso_c_binding
      character(c_char), value :: side
      character(c_char), value :: uplo
      character(c_char), value :: trans
      character(c_char), value :: diag
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
    end function hipblasDtrsm

    integer(c_int) function &
        & hipblasDtrsm(handle, side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb) &
        & bind(c, name="hipblasDtrsm")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: side
      integer(c_int), value :: uplo
      integer(c_int), value :: trans
      integer(c_int), value :: diag
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
    end function hipblasDtrsm

    integer(c_int) function &
        & hipblasDgemm(transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc) &
        & bind(c, name="hipblasDgemm")
      use, intrinsic :: iso_c_binding
      character(c_char), value :: transa
      character(c_char), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double), value :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double), value :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function hipblasDgemm

    integer(c_int) function &
        & hipblasDgemm(handle, transa, transb, m, n, k, alpha, dA, ldda, dB, lddb, beta, dC, lddc) &
        & bind(c, name="hipblasDgemm")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      real(c_double) :: beta
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function hipblasDgemm

    integer(c_int) function &
        & hipblasDgetrsBatched(handle, trans, n, nrhs, dA, ldda, dP, dB, lddb, hInfo, nbatch) &
        & bind(c, name="hipblasDgetrsBatched")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: trans
      integer(c_int), value :: n
      integer(c_int), value :: nrhs
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dP
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      type(c_ptr), value :: hInfo
      integer(c_int), value :: nbatch
    end function hipblasDgetrsBatched

    integer(c_int) function &
        & hipblasDgeam(handle, transa, transb, m, n, alpha, dA, ldda, beta, dB, lddb, dC, lddc) &
        & bind(c, name="hipblasDgeam")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      real(c_double) :: alpha
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      real(c_double) :: beta
      type(c_ptr), value :: dB
      integer(c_int), value :: lddb
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function hipblasDgeam

    integer(c_int) function &
        & hipblasDdgmm(handle, mode, m, n, dA, ldda, dx, incx, dC, lddc) &
        & bind(c, name="hipblasDdgmm")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int), value :: mode
      integer(c_int), value :: m
      integer(c_int), value :: n
      type(c_ptr), value :: dA
      integer(c_int), value :: ldda
      type(c_ptr), value :: dx
      integer(c_int), value :: incx
      type(c_ptr), value :: dC
      integer(c_int), value :: lddc
    end function hipblasDdgmm

  end interface
end module hipblasmodule
