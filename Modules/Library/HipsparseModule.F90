!***************************************************************************************************
! HipsparseModule.f90 10/18/17
! this file contains the module defining fortran interfaces for the hipblas library
!***************************************************************************************************

module HipsparseModule
  !-------------------------------------------------------------------------------------------------
  ! Interface to hipsparse routines
  !-------------------------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding

  implicit none

  type(c_ptr) :: hipsparse_handle
  !!$omp threadprivate(hipsparse_handle)

  enum, bind(c) !:: hipsparseIndexBase_t
    enumerator :: HIPSPARSE_INDEX_BASE_ZERO = 0
    enumerator :: HIPSPARSE_INDEX_BASE_ONE = 1
  end enum !hipsparseIndexBase_t

  enum, bind(c) !:: hipsparseDiagType_t
    enumerator :: HIPSPARSE_DIAG_TYPE_NON_UNIT = 0
    enumerator :: HIPSPARSE_DIAG_TYPE_UNIT = 1
  end enum !hipsparseDiagType_t

  enum, bind(c) !:: hipsparseFillMode_t
    enumerator :: HIPSPARSE_FILL_MODE_LOWER = 0
    enumerator :: HIPSPARSE_FILL_MODE_UPPER = 1
  end enum !hipsparseFillMode_t

  enum, bind(c) !:: hipsparseOperation_t
    enumerator :: HIPSPARSE_OPERATION_NON_TRANSPOSE = 0
    enumerator :: HIPSPARSE_OPERATION_TRANSPOSE = 1
    enumerator :: HIPSPARSE_OERATIONP_CONJUGATE_TRANSPOSE = 2
  end enum !hipsparseOperation_t

  enum, bind(c) !:: hipsparseStatus_t
    enumerator :: HIPSPARSE_STATUS_SUCCESS = 0
    enumerator :: HIPSPARSE_STATUS_NOT_INITIALIZED  = 1
    enumerator :: HIPSPARSE_STATUS_ALLOC_FAILED = 2
    enumerator :: HIPSPARSE_STATUS_INVALID_VALUE = 3
    enumerator :: HIPSPARSE_STATUS_ARCH_MISMATCH = 4
    enumerator :: HIPSPARSE_STATUS_MAPPING_ERROR = 5
    enumerator :: HIPSPARSE_STATUS_EXECUTION_FAILED = 6
    enumerator :: HIPSPARSE_STATUS_INTERNAL_ERROR = 7
    enumerator :: HIPSPARSE_MATRIX_TYPE_NOT_SUPPORTED = 8
    enumerator :: HIPSPARSE_STATUS_ZERO_PIVOT = 9
  end enum !hipsparseStatus_t

  interface

    function hipsparseCreate(handle) &
        & bind(c, name="hipsparseCreate")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: hipsparseCreate
      type(c_ptr) :: handle
    end function hipsparseCreate

    function hipsparseDestroy(handle) &
        & bind(c, name="hipsparseDestroy")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: hipsparseDestroy
      type(c_ptr), value :: handle
    end function hipsparseDestroy

    function hipsparseSetStream(handle, stream) &
        & bind(c, name="hipsparseSetStream")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: hipsparseSetStream
      type(c_ptr), value :: handle
      type(c_ptr), value :: stream
    end function hipsparseSetStream

    function hipsparseDgthr(handle, nnz, y, xval, xind, idxbase ) &
        & bind(c, name="hipsparseDgthr")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: hipsparseDgthr
      type(c_ptr), value :: handle
      integer(c_int), value :: nnz
      type(c_ptr), value :: y
      type(c_ptr), value :: xval
      type(c_ptr), value :: xind
      integer(c_int), value :: idxbase
    end function hipsparseDgthr

  end interface

end module HipsparseModule
