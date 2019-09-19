module hipmodule
  use, intrinsic :: iso_c_binding

  implicit none

  integer :: mystream
  type(c_ptr) :: stream, event
  !$omp threadprivate(stream, mystream, event)

  enum, bind(c)
     enumerator :: hipMemcpyHostToHost = 0, hipMemcpyHostToDevice = 1
     enumerator :: hipMemcpyDeviceToHost = 2, hipMemcpyDeviceToDevice =3
     enumerator :: hipMemcpyDefault = 4
  end enum

  enum, bind(c)
     enumerator :: hipSuccess = 0, hipErrorOutOfMemory = 2
     enumerator :: hipErrorNotInitialized = 3, hipErrorDeinitialized = 4
     enumerator :: hipErrorProfilerDisabled = 5, hipErrorProfilerNotInitialized = 6
     enumerator :: hipErrorProfilerAlreadyStarted = 7, hipErrorProfilerAlreadyStopped = 8
     enumerator :: hipErrorInvalidImage = 200
     enumerator :: hipErrorInvalidContext = 201, hipErrorContextAlreadyCurrent = 202
     enumerator :: hipErrorMapFailed = 205, hipErrorUnmapFailed = 206
     enumerator :: hipErrorArrayIsMapped = 207, hipErrorAlreadyMapped = 208
     enumerator :: hipErrorNoBinaryForGpu = 209, hipErrorAlreadyAcquired = 210
     enumerator :: hipErrorNotMapped = 211, hipErrorNotMappedAsArray = 212
     enumerator :: hipErrorNotMappedAsPointer = 213, hipErrorECCNotCorrectable = 214
     enumerator :: hipErrorUnsupportedLimit = 215, hipErrorContextAlreadyInUse = 216
     enumerator :: hipErrorPeerAccessUnsupported = 217, hipErrorInvalidKernelFile = 218
     enumerator :: hipErrorInvalidGraphicsContext = 219, hipErrorInvalidSource = 300
     enumerator :: hipErrorFileNotFound = 301, hipErrorSharedObjectSymbolNotFound = 302
     enumerator :: hipErrorSharedObjectInitFailed = 303, hipErrorOperatingSystem = 304
     enumerator :: hipErrorSetOnActiveProcess = 305, hipErrorInvalidHandle = 400
     enumerator :: hipErrorNotFound = 500, hipErrorIllegalAddress = 700
     enumerator :: hipErrorInvalidSymbol = 701, hipErrorMissingConfiguration = 1001
     enumerator :: hipErrorMemoryAllocation = 1002, hipErrorInitializationError = 1003
     enumerator :: hipErrorLaunchFailure = 1004, hipErrorPriorLaunchFailure = 1005
     enumerator :: hipErrorLaunchTimeOut = 1006, hipErrorLaunchOutOfResources = 1007
     enumerator :: hipErrorInvalidDeviceFunction = 1008, hipErrorInvalidConfiguration = 1009
     enumerator :: hipErrorInvalidDevice = 1010, hipErrorInvalidValue = 1011
     enumerator :: hipErrorInvalidDevicePointer = 1017, hipErrorInvalidMemcpyDirection = 1021
     enumerator :: hipErrorUnknown = 1030, hipErrorInvalidResourceHandle = 1033
     enumerator :: hipErrorNotReady = 1034, hipErrorNoDevice = 1038
     enumerator :: hipErrorPeerAccessAlreadyEnabled = 1050, hipErrorPeerAccessNotEnabled = 1051
     enumerator :: hipErrorRuntimeMemory = 1052, hipErrorRuntimeOther = 1053
     enumerator :: hipErrorHostMemoryAlreadyRegistered = 1061
     enumerator :: hipErrorHostMemoryNotRegistered = 1062, hipErrorMapBufferObjectFailed = 1071
     enumerator :: hipErrorTbd
  end enum

  interface

     function hipMalloc(ptr,sizeBytes) bind(c, name="hipMalloc")
       use iso_c_binding
       implicit none
       integer(c_int) :: hipMalloc
       type(c_ptr) :: ptr
       integer(c_size_t), value :: sizeBytes
     end function hipMalloc

     function hipFree(ptr) bind(c, name="hipFree")
       use iso_c_binding
       implicit none
       integer(c_int) :: hipFree
       type(c_ptr),value :: ptr
     end function hipFree

     function hipMemcpy(dst,src,sizeBytes,cpykind) bind(c,name="hipMemcpy")
       use iso_c_binding
       implicit none
       integer(c_int) :: hipMemcpy
       type(c_ptr),value :: dst
       type(c_ptr),value :: src
       integer(c_size_t), value :: sizeBytes

       ! We want to make sure we get the right integer for the enum
       integer(c_int), value :: cpykind

     end function hipMemcpy

     function hipDeviceSynchronize() bind(c, name="hipDeviceSynchronize")
       use iso_c_binding
       implicit none
       integer(c_int) :: hipDeviceSynchronize
     end function hipDeviceSynchronize

    integer(c_int) function &
        & hipStreamCreate(stream) &
        & bind(c, name="hipStreamCreate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: stream
    end function hipStreamCreate

    integer(c_int) function &
        & hipStreamDestroy(stream) &
        & bind(c, name="hipStreamDestroy")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
    end function hipStreamDestroy

    integer(c_int) function &
        & hipStreamSynchronize(stream) &
        & bind(c, name="hipStreamSynchronize")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: stream
    end function hipStreamSynchronize

    integer(c_int) function &
        & hipSetDevice(device) &
        & bind(c, name="hipSetDevice")
      use, intrinsic :: iso_c_binding
      integer(c_int), value :: device
    end function hipSetDevice

    integer(c_int) function &
        & hipGetDeviceCount(deviceCount) &
        & bind(c, name="hipGetDeviceCount")
      use, intrinsic :: iso_c_binding
      integer(c_int) :: deviceCount
    end function hipGetDeviceCount

  end interface
end module hipmodule
