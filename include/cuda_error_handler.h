#ifndef CUDA_ERROR_HANDLER_H
#define CUDA_ERROR_HANDLER_H

#include <string>

//---------------------------------------------------------------

#define SAFE_CALL( CallInstruction ) \
{ \
    cudaError_t cuerr = CallInstruction; \
    if(cuerr != cudaSuccess) \
        throw std::string( "CUDA error: " ).append( cudaGetErrorString(cuerr) ).append( " at call " ).append( #CallInstruction ); \
}

//---------------------------------------------------------------

#define SAFE_CURAND_CALL( CallInstruction ) \
{ \
    curandStatus_t cuerr = CallInstruction; \
    if(cuerr != CURAND_STATUS_SUCCESS) \
        throw std::string( "CURAND error at call " ).append( #CallInstruction ); \
}

//---------------------------------------------------------------

#define SAFE_KERNEL_CALL( KernelCallInstruction ) \
{ \
    KernelCallInstruction; \
    cudaError_t cuerr = cudaGetLastError(); \
    if(cuerr != cudaSuccess) \
        throw std::string( "CUDA error in kernel launch: " ).append( cudaGetErrorString(cuerr) ).append( " at kernel " ).append( #KernelCallInstruction ); \
    \
    cuerr = cudaDeviceSynchronize(); \
    if(cuerr != cudaSuccess) \
        throw std::string( "CUDA error in kernel execution: " ).append( cudaGetErrorString(cuerr) ).append( " at kernel " ).append( #KernelCallInstruction ); \
}

//---------------------------------------------------------------

#endif
