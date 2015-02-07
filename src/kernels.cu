#include "cuda_runtime.h"
#include "cuda_error_handler.h"
#include "defs_gpu.h"

//------------------------------------------------------------

__global__ void setInitialKernel( GPUQbit* data, const BASETYPE* rands, long long size )
{
    long long index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index >= size )
        return;

    GPUQbit* targetElement = data + index;

    BASETYPE randVal = rands[ index ];
    targetElement->aReal = randVal;
    targetElement->aImag = BASETYPE(0);
    targetElement->bReal = sqrt( BASETYPE(1) - randVal );
    targetElement->bImag = BASETYPE(0);
}

extern "C" void launchSetInitialKernel( GPUQbit* data, const BASETYPE* rands, long long size )
{
    dim3 block = dim3 ( CUDA_BLOCK_SIZE );
    dim3 grid = dim3 ( size / CUDA_BLOCK_SIZE + 1 );
    SAFE_KERNEL_CALL( ( setInitialKernel<<< grid, block >>>( data, rands, size ) ) );
}

//------------------------------------------------------------

__global__ void observeKernel( const GPUQbit* data, const BASETYPE* rands, bool* buf, long long size )
{
    long long index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index >= size )
        return;

    const GPUQbit* targetElement = data + index;
    const BASETYPE real = targetElement->aReal;
    const BASETYPE imag = targetElement->aImag;
    const BASETYPE mod = abs( real * real + imag * imag );
    buf[ index ] = rands[ index ] >= mod;
}

extern "C" void launchObserveKernel( const GPUQbit* data, const BASETYPE* rands, bool* buf, long long size )
{
    dim3 block = dim3 ( CUDA_BLOCK_SIZE );
    dim3 grid = dim3 ( size / CUDA_BLOCK_SIZE + 1 );
    SAFE_KERNEL_CALL( ( observeKernel<<< grid, block >>>( data, rands, buf, size ) ) );
}

//------------------------------------------------------------

__global__ void thetaBufferKernel( const GPUQbit* data, const GPUQbit* best, const bool* obs, 
                                  const bool* bestObs, const BASETYPE* theta, const bool* isBetter,
                                  BASETYPE* out, long long size )
{
    long long index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index >= size )
        return;

    const bool curIndBit  = obs[ index ];
    const bool bestIndBit = bestObs[ index ];
    const bool betterThanBest = isBetter[ index ];

    BASETYPE frac = GET_THETA( theta, curIndBit ? 1:0, bestIndBit ? 1:0, betterThanBest ? 1:0 );
    const GPUQbit* targetElement = data + index;    
    BASETYPE prodRealPart = targetElement->aReal * targetElement->bReal - 
                            targetElement->aImag * targetElement->bImag;

    if ( prodRealPart >= 0 )
        out[ index ] = curIndBit ? frac : -frac;
    else if ( prodRealPart < 0 )
        out[ index ] = curIndBit ? -frac : +frac;
}

extern "C" void launchThetaBufferKernel( const GPUQbit* data, const GPUQbit* best, const bool* obs,
                                        const bool* bestObs, const BASETYPE* theta, const bool* isBetter,
                                        BASETYPE* out, long long size )
{
    dim3 block = dim3 ( CUDA_BLOCK_SIZE );
    dim3 grid = dim3 ( size / CUDA_BLOCK_SIZE + 1 );
    SAFE_KERNEL_CALL( ( thetaBufferKernel<<< grid, block >>>( data, best, obs, bestObs, theta, isBetter, out, size ) ) );
}

//------------------------------------------------------------

__global__ void evolveKernel( GPUQbit* data, const BASETYPE* theta, long long size )
{
    long long index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index >= size )
        return;

    const BASETYPE angle = theta[ index ];
    BASETYPE matrix[4];
    matrix[0] = cosf( angle );
    matrix[2] = sinf( angle );
    matrix[1] = -matrix[2];
    matrix[3] =  matrix[0];

    GPUQbit* targetElement = data + index; 

    BASETYPE resAReal = matrix[0] * targetElement->aReal + matrix[1] * targetElement->bReal;
    BASETYPE resAImag = matrix[0] * targetElement->aImag + matrix[1] * targetElement->bImag;
    BASETYPE resBReal = matrix[2] * targetElement->aReal + matrix[2] * targetElement->bReal;
    BASETYPE resBImag = matrix[3] * targetElement->aImag + matrix[3] * targetElement->bImag;
   
    targetElement->aReal = resAReal;
    targetElement->aImag = resAImag;
    targetElement->bReal = resBReal;
    targetElement->bImag = resBImag;          
}

extern "C" void launchEvolveKernel( GPUQbit* data, const BASETYPE* theta, long long size )
{
    dim3 block = dim3 ( CUDA_BLOCK_SIZE );
    dim3 grid = dim3 ( size / CUDA_BLOCK_SIZE + 1 );
    SAFE_KERNEL_CALL( ( evolveKernel<<< grid, block >>>( data, theta, size ) ) );
}

//------------------------------------------------------------
