#include "cuda_runtime.h"
#include "cuda_error_handler.h"
#include "defs_gpu.h"

//------------------------------------------------------------

__global__ void setInitialKernel( GPUQbit* data, long long size )
{
    long long index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index >= size )
        return;

    printf( "%d\r\n", size );

    GPUQbit* targetElement = data + index;

    BASETYPE randVal = 0.0; // CRAP
    targetElement->aReal = randVal;
    targetElement->aImag = BASETYPE(0);
    targetElement->bReal = sqrt( BASETYPE(1) - randVal );
    targetElement->bImag = BASETYPE(0);
}

void launchSetInitialKernel( GPUQbit* data, long long size )
{
    dim3 block = dim3 ( CUDA_BLOCK_SIZE );
    dim3 grid = dim3 ( size / CUDA_BLOCK_SIZE + 1 );
    SAFE_KERNEL_CALL( ( setInitialKernel<<< grid, block >>>( ( GPUQbit* )data, size ) ) );
}

//------------------------------------------------------------
