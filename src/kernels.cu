#include "cuda_runtime.h"
#include "cuda_error_handler.h"
#include "defs.h"

//------------------------------------------------------------

__global__ void setInitialKernel( GPUQbit* data, long long size )
{
    long long index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index >= size )
        return;

    GPUQbit* targetElement = data + index;

    BASETYPE randVal = 0.0; // CRAP
    targetElement->aReal = randVal;
    targetElement->aImag = BASETYPE(0);
    targetElement->bReal = sqrt( BASETYPE(1) - randVal );
    targetElement->bimag = BASETYPE(0);    
}

void launchSetInitialKernel( QBit* data, long long size )
{
    dim3 block = dim3 ( 1, CUDA_BLOCK_SIZE );
    dim3 grid = dim3 ( 1, size / CUDA_BLOCK_SIZE + 1 );
    SAFE_KERNEL_CALL( ( setInitialKernel<<< grid, block >>>( ( GPUQbit* )data, size ) ) );
}

//------------------------------------------------------------
