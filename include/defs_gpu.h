#ifndef DEFS_GPU_H
#define DEFS_GPU_H
//------------------------------------------------------------

#define BASETYPE float

//------------------------------------------------------------

struct GPUQbit
{
    BASETYPE aReal;
    BASETYPE aImag;
    BASETYPE bReal;
    BASETYPE bImag;
};

//------------------------------------------------------------

static const int CUDA_BLOCK_SIZE = 64;

//------------------------------------------------------------
#endif
