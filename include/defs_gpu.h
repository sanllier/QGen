#ifndef DEFS_GPU_H
#define DEFS_GPU_H
//------------------------------------------------------------

#ifndef BASETYPE
    #define BASETYPE float
#endif

//------------------------------------------------------------

#define GET_THETA( _arr_, _x_, _y_, _z_ ) ( *(const BASETYPE *)( (_arr_) + (_x_) * 4 + (_y_) * 2 + (_z_) ) )

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
