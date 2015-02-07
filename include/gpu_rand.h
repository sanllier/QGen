#if !defined( GPU_RAND_H ) && defined( GPU )
#define GPU_RAND_H

#include "defs_gpu.h"
#include "curand.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class GPURand
{
public:
    GPURand( unsigned long long seed, long long size = 0 );
    ~GPURand();

    void generate();
    void resize( long long size );
    long long size() const { return m_size; }

    const BASETYPE* getGPUData() const { return m_randomBuf; }
    void copyInCPUData( BASETYPE* data ) const;

private:
    curandGenerator_t m_gen;

    long long m_size;
    BASETYPE* m_randomBuf;    
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
