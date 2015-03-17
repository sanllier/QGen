#if !defined( RANDOM_CURAND_H ) && defined( CURAND )
#define RANDOM_CURAND_H

#include "random.h"
#include <stdlib.h>

//------------------------------------------------------------

class RandomCURand: public IRandom
{
public:
    RandomCURand();
    ~RandomCURand();

    void setSize( long long size );
    long long size() const { return m_size; }
    void setSeed( unsigned seed ) OVERRIDE;

    BASETYPE next() OVERRIDE;

private:
    void generate() const;

private:
    void* m_gen;

    long long m_size;
    long long m_bufPos;
    BASETYPE* m_randomBufGPU;
    BASETYPE* m_randomBufCPU;
};

//------------------------------------------------------------
#endif
