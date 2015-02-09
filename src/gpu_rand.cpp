#ifdef GPU

#include "gpu_rand.h"
#include "curand.h"
#include "cuda_error_handler.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

GPURand::GPURand( unsigned long long seed, long long size/* = 0*/ )
    : m_gen(0)
    , m_randomBuf(0)
    , m_size(0)
{
    m_gen = new curandGenerator_t;

    SAFE_CURAND_CALL( curandCreateGenerator( ( curandGenerator_t * )m_gen, CURAND_RNG_PSEUDO_DEFAULT ) );
    SAFE_CURAND_CALL( curandSetPseudoRandomGeneratorSeed( *( curandGenerator_t * )m_gen, seed ) );
    if ( size > 0 )
        resize( size );
}

//------------------------------------------------------------

GPURand::~GPURand()
{
    SAFE_CURAND_CALL( curandDestroyGenerator( *( curandGenerator_t * )m_gen ) );
    SAFE_CALL( cudaFree( m_randomBuf ) );
    delete ( curandGenerator_t * )m_gen;
}

//------------------------------------------------------------

void GPURand::generate() const
{
   if ( m_size <= 0 )
       throw std::string( "GPURand trying to generate with unsetted size. " ).append( __FUNCTION__ );

    SAFE_CURAND_CALL( curandGenerateUniform( *( curandGenerator_t * )m_gen, m_randomBuf, size_t( m_size ) ));
}

//------------------------------------------------------------

void GPURand::resize( long long size )
{
   if ( size <= 0 )
       throw std::string( "GPURand trying to resize with incorrect parameter. " ).append( __FUNCTION__ );

    if ( m_randomBuf )
        SAFE_CALL( cudaFree( m_randomBuf ) );

    m_size = size;
    SAFE_CALL( cudaMalloc( &m_randomBuf, size_t( m_size * sizeof( BASETYPE ) ) ) );
}

//------------------------------------------------------------

void GPURand::copyInCPUData( BASETYPE* data ) const
{
    SAFE_CALL( cudaMemcpy( data, m_randomBuf, size_t( m_size * sizeof( BASETYPE ) ), cudaMemcpyDeviceToHost ) );
}

//------------------------------------------------------------
}

//------------------------------------------------------------

#endif
