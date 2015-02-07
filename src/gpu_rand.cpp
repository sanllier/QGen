#include "gpu_rand.h"
#include "cuda_error_handler.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

GPURand::GPURand( unsigned long long seed, long long size/* = 0*/ )
    : m_randomBuf(0)
    , m_size(0)
{
    SAFE_CURAND_CALL( curandCreateGenerator( &m_gen, CURAND_RNG_PSEUDO_DEFAULT ) );
    SAFE_CURAND_CALL( curandSetPseudoRandomGeneratorSeed( m_gen, seed ) );
    if ( size > 0 )
        resize( size );
}

//------------------------------------------------------------

GPURand::~GPURand()
{
    SAFE_CURAND_CALL( curandDestroyGenerator( m_gen ) );
    SAFE_CALL( cudaFree( m_randomBuf ) );
}

//------------------------------------------------------------

void GPURand::generate()
{
   if ( m_size <= 0 )
       throw std::string( "GPURand trying to generate with unsetted size. " ).append( __FUNCTION__ );

    SAFE_CURAND_CALL( curandGenerateUniform( m_gen, m_randomBuf, m_size ));
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
