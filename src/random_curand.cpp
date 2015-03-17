#ifdef CURAND
//------------------------------------------------------------

#include "random_curand.h"
#include "curand.h"
#include "cuda_error_handler.h"

//------------------------------------------------------------

RandomCURand::RandomCURand()
    : m_size(0)
    , m_bufPos(0)
    , m_gen(0)
    , m_randomBufGPU(0)
    , m_randomBufCPU(0)
{
    m_gen = new curandGenerator_t;
    SAFE_CURAND_CALL( curandCreateGenerator( ( curandGenerator_t * )m_gen, CURAND_RNG_PSEUDO_DEFAULT ) );
}

//------------------------------------------------------------

RandomCURand::~RandomCURand()
{
    SAFE_CURAND_CALL( curandDestroyGenerator( *( curandGenerator_t * )m_gen ) );
    SAFE_CALL( cudaFree( m_randomBufGPU ) );
    delete[] m_randomBufCPU;
    delete ( curandGenerator_t * )m_gen;
}

//------------------------------------------------------------

void RandomCURand::setSize( long long size )
{
   if ( size <= 0 )
       throw std::string( "GPURand trying to resize with incorrect parameter. " ).append( __FUNCTION__ );

    if ( m_randomBufGPU )
        SAFE_CALL( cudaFree( m_randomBufGPU ) );
    delete[] m_randomBufCPU;

    m_size = size;
    m_bufPos = size;
    SAFE_CALL( cudaMalloc( &m_randomBufGPU, size_t( m_size * sizeof( BASETYPE ) ) ) );
    m_randomBufCPU = new BASETYPE[ size_t( m_size ) ];
}

//------------------------------------------------------------

void RandomCURand::setSeed( unsigned seed ) 
{
    // TODO
}

//------------------------------------------------------------

BASETYPE RandomCURand::next()
{
   if ( m_size <= 0 )
       throw std::string( "GPURand trying to generate with unsetted size. " ).append( __FUNCTION__ );

   if (m_size == m_bufPos)
   {
       generate();
       m_bufPos = 0;
   }

   return m_randomBufCPU[ m_bufPos++ ];
}

//------------------------------------------------------------

void RandomCURand::generate() const
{
   if ( m_size <= 0 )
       throw std::string( "GPURand trying to generate with unsetted size. " ).append( __FUNCTION__ );

    SAFE_CURAND_CALL( curandGenerateUniform( *( curandGenerator_t * )m_gen, m_randomBufGPU, size_t( m_size ) ));
    SAFE_CALL( cudaMemcpy( m_randomBufCPU, m_randomBufGPU, size_t( m_size * sizeof( BASETYPE ) ), cudaMemcpyDeviceToHost ) );
}

//------------------------------------------------------------

#endif

//------------------------------------------------------------
