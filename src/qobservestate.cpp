#include "qobservestate.h"
#include "qindivid_cpu.h"
#include "sharedmtrand.h"
#include "mpicheck.h"
#include <string.h>
#ifdef GPU
    #include "qindivid_gpu.h"
    #include "cuda_runtime.h"
    #include "cuda_error_handler.h"
#ifdef CURAND
    #include <time.h>
#endif
#endif

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QObserveState::QObserveState()
    : m_state(0)
    , m_stateSize(0)
#ifdef GPU
    , m_gpuBuf(0)
#ifdef CURAND
    , m_gpuRand( time(0) )
    , m_randBuf(0)
#endif
#endif
{

}

//------------------------------------------------------------

QObserveState::~QObserveState()
{
    clear();
}

//------------------------------------------------------------

void QObserveState::clear()
{ 
    delete[] m_state;
    m_state = 0;

#ifdef GPU
    if ( m_gpuBuf )
        SAFE_CALL( cudaFree( m_gpuBuf ) );
    m_gpuBuf = 0;

#ifdef CURAND
    if ( m_randBuf )
        delete[] m_randBuf;
    m_randBuf = 0;
#endif

#endif
}

//------------------------------------------------------------

void QObserveState::observe( const QBaseIndivid& ind )
{
    const long long localStateSize = ind.locaQSize();
    resize( localStateSize );

    switch ( ind.getType() )
    {
        case INDIVID_TYPE_CPU:
        {
        #if defined( GPU ) && defined( CURAND )
            m_gpuRand.generate();
            m_gpuRand.copyInCPUData( m_randBuf );
        #endif
            for ( long long i = 0; i < localStateSize; ++i )
            {
                const QCPUIndivid& castedIndivid = ( const QCPUIndivid& )ind;
            #ifdef STDRAND
                const BASETYPE randVal = BASETYPE( rand() ) / RAND_MAX;
            #elif defined( GPU ) && defined( CURAND )
                const BASETYPE& randVal = m_randBuf[i];
            #else
                const BASETYPE randVal = BASETYPE( SharedMTRand::getClosedInstance()() );
            #endif
                const QBit& qbit = castedIndivid.localAt(i);
                const BASETYPE real = qbit.a.real();
                const BASETYPE imag = qbit.a.imag();
                const BASETYPE mod = std::abs( real * real + imag * imag );
                m_state[i] = randVal >= mod;
            }
            break;
        }

    #ifdef GPU
        case INDIVID_TYPE_GPU:
        {
            const QGPUIndivid& castedIndivid = ( const QGPUIndivid& )ind;            
            castedIndivid.runObserveKernel( m_gpuBuf );
            SAFE_CALL( cudaMemcpy( m_state, m_gpuBuf, size_t( m_stateSize ), cudaMemcpyDeviceToHost ) );
            break;
        }
    #endif

        default:
        {
            throw std::string( "QObserveState is trying to observe individ with unknown type. " ).append( __FUNCTION__ );  
            break;
        }
    }
}

//------------------------------------------------------------

void QObserveState::resize( long long size )
{
    if ( m_stateSize != size )
    {
        clear();
        m_state = new bool[ size_t(size) ];
#ifdef GPU
        SAFE_CALL( cudaMalloc( &m_gpuBuf, size_t(size) ) );
#ifdef CURAND
        m_gpuRand.resize( size );
        m_randBuf = new BASETYPE[ size_t(size) ];
#endif
#endif
        m_stateSize = size;
    }
}

//------------------------------------------------------------

bool QObserveState::at( long long pos ) const
{ 
    if ( pos < 0 || pos >= m_stateSize )
        throw std::string( "QObservState out of bounds. " ).append( __FUNCTION__ ); 

    return m_state[ pos ]; 
}

//------------------------------------------------------------

void QObserveState::set( long long pos, bool val )
{
    if ( pos < 0 || pos >= m_stateSize )
        throw std::string( "QObservState out of bounds. " ).append( __FUNCTION__ );

    m_state[ pos ] = val;
}

//------------------------------------------------------------

void QObserveState::bcast( int root, MPI_Comm comm )
{
    CHECK( MPI_Bcast( m_state, int(m_stateSize), MPI_CHAR, root, comm ) );
}

//------------------------------------------------------------

QObserveState& QObserveState::operator=( const QObserveState& rState )
{
    resize( rState.m_stateSize );
    memcpy( m_state, rState.m_state, size_t(m_stateSize) );
    return *this;
}

//------------------------------------------------------------
}
