#include "qobservestate.h"
#include "mpicheck.h"
#include "qindivid_cpu.h"
#include <string.h>
#include <sstream>
#ifdef GPU
    #include "qindivid_gpu.h"
    #include "cuda_runtime.h"
    #include "cuda_error_handler.h"
#endif

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QObserveState::QObserveState()
    : m_state(0)
    , m_stateSize(0)
    , m_rand(0)
#ifdef GPU
    , m_gpuBuf(0)
#endif
{
    m_rand = Randomizer::Get();
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
            for ( long long i = 0; i < localStateSize; ++i )
            {
                const QCPUIndivid& castedIndivid = ( const QCPUIndivid& )ind;
                const BASETYPE randVal = m_rand->next();

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
#endif
#ifdef CURAND
        ( (RandomCURand *)m_rand )->setSize( size );
#endif
        m_stateSize = size;
    }
}

//------------------------------------------------------------

bool QObserveState::at( long long pos ) const
{ 
    if ( pos < 0 || pos >= m_stateSize )
    {
        std::stringstream sStr;
        sStr << "QObservState out of bounds. " << pos << "/" << m_stateSize << " " << __FUNCTION__;
        throw sStr.str();
    }
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
