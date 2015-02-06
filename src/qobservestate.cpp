#include "qobservestate.h"
#include "qindivid_cpu.h"
#include "qindivid_gpu.h"
#include "sharedmtrand.h"
#include "mpicheck.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QObserveState::QObserveState()
    : m_state(0)
    , m_stateSize(0)
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
            #ifdef STDRAND
                const BASETYPE randVal = BASETYPE( rand() ) / RAND_MAX;
            #else
                const BASETYPE randVal = BASETYPE( SharedMTRand::getClosedInstance()() );
            #endif
                const BASETYPE mod = std::abs( castedIndivid.localAt(i).a );
                m_state[i] = randVal >= mod * mod;
            }

            break;
        }

    #ifdef GPU
        case INDIVID_TYPE_GPU:
        {
            // GPU code here
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
    std::memcpy( m_state, rState.m_state, size_t(m_stateSize) );
    return *this;
}

//------------------------------------------------------------
}
