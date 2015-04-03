#include "qindivid_cpu.h"
#include "qgen.h"
#include "qrotoperator.h"
#include "mpicheck.h"
#include <string.h>
#include <time.h>
#ifdef GPU
    #include "qindivid_gpu.h"
    #include "cuda_runtime.h"
    #include "cuda_error_handler.h"
#endif

#if defined( CURAND )
    #include "random_curand.h"
#elseif defined( MTRAND )
    #include "random_mtrand.h"
#else 
    #include "random_def.h"   
#endif

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QCPUIndivid::QCPUIndivid( long long size, int coords[2], MPI_Comm comm )
    : QBaseIndivid( size, coords, comm )
{
    resize( size );
    setInitial();
}

//------------------------------------------------------------

QCPUIndivid::~QCPUIndivid()
{
    delete[] m_data;
}

//------------------------------------------------------------

bool QCPUIndivid::resize( long long newSize )
{
    if ( !QBaseIndivid::resize( newSize ) )
        return false;

    delete[] m_data;
    m_data = new QBit[ size_t( m_localLogicSize ) ];
    return true;
}

//------------------------------------------------------------

void QCPUIndivid::setInitial()
{
    IRandom* random = 0;
#if defined( CURAND )
    random = new RandomCURand();
    ( (RandomCURand*)random )->setSize( m_localLogicSize );
#elseif defined( MTRAND )
    random = new RandomMTRand();
#else 
    random = new RandomDefault();
#endif

    random->setSeed( (unsigned)time(0) ^ ( m_context.coords[0] + m_context.coords[1] ) );  

    for ( long long i = 0; i < m_localLogicSize; ++i )
    {
        m_data[i].a = random->next();
        m_data[i].b = std::sqrt( QComplex( BASETYPE(1) ) - m_data[i].a );
    }

    m_needRecalcFitness = true;
    delete random;
}

//------------------------------------------------------------

void QCPUIndivid::evolve( const QBaseIndivid& bestInd )
{
    QRotOperator op;
    const QCPUIndivid& castedIndivid = ( const QCPUIndivid& )bestInd;
    for ( long long i = 0; i < m_localLogicSize; ++i )
    {
        op.compute( getThetaForQBit( castedIndivid, i ) );
        m_data[i] = op * m_data[i];
    }
}

//------------------------------------------------------------

bool QCPUIndivid::bcast( int root )
{
    if ( !QBaseIndivid::bcast( root ) )
        return false;

    CHECK( MPI_Bcast( m_data, int( m_localLogicSize ), MPI_QBIT, root, m_context.rowComm ) );
    return true;
}

//-----------------------------------------------------------

QBaseIndivid& QCPUIndivid::operator=( const QBaseIndivid& rInd )
{
    // CRAP
    const QCPUIndivid& castedIndivid = ( const QCPUIndivid& )rInd;
    if ( m_localLogicSize != castedIndivid.m_localLogicSize )
        throw std::string( "Trying to assing individs with different topologies. " ).append( __FUNCTION__ );

    switch ( rInd.getType() )
    {
        case INDIVID_TYPE_CPU:
        {
            memcpy( m_data, castedIndivid.m_data, size_t( m_localLogicSize * sizeof( m_data[0] ) ) );
            break;
        }

    #ifdef GPU
        case INDIVID_TYPE_GPU:
        {
            SAFE_CALL( cudaMemcpy( m_data, castedIndivid.m_data, 
                size_t( m_localLogicSize * sizeof( m_data[0] ) ), cudaMemcpyDeviceToHost ) );
            break;
        }
    #endif

        default:
        {
            throw std::string( "QCPUIndivid is trying to assign individ with unknown type. " ).append( __FUNCTION__ );  
            break;
        }
    }

    *m_observeState     = *castedIndivid.m_observeState;
    m_fitness           = castedIndivid.m_fitness;
    m_needRecalcFitness = castedIndivid.m_needRecalcFitness;

    return *this;
}

//-----------------------------------------------------------

BASETYPE QCPUIndivid::getThetaForQBit( const QCPUIndivid& bestInd, long long qbitIndex ) const
{
    const bool curIndBit      = m_observeState->at( qbitIndex );
    const bool bestIndBit     = bestInd.m_observeState->at( qbitIndex );
    const bool betterThanBest = m_fitness > bestInd.m_fitness;

    BASETYPE theta = m_thetaField[ curIndBit ? 1:0 ][ bestIndBit ? 1:0 ][ betterThanBest ? 1:0 ];
    const QBit curIndQBit = m_data[ qbitIndex ];
    BASETYPE prodRealPart = std::real( curIndQBit.a * curIndQBit.b );

    if ( prodRealPart >= 0 )
        return curIndBit ? theta : -theta;        

    return curIndBit ? -theta : theta;
}

//-----------------------------------------------------------

const QBit& QCPUIndivid::localAt( long long pos ) const
{
    if ( pos < 0 || pos >= m_localLogicSize ) 
        throw std::string( "QXIndivid out of bounds. " ).append( __FUNCTION__ ); 

    return m_data[ pos ]; 
}

//------------------------------------------------------------
}
