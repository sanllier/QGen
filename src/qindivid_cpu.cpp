#include "qindivid_cpu.h"
#include "qgen.h"
#include "qrotoperator.h"
#include "sharedmtrand.h"
#include "mpicheck.h"
#ifdef GPU
    #include "qindivid_gpu.h"
    #include "cuda_runtime.h"
    #include "cuda_error_handler.h"
#endif

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QCPUIndivid::QCPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] )
    : QBaseIndivid( size, generalComm,  rowComm, coords )
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

void QCPUIndivid::resize( long long newSize )
{
    QBaseIndivid::resize( newSize );

    delete[] m_data;
    m_data = new QBit[ size_t( m_localLogicSize ) ];
}

//------------------------------------------------------------

void QCPUIndivid::setInitial()
{
    for ( long long i = 0; i < m_localLogicSize; ++i )
    {
        m_data[i].a = ( BASETYPE )( SharedMTRand::getClosedInstance()() );
        m_data[i].b = std::sqrt( QComplex( BASETYPE(1) ) - m_data[i].a );
    }

    m_needRecalcFitness = true;
}

//------------------------------------------------------------

void QCPUIndivid::evolve( const QBaseIndivid& bestInd )
{
    QRotOperator op;
    for ( long long i = 0; i < m_localLogicSize; ++i )
    {
        op.compute( getThetaForQBit( bestInd, i ) );
        m_data[i] = op * m_data[i];
    }
}

//------------------------------------------------------------

void QCPUIndivid::bcast( int root )
{
    QBaseIndivid::bcast( root );

    if ( m_context.rowComm != MPI_COMM_NULL )
        CHECK( MPI_Bcast( m_data, int( m_localLogicSize ), MPI_QBIT, root, m_context.rowComm ) );
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
            std::memcpy( m_data, castedIndivid.m_data, size_t( m_localLogicSize * sizeof( m_data[0] ) ) );
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

    m_observeState      = castedIndivid.m_observeState;
    m_fitness           = castedIndivid.m_fitness;
    m_needRecalcFitness = castedIndivid.m_needRecalcFitness;

    return *this;
}

//-----------------------------------------------------------

BASETYPE QCPUIndivid::getThetaForQBit( const QBaseIndivid& bestInd, long long qbitIndex ) const
{
    const static BASETYPE PI = BASETYPE( 3.14159265359 );
    const static BASETYPE thetaField[2][2][2] = { { { BASETYPE( 0.01 ) * PI, BASETYPE( 0.01 ) * PI }, 
                                                    { BASETYPE( 0.8 )  * PI, BASETYPE( 0.01 ) * PI } }, 
                                                  { { BASETYPE( 0.8 )  * PI, BASETYPE( 0.01 ) * PI }, 
                                                    { BASETYPE( 0.01 ) * PI, BASETYPE( 0.01 ) * PI } } };

    const QCPUIndivid& bestIndividCasted = ( const QCPUIndivid& )bestInd;
    const bool curIndBit      = m_observeState.at( qbitIndex );
    const bool bestIndBit     = bestIndividCasted.m_observeState.at( qbitIndex );
    const bool betterThanBest = m_fitness > bestIndividCasted.m_fitness;

    BASETYPE theta = thetaField[ curIndBit ? 1:0 ][ bestIndBit ? 1:0 ][ betterThanBest ? 1:0 ];
    const QBit curIndQBit = m_data[ qbitIndex ];
    BASETYPE prodRealPart = std::real( curIndQBit.a * curIndQBit.b );

    if ( prodRealPart > 0 )
        return curIndBit ? theta : -theta;        
    else if ( prodRealPart < 0 )
        return curIndBit ? -theta : theta;

    return theta;
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
