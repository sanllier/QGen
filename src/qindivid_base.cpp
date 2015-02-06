#include "qindivid_base.h"
#include "mpicheck.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

MPI_Datatype QBaseIndivid::MPI_QBIT = MPI_DATATYPE_NULL;

//------------------------------------------------------------

QBaseIndivid::QBaseIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] ) 
    : m_fitness( BASETYPE(0) )
    , m_needRecalcFitness( true )
{
    if ( size < 0 )
        throw std::string( "Invalid individ size. " ).append( __FUNCTION__ );
    if ( generalComm == MPI_COMM_NULL )
        throw std::string( "Invalid general communicator. " ).append( __FUNCTION__ );

    std::memcpy( m_context.coords, coords, 2 * sizeof( coords[0] ) );
    
    int remainIndDims[2] = { 1, 0 };
    CHECK( MPI_Cart_sub( generalComm, remainIndDims, &m_context.indComm ) ); 

    if ( MPI_QBIT == MPI_DATATYPE_NULL )
    {
        CHECK( MPI_Type_vector( 1, 4, sizeof( BASETYPE ), MPI_BASETYPE, &MPI_QBIT ) );
        CHECK( MPI_Type_commit( &MPI_QBIT ) );
    }
}

//------------------------------------------------------------

QBaseIndivid::~QBaseIndivid()
{

}

//------------------------------------------------------------

void QBaseIndivid::resize( long long newSize )
{
    int commSize = 0;
    CHECK( MPI_Comm_size( m_context.indComm, &commSize ) );

    const int procGridHeight = m_context.coords[0];

    m_globalLogicSize = newSize;
    m_localLogicSize = newSize / commSize + ( procGridHeight < newSize % commSize ? 1 : 0 );
    m_firstQbit = ( newSize / commSize ) * procGridHeight + std::min( newSize % commSize, ( long long )procGridHeight ); 

    // CRAP
    //m_observeState.requestMemory( m_localLogicSize );
    m_needRecalcFitness = true;
}

//------------------------------------------------------------

BASETYPE QBaseIndivid::calculateFitness( QFitnessClass* fClass )
{
    if ( !fClass )
        throw std::string( "QXIndivid is trying to calculate fitness with (NULL) func" ).append( __FUNCTION__ );  

    bool needRecalc = false;
    CHECK( MPI_Allreduce( &m_needRecalcFitness, &needRecalc, 1, MPI_CHAR, MPI_SUM, m_context.indComm ) );

    if ( needRecalc )
    {
        m_fitness = (*fClass)( m_context.indComm, m_observeState, m_firstQbit, m_context.coords[0] );
        m_needRecalcFitness = false;
    }

    return m_fitness;
}

//------------------------------------------------------------

void QBaseIndivid::calculateObservState()
{
    m_observeState.observe( *this );
    m_needRecalcFitness = true;
}

//------------------------------------------------------------

QObserveState& QBaseIndivid::getObservState()
{
    m_needRecalcFitness = true;
    return m_observeState;
}

//------------------------------------------------------------

const QObserveState& QBaseIndivid::getObservState() const
{
    return m_observeState;
}

//-----------------------------------------------------------

void QBaseIndivid::repair( QRepairClass* repClass )
{
    if ( !repClass )
        throw std::string( "QXIndivid is trying to repair with (NULL) func" ).append( __FUNCTION__ );

    (*repClass)( m_context.indComm, m_observeState, m_firstQbit, m_context.coords[0] );
    m_needRecalcFitness = true;
}

//------------------------------------------------------------

void QBaseIndivid::bcast( int root )
{
    if ( m_context.rowComm == MPI_COMM_NULL )
        return;

    int rowCommSize = 0;
    CHECK( MPI_Comm_size( m_context.rowComm, &rowCommSize ) );
    if ( root < 0 || root >= rowCommSize )
        throw std::string( "QXIndivid is trying to bcast with invalid params" ).append( __FUNCTION__ ); 
    
    m_observeState.bcast( root, m_context.rowComm );
    CHECK( MPI_Bcast( &m_fitness, 1, MPI_BASETYPE, root, m_context.rowComm ) );
    CHECK( MPI_Bcast( &m_needRecalcFitness, 1, MPI_CHAR, root, m_context.rowComm ) );    
}

//------------------------------------------------------------
}
