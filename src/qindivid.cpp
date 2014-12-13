#include "qindivid.h"
#include "qgen.h"
#include "qobservstate.h"
#include "sharedmtrand.h"
#include "mpicheck.h"

#include <string>
#include <cstring>
#include <cmath>
#include <iostream>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QIndivid::QIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] )
    : m_data(0)
    , m_fitness( BASETYPE(0) )
    , m_rowComm( rowComm )
{
    if ( size < 0 )
        throw std::string( "Invalid individ size. " ).append( __FUNCTION__ );
    if ( generalComm == MPI_COMM_NULL )
        throw std::string( "Invalid general communicator. " ).append( __FUNCTION__ );

    std::memcpy( m_coords, coords, 2 * sizeof( *coords ) );
    
    int remainIndDims[2] = { 1, 0 };
    CHECK( MPI_Cart_sub( generalComm, remainIndDims, &m_indComm ) ); 

    resize( size );
    setInitial();
}

QIndivid::~QIndivid()
{
    delete[] m_data;
}

//------------------------------------------------------------

void QIndivid::resize( long long newSize )
{
    int commSize = 0;
    CHECK( MPI_Comm_size( m_indComm, &commSize ) );

    m_globalLogicSize = newSize;
    m_localLogicSize = newSize / commSize + ( m_coords[0] < newSize % commSize ? 1 : 0 );
    delete[] m_data;
    m_data = new QBit[ m_localLogicSize ];

    m_startQbit = ( newSize / commSize ) * m_coords[0] + std::min( newSize % commSize, ( long long )m_coords[0] ); 

    m_obsState.requestMemory( m_localLogicSize );
    m_needRecalcFitness = true;
}

void QIndivid::setInitial()
{
    for ( size_t i = 0; i < m_localLogicSize; ++i )
    {
        m_data[i].a = ( BASETYPE )( SharedMTRand::getClosedInstance()() );
        m_data[i].b = std::sqrt( QComplex( BASETYPE(1) ) - m_data[i].a );
    }

    m_needRecalcFitness = true;
}

QIndivid& QIndivid::operator=( const QIndivid& ind )
{
    if ( m_localLogicSize != ind.m_localLogicSize )
        throw std::string( "Different topologies. " ).append( __FUNCTION__ );

    std::memcpy( m_data, ind.m_data, m_localLogicSize * sizeof( *m_data ) );
    m_obsState = ind.m_obsState;
    m_fitness = ind.m_fitness;
    m_needRecalcFitness = ind.m_needRecalcFitness;

    return *this;
}

//-----------------------------------------------------------

BASETYPE getThetaForQBit( const QIndivid& curInd, const QIndivid& bestInd, size_t qbitIdx )
{
    const static BASETYPE PI = BASETYPE( 3.14159265359 );
    const static BASETYPE thetaField[2][2][2] = { { { BASETYPE( 0.01 ) * PI, BASETYPE( 0.01 ) * PI }, 
                                                    { BASETYPE( 0.8 )  * PI, BASETYPE( 0.01 ) * PI } }, 
                                                  { { BASETYPE( 0.8 )  * PI, BASETYPE( 0.01 ) * PI }, 
                                                    { BASETYPE( 0.01 ) * PI, BASETYPE( 0.01 ) * PI } } };

    const bool curIndBit  = curInd.getObsState().at( qbitIdx );
    const bool bestIndBit = bestInd.getObsState().at( qbitIdx );
    const bool betterThanBest = curInd.getFitnessUNSAFE() > bestInd.getFitnessUNSAFE();

    BASETYPE theta = thetaField[ curIndBit ? 1:0 ][ bestIndBit ? 1:0 ][ betterThanBest ? 1:0 ];
    const QBit curIndQBit = curInd.localAt( qbitIdx );
    BASETYPE prodRealPart = std::real( curIndQBit.a * curIndQBit.b );

    if ( prodRealPart > 0 )
        return curIndBit ? theta : -theta;        
    else if ( prodRealPart < 0 )
        return curIndBit ? -theta : theta;

    return theta;
}

void QIndivid::tick( const QIndivid& bestInd )
{
    QRotOperator op;
    for ( size_t i = 0; i < m_localLogicSize; ++i )
    {
        op.compute( getThetaForQBit( *this, bestInd, i ) );
        m_data[i] = op * m_data[i];
    }
}

//------------------------------------------------------------

BASETYPE QIndivid::getFitness( QFitnessClass* fClass )
{
    bool needRecalc = false;
    CHECK( MPI_Allreduce( &m_needRecalcFitness, &needRecalc, 1, MPI_CHAR, MPI_SUM, m_indComm ) );

    if ( needRecalc && !fClass )
        throw std::string( "QIndivid is trying to recalc fitness with (NULL) func" ).append( __FUNCTION__ );  

    if ( needRecalc )
    {
        m_fitness = (*fClass)( m_indComm, m_obsState, m_startQbit, m_coords[0] );
        m_needRecalcFitness = false;
    }

    return m_fitness;
}

//------------------------------------------------------------

void QIndivid::repair( QRepairClass* repClass )
{
    if ( !repClass )
        throw std::string( "QIndivid is trying to repair with (NULL) func" ).append( __FUNCTION__ ); 

    (*repClass)( m_indComm, m_obsState, m_startQbit, m_coords[0] );
    m_needRecalcFitness = true;
}

//------------------------------------------------------------

void QIndivid::bcast( int rootInd )
{
    if ( m_rowComm == MPI_COMM_NULL )
        return;

    static int rowCommSize = -1;
    if ( rowCommSize == -1 )
    {
        CHECK( MPI_Comm_size( m_rowComm, &rowCommSize ) );
        if ( rootInd < 0 || rootInd >= rowCommSize )
            throw std::string( "QIndivid is trying to bcast with invalid params" ).append( __FUNCTION__ ); 
    }

    CHECK( MPI_Bcast( m_data, m_localLogicSize, QGenProcess::getQbitType(), rootInd, m_rowComm ) );
    CHECK( MPI_Bcast( m_obsState.data(), m_obsState.size(), MPI_C_BOOL, rootInd, m_rowComm ) );
    CHECK( MPI_Bcast( &m_fitness, 1, MPI_BASETYPE, rootInd, m_rowComm ) );
    m_needRecalcFitness = false;
}

//------------------------------------------------------------
}
