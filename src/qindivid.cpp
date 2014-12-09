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

QIndivid::QIndivid( size_t size/* = 0*/ )
    : m_data(0)
    , m_fitness( BASETYPE(0) )
    , m_needRecalcFitness( true )
{
    if ( size < 0 )
        throw std::string( "Invalid individ size. " ).append( __FUNCTION__ );

    m_data = new QBit[ size ];
    m_dataLogicSize = size;

    setInitial();

    m_obsState.requestMemory( size );
}

QIndivid::QIndivid( const QIndivid& ind )    
    : m_data(0)
    , m_fitness( BASETYPE(0) )
    , m_needRecalcFitness( true )
{
    m_data = new QBit[ ind.m_dataLogicSize ];
    m_dataLogicSize = ind.m_dataLogicSize;
    std::memcpy( m_data, ind.m_data, m_dataLogicSize * sizeof( QBit ) );

    m_obsState = ind.m_obsState;
}

QIndivid::~QIndivid()
{
    delete[] m_data;
}

//------------------------------------------------------------

void QIndivid::resize( size_t size )
{
    delete[] m_data;
    m_data = new QBit[ size ];
    m_dataLogicSize = size;

    m_obsState.requestMemory( size );

    m_needRecalcFitness = true;
}

void QIndivid::setInitial()
{
    for ( size_t i = 0; i < m_dataLogicSize; ++i )
    {
        m_data[i].a = ( BASETYPE )( SharedMTRand::getClosedInstance()() );
        m_data[i].b = std::sqrt( QComplex( BASETYPE(1) ) - m_data[i].a );
    }

    m_needRecalcFitness = true;
}

QIndivid& QIndivid::operator=( const QIndivid& ind )
{
    if ( ind.qsize() != this->qsize() )
    {
        delete[] m_data;
        m_data = new QBit[ ind.m_dataLogicSize ];
        m_dataLogicSize = ind.m_dataLogicSize;
    }
    std::memcpy( m_data, ind.m_data, m_dataLogicSize * sizeof( QBit ) );
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
    const QBit curIndQBit = curInd.at( qbitIdx );
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
    for ( size_t i = 0; i < m_dataLogicSize; ++i )
    {
        op.compute( getThetaForQBit( *this, bestInd, i ) );
        m_data[i] = op * m_data[i];
    }
}

//------------------------------------------------------------

void QIndivid::repair( QRepairClass* repClass )
{
    if ( m_needRecalcFitness && !repClass )
        throw std::string( "QIndivid is trying to repair with (NULL) func" ).append( __FUNCTION__ ); 

    (*repClass)( m_obsState );
    m_needRecalcFitness = true;
}

//------------------------------------------------------------

void QIndivid::bcast( int root, MPI_Comm comm )
{
    if ( root < 0 || comm ==  MPI_COMM_NULL )
        throw std::string( "QIndivid is trying to bcast with invalid params" ).append( __FUNCTION__ ); 

    //CHECK( MPI_Bcast( m_data, m_dataLogicSize, QGenProcess::getQbitType(), root, comm ) );
    CHECK( MPI_Bcast( m_obsState.data(), m_obsState.size(), MPI_C_BOOL, root, comm ) );
    CHECK( MPI_Bcast( &m_fitness, 1, MPI_BASETYPE, root, comm ) );
    m_needRecalcFitness = false;
}

//------------------------------------------------------------
}
