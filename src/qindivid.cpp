#include "qindivid.h"

#include <string>
#include <cstring>
#include <cmath>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QIndivid::QIndivid( size_t size/* = 0*/ )
{
    if ( size < 0 )
        throw std::string( "Invalid individ size. " ).append( __FUNCTION__ );

    resize( size );
    setInitial();
}

QIndivid::QIndivid( const QIndivid& ind )
{
    m_data = new QBit[ ind.m_dataLogicSize ];
    m_dataLogicSize = ind.m_dataLogicSize;

    std::memcpy( m_data, ind.m_data, m_dataLogicSize * sizeof( QBit ) );
}

QIndivid::~QIndivid()
{
    delete[] m_data;
}

void QIndivid::resize( size_t size )
{
    m_data = new QBit[ size ];
    m_dataLogicSize = size;
}

void QIndivid::setInitial()
{
    const BASETYPE defVal = (BASETYPE)( 1.0 / std::sqrt( 2 ) );
    for ( size_t i = 0; i < m_dataLogicSize; ++i )
    {
        m_data[i].a = defVal;
        m_data[i].b = defVal;
    }
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

    return *this;
}

//-----------------------------------------------------------

inline BASETYPE func( const QBit& myQbit, const QBit& bestIndQbit )
{
    const BASETYPE dBest = bestIndQbit.a * bestIndQbit.b;
    const BASETYPE dMy   = myQbit.a * myQbit.b;
    const float sigmaBest = std::fabsf( std::atan2( bestIndQbit.b, bestIndQbit.a ) );
    const float sigmaMy   = std::fabsf( std::atan2( myQbit.b, myQbit.a ) );

    if ( ( dBest > 0 && dMy > 0 ) || ( dBest <= 0 && dMy <= 0 ) )
        return sigmaBest > sigmaMy ? +1.0f : -1.0f;
    else
        return sigmaBest > sigmaMy ? -1.0f : +1.0f;

    return 0;
}

void QIndivid::tick( const QIndivid& globalBestInd, BASETYPE factor )
{
    QBit temp;
    QRotOperator op;

    for ( size_t i = 0; i < m_dataLogicSize; ++i )
    {
        op.compute( factor * func( m_data[i], globalBestInd.m_data[i] ) );

        temp.a = m_data[i].a * op[0][0] + m_data[i].b * op[0][1]; 
        temp.b = m_data[i].a * op[1][0] + m_data[i].b * op[1][1]; 
        m_data[i] = temp;
    }
}

//------------------------------------------------------------
}
