#include "qrotoperator.h"

#include <cmath>
#include <cstring>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

QRotOperator::QRotOperator( BASETYPE angle/* = BASETYPE(0)*/ )
{
    compute( angle );
}

//------------------------------------------------------------

void QRotOperator::compute( BASETYPE angle )
{
    m_matrix[0][0] = std::cos( angle );
    m_matrix[1][0] = std::sin( angle );

    m_matrix[1][1] =  m_matrix[0][0];    
    m_matrix[0][1] = -m_matrix[1][0];
}

//------------------------------------------------------------

QBit QRotOperator::operator*( const QBit& ampl ) const
{
    QBit res;
    res.a = m_matrix[0][0] * ampl.a + m_matrix[0][1] * ampl.b;
    res.b = m_matrix[1][0] * ampl.a + m_matrix[1][1] * ampl.b; 
    return res;
}

//------------------------------------------------------------
}
