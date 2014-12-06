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

QRotOperator::QRotOperator( const QRotOperator& op )
{
    std::memcpy( m_matrix[0], op.m_matrix[0], 2 * sizeof( BASETYPE ) );
    std::memcpy( m_matrix[1], op.m_matrix[1], 2 * sizeof( BASETYPE ) );
}

void QRotOperator::compute( BASETYPE angle )
{
    m_matrix[0][0] = std::cos( angle );
    m_matrix[1][0] = std::sin( angle );

    m_matrix[1][1] =  m_matrix[0][0];    
    m_matrix[0][1] = -m_matrix[1][0];
}

//------------------------------------------------------------
}
