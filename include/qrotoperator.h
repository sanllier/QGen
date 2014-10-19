#ifndef QROTOPERATOR_H
#define QROTOPERATOR_H

#include <complex>

#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QRotOperator
{
public:
    QRotOperator( BASETYPE angle = BASETYPE(0) );
    QRotOperator( const QRotOperator& op );
    virtual ~QRotOperator() {}

    void compute( BASETYPE angle );
    inline const BASETYPE* operator[]( size_t row ) const { return m_matrix[ row ]; } 

private:
    BASETYPE m_matrix[2][2];
    BASETYPE m_curAngle;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
