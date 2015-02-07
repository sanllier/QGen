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
    virtual ~QRotOperator() {}

    void compute( BASETYPE angle );

    QBit operator*( const QBit& ampl ) const;

private:
    BASETYPE m_matrix[2][2];
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
