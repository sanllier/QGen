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
    inline QBit operator*( const QBit& ampl ) const
    {
        QBit res;
        res.a = m_matrix[0][0] * ampl.a + m_matrix[0][1] * ampl.b;
        res.b = m_matrix[1][0] * ampl.a + m_matrix[1][1] * ampl.b; 
        return res;
    }

private:
    BASETYPE m_matrix[2][2];
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
