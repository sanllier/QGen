#ifndef QINDIVID_CPU_H
#define QINDIVID_CPU_H

#include "qindivid_base.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QCPUIndivid: public QBaseIndivid
{
public:
    QCPUIndivid( long long size, int coords[2], MPI_Comm comm );
    ~QCPUIndivid() OVERRIDE;

        // QBaseIndivid

    EIndividType getType() const OVERRIDE { return INDIVID_TYPE_CPU; };

    bool resize( long long newSize ) OVERRIDE;
    void setInitial() OVERRIDE;
    
    void evolve( const QBaseIndivid& bestInd ) OVERRIDE;

    bool bcast( int root ) OVERRIDE;

    QBaseIndivid& operator=( const QBaseIndivid& rInd ) OVERRIDE;

private:
    BASETYPE getThetaForQBit( const QCPUIndivid& bestInd, long long qbitIndex ) const;

public:
    const QBit& localAt( long long pos ) const;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
