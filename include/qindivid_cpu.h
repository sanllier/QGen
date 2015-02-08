#ifndef QINDIVID_CPU_H
#define QINDIVID_CPU_H

#include "qindivid_base.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QCPUIndivid: public QBaseIndivid
{
public:
    QCPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] );
    ~QCPUIndivid() override;

        // QBaseIndivid

    EIndividType getType() const override { return INDIVID_TYPE_CPU; };

    bool resize( long long newSize ) override;
    void setInitial() override;
    
    void evolve( const QBaseIndivid& bestInd ) override;

    bool bcast( int root ) override;

    QBaseIndivid& operator=( const QBaseIndivid& rInd ) override;

private:
    BASETYPE getThetaForQBit( const QCPUIndivid& bestInd, long long qbitIndex ) const;

public:
    const QBit& localAt( long long pos ) const;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
