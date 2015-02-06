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

    void resize( long long newSize ) override;
    void setInitial() override;
    
    void evolve( const QBaseIndivid& bestInd ) override;

    void bcast( int root ) override;

    QBaseIndivid& operator=( const QBaseIndivid& rInd ) override;

protected:
    BASETYPE getThetaForQBit( const QBaseIndivid& bestInd, long long qbitIndex ) const override;

public:
    const QBit& localAt( long long pos ) const;

private:
    QBit* m_data;

};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
