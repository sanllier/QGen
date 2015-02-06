#if !defined( QINDIVID_GPU_H ) && defined( GPU )
#define QINDIVID_GPU_H

#include "qindivid_base.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

class QGPUIndivid: public QBaseIndivid
{
public:
    QGPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] );
    ~QGPUIndivid() override;

        // QBaseIndivid

    EIndividType getType() const override { return INDIVID_TYPE_GPU; };

    void resize( long long newSize ) override;
    void setInitial() override;
    
    void evolve( const QBaseIndivid& bestInd ) override;

    void bcast( int root ) override;

    QBaseIndivid& operator=( const QBaseIndivid& rInd ) override;

protected:
    BASETYPE getThetaForQBit( const QBaseIndivid& bestInd, long long qbitIndex ) const override;

};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
