#ifndef QINDIVID_BASE_H
#define QINDIVID_BASE_H

#include "mpi.h"
#include "qbit.h"
#include "qobservestate.h"
#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

enum EIndividType
{
    INDIVID_TYPE_CPU,
    INDIVID_TYPE_GPU
};

//------------------------------------------------------------

class IFitness;
class IRepair;

//------------------------------------------------------------

class QBaseIndivid
{
public:
    QBaseIndivid( long long size, int coords[2], MPI_Comm individComm, MPI_Comm rowComm, BASETYPE thetaFrac = BASETYPE(1) );
    virtual ~QBaseIndivid();

    virtual EIndividType getType() const = 0;

    virtual bool resize( long long newSize );
    virtual void setInitial() = 0;
    
    inline long long qSize() const { return m_globalLogicSize; }
    inline long long locaQSize() const { return m_localLogicSize; }
    inline long long firstQBit() const { return m_firstQbit; }

    BASETYPE calculateFitness( IFitness* fClass );
    BASETYPE getFitness() const { return m_fitness; } // unsafe, m_fitness may be dirty

    void calculateObservState();
    QObserveState& getObservState();
    const QObserveState& getObservState() const;

    virtual void evolve( const QBaseIndivid& bestInd ) = 0;
    void repair( IRepair* repClass );

    virtual bool bcast( int root );

    virtual QBaseIndivid& operator=( const QBaseIndivid& rInd ) = 0;

protected:
    static MPI_Datatype MPI_QBIT;

    QBit* m_data;

    long long m_localLogicSize;
    long long m_globalLogicSize;
    long long m_firstQbit;

    QObserveState* m_observeState;

    BASETYPE m_fitness;
    bool m_needRecalcFitness;

    static BASETYPE m_thetaField[2][2][2];

    struct SBaseIndividContext
    {
        MPI_Comm indComm;
        MPI_Comm rowComm;
        int rowSize;
        int coords[2];

        SBaseIndividContext()
            : indComm( MPI_COMM_NULL )
            , rowComm( MPI_COMM_NULL )
            , rowSize(0)
        {}
    };
    SBaseIndividContext m_context;
};

//------------------------------------------------------------
}

//------------------------------------------------------------
#endif
