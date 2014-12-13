#ifndef QGEN_H
#define QGEN_H

#include <complex>
#include <vector>

#include "mpi.h"

#include "qindivid.h"
#include "qobservstate.h"
#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

struct QGenProcessSettings
{
    long long cycThreshold;
    double timeThreshold;
    int individsNum;
    int indSize;
    int topoRows;
    int topoCols;
    long long catastropheThreshold;
    long long immigrationThreshold;
    long long immigrationSize;
    BASETYPE targetFitness;
    BASETYPE accuracy;
    QFitnessClass* fClass;
    QRepairClass*  repClass;

    QGenProcessSettings()
        : cycThreshold(0)
        , timeThreshold( 0.0 )
        , individsNum(0)
        , indSize(0)
        , topoRows(1)
        , topoCols(1)
        , catastropheThreshold(0)
        , immigrationThreshold(0)
        , immigrationSize(0)
        , targetFitness( BASETYPE(0) )
        , accuracy( BASETYPE(0) )
        , fClass(0) 
        , repClass(0) {}
};

//------------------------------------------------------------

struct QGenProcessContext
{
    MPI_Comm generalComm;
    MPI_Comm rowComm;

    int generalRank;
    int rowRank;

    int generalSize;
    int rowSize;

    int coords[2];

    QGenProcessContext()
        : generalComm( MPI_COMM_NULL )
        , rowComm( MPI_COMM_NULL )
        , generalRank( -1 )
        , rowRank( -1 )
        , generalSize(0)
        , rowSize(0)
    {}
};

//------------------------------------------------------------

class QGenProcess
{
public:
    QGenProcess( const QGenProcessSettings& settings, MPI_Comm comm = MPI_COMM_WORLD );
    ~QGenProcess();

    double process();
    const QIndivid& getBestIndivid() const { return *m_totalBest.ind; }

    inline bool isMaster() const { return m_ctx.generalRank == ROOT_ID; }

    inline static MPI_Datatype getQbitType() { return MPI_QBIT; }

private:
    BASETYPE findIterationBestInd();
    bool immigration();

    inline bool active() const { return m_ctx.generalComm != MPI_COMM_NULL; }

private:
    static int m_instancesCount;

    QGenProcessContext m_ctx;
    QGenProcessSettings m_settings;

    std::vector< QIndivid* > m_individs;
    
    struct BestSolution
    {
        int procRank;
        int localIdx;
        QIndivid* ind;

        BestSolution()
            : procRank(-1)
            , localIdx(-1)
            , ind(0)
        {}
    };
    BestSolution m_totalBest;
    BestSolution m_iterBest;

    static MPI_Datatype MPI_QBIT;
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
