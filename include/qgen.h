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
        , catastropheThreshold(0)
        , immigrationThreshold(0)
        , immigrationSize(0)
        , targetFitness( BASETYPE(0) )
        , accuracy( BASETYPE(0) )
        , fClass(0) 
        , repClass(0) {}
};

class QGenProcess
{
public:
    QGenProcess( const QGenProcessSettings& settings, MPI_Comm comm = MPI_COMM_WORLD );
    ~QGenProcess();

    void process();
    const QIndivid& getBestIndivid() const { return m_totalBest.ind; }

    inline bool isMaster() const { return m_myID == ROOT_ID; }

private:
    BASETYPE findIterationBestInd();
    bool immigration();

    inline bool active() const { return m_comm != MPI_COMM_NULL; }

private:
    static int m_instancesCount;
    MPI_Comm m_comm;
    int m_myID;
    int m_commSize;

    QGenProcessSettings m_settings;

    std::vector< QIndivid > m_individs;
    
    struct BestSolution
    {
        int procRank;
        int localIdx;
        QIndivid ind;

        BestSolution()
            : procRank(-1)
            , localIdx(-1)
        {}
    };
    BestSolution m_totalBest;
    BestSolution m_iterBest;
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
