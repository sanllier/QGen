#ifndef QGEN_H
#define QGEN_H

#include <complex>
#include <vector>

#include "mpi.h"
#include "mpicheck.h"
#include "qindivid.h"
#include "qobservstate.h"
#include "defs.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------


class QProcessScreen
{
public:
    virtual void operator()( long long cycle, 
                             const int coords[2], 
                       const QIndivid& totalBest, 
                       const QIndivid& iterBest ) = 0;
};

//------------------------------------------------------------

struct QGenProcessSettings
{
    long long cycThreshold;
    int individsNum;
    int indSize;
    int topoRows;
    int topoCols;
    BASETYPE targetFitness;
    BASETYPE accuracy;
    QFitnessClass* fClass;
    QRepairClass*  repClass;
    QProcessScreen* screenClass;

    QGenProcessSettings()
        : cycThreshold(0)
        , individsNum(0)
        , indSize(0)
        , topoRows(1)
        , topoCols(1)
        , targetFitness( BASETYPE(0) )
        , accuracy( BASETYPE(0) )
        , fClass(0) 
        , repClass(0)
        , screenClass(0)
    {}
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
    const QIndivid* getBestIndivid() const { return m_totalBest->ind; }

    inline bool isMaster() const { return m_ctx.generalRank == ROOT_ID; }
    inline bool isMasterInd() const { return m_ctx.coords[1] == 0; }

private:
    BASETYPE findIterationBestInd();

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
        MPI_Comm rowComm;
        QIndivid* ind;

        BestSolution()
            : procRank(-1)
            , localIdx(-1)
            , rowComm( MPI_COMM_NULL )
            , ind(0)
        {}
        ~BestSolution()
        {
            delete ind;
            //FIXME: CHECK( MPI_Comm_free( &rowComm ) );
        }

        BestSolution& operator=( const BestSolution& rSol )
        {
            procRank = rSol.procRank;
            localIdx = rSol.localIdx;
            *ind = *rSol.ind;

            return *this;
        }
    };
    BestSolution* m_totalBest;
    BestSolution* m_iterBest;
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
