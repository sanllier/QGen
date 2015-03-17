#ifndef SPARAMS_H
#define SPARAMS_H

#include "mpi.h"
#include "mpicheck.h"
#include "interfaces.h"
#include "defs.h"

#include <string>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

struct SParams
{
#ifdef GPU
    bool gpu;
#endif
    long long cycThreshold;
    int individsNum;
    int indSize;
    int topoRows;
    int topoCols;
    BASETYPE targetFitness;
    BASETYPE accuracy;
    IFitness*  fClass;
    IRepair*   repClass;
    IScreen* screenClass;
    std::string outFile;

    SParams( MPI_Comm comm = MPI_COMM_NULL, const char* file = 0, IFitness* fC = 0, IRepair* rC = 0, IScreen* sC = 0 );

    void initWithFile( MPI_Comm comm, const char* file, IFitness* fC, IRepair* rC = 0, IScreen* sC = 0 );
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
