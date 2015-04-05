#ifndef SPARAMS_H
#define SPARAMS_H

#include "mpi.h"
#include "mpicheck.h"
#include "interfaces.h"
#include "defs.h"

#include <map>
#include <string>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

struct SParams
{
#ifdef GPU
    bool gpu;
#endif
    int problemSize;
    long long cycThreshold;
    int individsNum;
    int indSize;
    int topoRows;
    int topoCols;
    BASETYPE thetaFrac;
    BASETYPE targetFitness;
    BASETYPE accuracy;
    IFitness*  fClass;
    IRepair*   repClass;
    IScreen* screenClass;
    std::string outFile;

    std::map< std::string, std::string > m_custom;

    SParams( MPI_Comm comm = MPI_COMM_NULL, const char* file = 0, IFitness* fC = 0, IRepair* rC = 0, IScreen* sC = 0 );

    void initWithFile( MPI_Comm comm, const char* file, IFitness* fC, IRepair* rC = 0, IScreen* sC = 0 );
    const char* getCustomParameter( const char* name ) const;
};

//-----------------------------------------------------------
}

//------------------------------------------------------------
#endif
