#include "mpicheck.h"

#include <string>

#include "mpi.h"

//------------------------------------------------------------

void checkMPIRes( int errCode, const char* location )
{
    if ( errCode != MPI_SUCCESS )
    {
        char errText[ MPI_MAX_ERROR_STRING ];
        int len;
        MPI_Error_string( errCode, errText, &len );
        throw std::string( errText ).append( location );
    }
}

//------------------------------------------------------------
