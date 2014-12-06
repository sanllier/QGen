#ifndef MPICHECK_H
#define MPICHECK_H

//------------------------------------------------------------

#define CHECK( __ERR_CODE__ ) checkMPIRes( __ERR_CODE__, __FUNCTION__ )

void checkMPIRes( int errCode, const char* location );

//------------------------------------------------------------
#endif
