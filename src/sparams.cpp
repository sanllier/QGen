#include "sparams.h"
#include "mpicheck.h"
#include "pugixml.hpp"

#include <vector>

//------------------------------------------------------------

namespace QGen {
//-----------------------------------------------------------

SParams::SParams( MPI_Comm comm/* = MPI_COMM_NULL*/, const char* file/* = 0*/, IFitness* fC/* = 0*/, IRepair* rC/* = 0*/, IProcessScreen* sC/* = 0*/ )
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
    #ifdef GPU
        , gpu(false)
    #endif
{
    if ( comm == MPI_COMM_NULL )
        throw std::string( "Invalid communicator. " ).append( __FUNCTION__ );

    if ( file )
        initWithFile( comm, file, fC, rC, sC );
}

//-----------------------------------------------------------

void SParams::initWithFile( MPI_Comm comm, const char* file, IFitness* fC, IRepair* rC/* = 0*/, IProcessScreen* sC/* = 0*/ )
{
    if ( !file || !file[0] )
        throw std::string( "Some problems with params file. " ).append( __FUNCTION__ );

    MPI_File fp = MPI_FILE_NULL;
    CHECK( MPI_File_open( comm, file, MPI_INFO_NULL, MPI_MODE_RDONLY, &fp ) );
    
    MPI_Offset fileSize = 0;
    CHECK( MPI_File_get_size( fp, &fileSize ) );
    
    std::vector<char> buf( (unsigned)fileSize );
    MPI_Status status;    
    MPI_File_read_at_all( fp, 0, &buf[0], int( fileSize ), MPI_CHAR, &status );

    CHECK( MPI_File_close( &fp ) );

    pugi::xml_document doc;
    doc.load_buffer( &buf[0], unsigned( fileSize ) );
    
    pugi::xml_node qgenNode = doc.child( "qgen" );    
    if ( !qgenNode )
        throw std::string( "Some problems with params file. " ).append( __FUNCTION__ );

    for ( pugi::xml_node node = qgenNode.child( "parameter" ); node; node = node.next_sibling() )
    {
        const char* name = node.attribute( "name" ).as_string();
        if ( !name )
            continue;

        if ( 0 == strcmp( "cycle-threshold", name ) )
        {
            cycThreshold = (long long)node.attribute( "value" ).as_uint(0);
        }
        else if ( 0 == strcmp( "individs-num", name ) )
        {
            individsNum = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "individ-size", name ) )
        {
            indSize = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "topology-rows", name ) )
        {
            topoRows = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "topology-cols", name ) )
        {
            topoCols = (int)node.attribute( "value" ).as_int(0);
        }
        else if ( 0 == strcmp( "target-fitness", name ) )
        {
            targetFitness = (BASETYPE)node.attribute( "value" ).as_double(0.0);
        }
        else if ( 0 == strcmp( "target-accuracy", name ) )
        {
            accuracy = (BASETYPE)node.attribute( "value" ).as_double(0.0);
        }
        else if ( 0 == strcmp( "out-file", name ) )
        {
            outFile = node.attribute( "value" ).as_string("");
        }
    #ifdef GPU
        else if ( 0 == strcmp( "use-gpu", name ) )
        {
            gpu = node.attribute( "value" ).as_bool( false );
        }
    #endif
    }

    fClass      = fC;
    repClass    = rC;
    screenClass = sC;
}

//------------------------------------------------------------
}
