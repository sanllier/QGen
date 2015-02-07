#ifdef GPU

#include "qindivid_gpu.h"
#include "qgen.h"
#include "qrotoperator.h"
#include "defs_gpu.h"

#include "cuda_runtime.h"
#include "cuda_error_handler.h"

#include "sharedmtrand.h"
#include "mpicheck.h"

#include <time.h>

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

int QGPUIndivid::deviceId = -1;

//------------------------------------------------------------

extern "C" void launchSetInitialKernel( GPUQbit* data, const BASETYPE* rands, long long size );

//------------------------------------------------------------

QGPUIndivid::QGPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] )
    : QBaseIndivid( size, generalComm,  rowComm, coords )
    , m_rand( unsigned( time(0) ) ^ unsigned( coords[0] ))
{
    if ( deviceId == -1 )
    {
        deviceId = selectDevice();
        SAFE_CALL( cudaSetDevice( deviceId ) );
    }

    resize( size );
    setInitial();    
}

//------------------------------------------------------------

QGPUIndivid::~QGPUIndivid()
{
    if ( m_data )
        SAFE_CALL( cudaFree( m_data ) );
}

//------------------------------------------------------------

void QGPUIndivid::resize( long long newSize )
{
    QBaseIndivid::resize( newSize );

    if ( m_data )
        SAFE_CALL( cudaFree( m_data ) );

    SAFE_CALL( cudaMalloc( &m_data, size_t( m_localLogicSize * sizeof( QBit ) ) ) );

    m_rand.resize( newSize );
}

//------------------------------------------------------------

void QGPUIndivid::setInitial()
{
    m_rand.generate();
    launchSetInitialKernel( ( GPUQbit* )m_data, m_rand.getGPUData(), m_localLogicSize );
    m_needRecalcFitness = true;
}

//------------------------------------------------------------

void QGPUIndivid::evolve( const QBaseIndivid& bestInd )
{

}

//------------------------------------------------------------

void QGPUIndivid::bcast( int root )
{

}

//-----------------------------------------------------------

QBaseIndivid& QGPUIndivid::operator=( const QBaseIndivid& rInd )
{
    // CRAP
    const QGPUIndivid& castedIndivid = ( const QGPUIndivid& )rInd;
    if ( m_localLogicSize != castedIndivid.m_localLogicSize )
        throw std::string( "Trying to assing individs with different topologies. " ).append( __FUNCTION__ );

    switch ( rInd.getType() )
    {
        case INDIVID_TYPE_CPU:
        {
            SAFE_CALL( cudaMemcpy( m_data, castedIndivid.m_data, 
                size_t( m_localLogicSize * sizeof( m_data[0] ) ), cudaMemcpyHostToDevice ) );
            break;
        }

        case INDIVID_TYPE_GPU:
        {
            SAFE_CALL( cudaMemcpy( m_data, castedIndivid.m_data, 
                size_t( m_localLogicSize * sizeof( m_data[0] ) ), cudaMemcpyDeviceToDevice ) );
            break;
        }

        default:
        {
            throw std::string( "QCPUIndivid is trying to assign individ with unknown type. " ).append( __FUNCTION__ );  
            break;
        }
    }

    m_observeState      = castedIndivid.m_observeState;
    m_fitness           = castedIndivid.m_fitness;
    m_needRecalcFitness = castedIndivid.m_needRecalcFitness;

    return *this;
}

//-----------------------------------------------------------

BASETYPE QGPUIndivid::getThetaForQBit( const QBaseIndivid& bestInd, long long qbitIndex ) const
{
    return 0.0;
}

//------------------------------------------------------------

int QGPUIndivid::selectDevice()
{
   int deviceCount = 0;
   int suitableDevice = -1;

   SAFE_CALL( cudaGetDeviceCount( &deviceCount ) );
   if ( deviceCount <= 0 )
       throw std::string( "Not found CUDA devices. " ).append( __FUNCTION__ );

   for ( int device = 0; device < deviceCount; ++device )
   {
       cudaDeviceProp devProp;
       SAFE_CALL( cudaGetDeviceProperties ( &devProp, device ) );

       if( devProp.major >= 2 )
          suitableDevice = device ;
   }

   if ( suitableDevice <= -1 )
       throw std::string( "Not found suitable CUDA device. " ).append( __FUNCTION__ );

   return suitableDevice;
}

//------------------------------------------------------------
}

#endif
