#ifdef GPU

#include "qindivid_gpu.h"
#include "qgen.h"
#include "qrotoperator.h"

#include "cuda_runtime.h"
#include "cuda_error_handler.h"

#include "sharedmtrand.h"
#include "mpicheck.h"

//------------------------------------------------------------

namespace QGen {
//------------------------------------------------------------

int QGPUIndivid::deviceId = -1;

//------------------------------------------------------------

void launchSetInitialKernel( QBit* data, long long size );

//------------------------------------------------------------

QGPUIndivid::QGPUIndivid( long long size, MPI_Comm generalComm, MPI_Comm rowComm, int coords[2] )
    : QBaseIndivid( size, generalComm,  rowComm, coords )
    , m_deviceData(0)
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
    if ( m_deviceData )
        CHECK( cudaFree( m_deviceData ) );
}

//------------------------------------------------------------

void QGPUIndivid::resize( long long newSize )
{
    QBaseIndivid::resize( newSize );

    if ( m_deviceData )
        CHECK( cudaFree( m_deviceData ) );

    CHECK( cudaMalloc( &m_deviceData, size_t( m_localLogicSize * sizeof( QBit ) ) ) );
}

//------------------------------------------------------------

void QGPUIndivid::setInitial()
{
    launchSetInitialKernel( m_deviceData, m_localLogicSize );
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
