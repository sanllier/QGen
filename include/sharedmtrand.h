#if !defined( SHAREDMTRAND_H ) && !defined( NOTMTRAND )
#define SHAREDMTRAND_H

#include "mtrand.h"

#include <ctime>

//------------------------------------------------------------

class SharedMTRand
{
private:
    SharedMTRand()  {}
    ~SharedMTRand() {}

public:
    inline static MTRand_closed& getClosedInstance( unsigned seed = 0 )
    {
        static MTRand_closed randomizer( unsigned( time(0) ) ^ seed );
        return randomizer;
    }

    inline static MTRand_int32& get32UnsignedInstance( unsigned seed = 0 )
    {
        static MTRand_int32 randomizer( unsigned( time(0) ) ^ seed );
        return randomizer;
    }
};

//------------------------------------------------------------
#endif
