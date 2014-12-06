#ifndef SHAREDMTRAND_H
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
    inline static MTRand_closed& getClosedInstance()
    {
        static MTRand_closed randomizer( ( unsigned long )( std::time(0) ) );
        return randomizer;
    }

    inline static MTRand_int32& get32UnsignedInstance()
    {
        static MTRand_int32 randomizer( ( unsigned long )( std::time(0) ) );
        return randomizer;
    }
};

//------------------------------------------------------------
#endif
