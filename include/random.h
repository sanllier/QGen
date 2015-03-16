#ifndef RANDOM_H
#define RANDOM_H

#include "defs.h"

//------------------------------------------------------------

class IRandom
{
public:
    virtual ~IRandom() {};

    virtual void setSeed( unsigned seed ) = 0;
    virtual BASETYPE next() = 0;
};

//------------------------------------------------------------
#endif
