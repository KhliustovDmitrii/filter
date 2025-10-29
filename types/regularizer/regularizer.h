#ifndef _REGULARIZER_H_
#define _REGULARIZER_H_

#include "../model/model.h"

class regularizer
{
public:
    model *m;

    virtual void prune() = 0;
}

#endif