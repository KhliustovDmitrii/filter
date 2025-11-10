#ifndef _REGULARIZER_H_
#define _REGULARIZER_H_

#include "../model/model.h"

// prunes model parameters unconditionally
class Regularizer
{

protected:
    Model *m;

public:

    virtual void prune() = 0;

    Regularizer(Model *model) : m(model){};
    virtual ~Regularizer() = default;
};

#endif