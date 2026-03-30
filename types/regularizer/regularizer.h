#ifndef REGULARIZER_H
#define REGULARIZER_H

#include "../model/model.h"

namespace filter::components
{
// prunes model parameters unconditionally
class Regularizer
{

protected:
    Model &m;

public:

    virtual void prune() = 0;

    Regularizer(Model &model) : m(model){};
    virtual ~Regularizer() = default;
};
}; // filter::components
#endif