#ifndef _UPDATER_H_
#define _UPDATER_H_
class updater
{
public:
    model *m;

    virtual double update() = 0;
}
#endif