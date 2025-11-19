#ifndef _EQUATOR_DATA_LOADER_H_
#define _EQUATOR_DATA_LOADER_H_

#include "../../../types/data_loader/data_loader.h"
#include "../model/EQUATOR_C/EQUATOR_C.h"

// data loader for equator model and data
class EQUATOR_data_loader : public Data_Loader
{
public:

    virtual void load(const std::vector<double> & data) override;

    EQUATOR_data_loader(EQUATOR_C *Eq, int hdp, int vdp, int ap, int frp, int fip, int cp):
    Data_Loader(Eq->forward_size), E(Eq), hor_dist_pos(hdp), ver_dist_pos(vdp), alt_pos(ap),
    freq_re_pos(frp), freq_im_pos(fip), chan_pos(cp){};

private:

    EQUATOR_C *E;

    int hor_dist_pos; // position of hor dist in data vector
    int ver_dist_pos; // position of ver dist in data vector
    int alt_pos;      // position of altitude in data vector

    int freq_re_pos;  // position of re fd measurements in data vector
    int freq_im_pos;  // position of im fd measurements in data vector
    int chan_pos;     // position of td measurements in data vector
};
#endif