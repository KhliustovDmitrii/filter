#ifndef EQUATOR_C_H
#define EQUATOR_C_H

// compressed representation of EQUATOR
// complies to model interface
#include "../EQUATOR/EQUATOR.h"
#include "types/model/model.h"

#include <memory>
#include <functional>

namespace filter::examples
{
class EQUATOR_data_loader;
class EQUATOR_aggregator;

class EQUATOR_C : public Model
{

public:

    friend class EQUATOR_data_loader;
    friend class EQUATOR_aggregator;

    EQUATOR_C(EQUATOR &EQ, std::vector<double> weights_) : 
    Model(0, 2 * EQ.num_freqs + EQ.num_channels), 
    E(EQ), weights(weights_) {};

    // forward function
    virtual void response(std::vector<double> &resp_arr) const override;

    // proximity in response between model and measurements
    virtual double residual(std::vector<double> &mes, std::vector<double> &resp) const override;

    // getters and setters
    virtual void set_param(int ind, double val) override;
    virtual double get_param(int ind) const override;

    // add differentiable param to model
    template <typename forward_scale, typename backward_scale>
    void add_param(std::vector<double> &model_part, size_t ind)
    {
        params.emplace_back(
            std::make_unique<scalable_parameter<forward_scale, backward_scale>>(model_part, ind));

        num_pars++;
    }

private:

    EQUATOR &E;

    // a class encapsulating model interaction logic
    class abstract_scalable_parameter
    {
    public:
      virtual double get() = 0;
      virtual void set(double val) = 0; 
      virtual ~abstract_scalable_parameter() = default;
    };

    // family of concrete implementations
    template <typename forward_scale, typename backward_scale>
    class scalable_parameter : public abstract_scalable_parameter
    {
    private:
      forward_scale fs;
      backward_scale bs;

      std::vector<double>& model_part;
      size_t model_ind;

    public:
      virtual double get() { return bs(model_part[model_ind]); }
      virtual void set(double val) { model_part[model_ind] = fs(val); }

      scalable_parameter(std::vector<double>& model_part_, size_t model_ind_) :
      abstract_scalable_parameter(), model_part(model_part_), model_ind(model_ind_)
      {};
    };

    std::vector<std::unique_ptr<abstract_scalable_parameter>> params;

    // weights to compute residual with
    std::vector<double> weights;

};
}; // filter::examples
#endif