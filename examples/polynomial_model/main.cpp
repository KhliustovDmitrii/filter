#include <vector>
#include <random>
#include <iostream>
#include <fstream>

#include "model/polynomial_model.h"

#include "core/filter/kalman_extended/kalman_extended.h"
#include "core/filter/kalman_unscented/kalman_unscented.h"
#include "core/updater/decay_updater/decay_updater.h"
#include "types/workspace/filter_workspace.h"

int main()
{
    int i, iter;
    double residual, residual_min;
    double factor;

    std::ofstream mes("polynomial_measurements.XYZ"); // measurements and model response
    std::ofstream coef("polynomial_coefs.XYZ"); // coefficients and variances

    // ------------ MODEL

    std::vector<double> points{0, -5, -10, 5, 10};

    // target model
    std::vector<double> coeffs{3, 2, 1, 2, 3};
    filter::examples::Polynomial_Model target_model(points, coeffs, 2);

    // model to fit
    std::vector<double> start_coeffs(coeffs.size(), 0);
    filter::examples::Polynomial_Model model(points, start_coeffs, 2);

    // ------------ FILTER

    // extended
    filter::Kalman_Unscented filter_ext(model);

    auto R = filter_ext.get_R();
    auto S = filter_ext.get_S();

    double err_par = 0.3;
    for(i=0; i<model.num_pars; i++)
        S[i*model.num_pars + i] = err_par;

    double err_data = 1.;
    for(i=0; i<model.forward_size; i++)
        R[i*model.forward_size + i] = err_data;

    filter_ext.set_R(R);
    filter_ext.set_S(S);

    auto extended_ws = filter_ext.allocate_workspace();

    // ------------ UPDATER

    filter::Decay_Updater updater(model, 3, 0.5);

    std::vector<double> measurements(target_model.forward_size, 0);
    std::vector<double> response(model.forward_size, 0);
    std::vector<double> upd_vec(model.num_pars, 0);
    std::vector<double> upd_cov(model.num_pars*model.num_pars, 0);

    // make fake data
    target_model.response(measurements);

    // add noise
    std::mt19937 generator(42);
    std::normal_distribution<double> distribution(0, err_data);

    for(i=0; i<target_model.forward_size; i++)
        measurements[i] = measurements[i] + distribution(generator);

    std::cout << "TARGET MODEL PARAMETERS:" << std::endl;
    for(i=0; i<target_model.num_pars; i++)
        std::cout << target_model.get_param(i) << " ";
    std::cout << std::endl;

    for(i=0; i<target_model.num_pars; i++)
        coef << target_model.get_param(i) << " ";
    coef << std::endl;

    std::cout << "MEASUREMENTS:" << std::endl;
    for(i=0; i<target_model.forward_size; i++)
        std::cout << measurements[i] << " ";
    std::cout << std::endl;

    for(i=0; i<target_model.forward_size; i++)
        mes << measurements[i] << " ";
    mes << 0 << std::endl;

    std::cout << "MODEL PARAMETERS:" << std::endl;
    for(i=0; i<model.num_pars; i++)
        std::cout << model.get_param(i) << " ";
    std::cout << std::endl;

    for(i=0; i<model.num_pars; i++)
        coef << model.get_param(i) << " ";

    for(i=0; i<model.num_pars; i++)
        coef << S[i*model.num_pars + i] << " ";    
    coef << std::endl;

    

    model.response(response);
    residual_min = model.residual(measurements, response);
    std::cout << "START RESIDUAL:   " << residual_min << std::endl;

    for(i=0; i<model.forward_size; i++)
        mes << response[i] << " ";
    mes << residual_min;
    mes << std::endl;

    for(iter=0; iter<20; iter++)
    {
        filter_ext.get_update(measurements, upd_vec, upd_cov, *extended_ws);
        factor = updater.update(upd_vec, measurements);
        filter_ext.update_covariance(*extended_ws, upd_cov, factor);


        model.response(response);
        std::cout << "+++++   ";
        for(i=0; i<model.forward_size; i++)
            std::cout << response[i] << " ";
        std::cout << std::endl;

        residual = model.residual(measurements, response);
        std::cout << "------   " << residual << std::endl;

        auto S_modified = filter_ext.get_S();
        filter_ext.set_S(S);

        for(i=0; i<model.num_pars; i++)
            coef << model.get_param(i) << " ";

        for(i=0; i<model.num_pars; i++)
            coef << S_modified[i*model.num_pars + i] << " ";    
        coef << std::endl;

        for(i=0; i<model.forward_size; i++)
            mes << response[i] << " ";

        mes << residual;
        mes << std::endl;

        if(residual > residual_min)
        {
            std::cout << "RESIDUAL INCREASED!!!" << std::endl;
            break;
        }

        residual_min = residual;

        if(residual < err_data) break;
    }

    std::cout << "MODEL PARAMETERS:" << std::endl;
    for(i=0; i<model.num_pars; i++)
        std::cout << model.get_param(i) << " ";
    std::cout << std::endl;

    std::cout << "RESIDUAL:   " << residual_min << std::endl;

    return 0;
}