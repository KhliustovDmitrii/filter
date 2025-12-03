#include <vector>
#include <random>
#include <iostream>

#include "model/polynomial_model.h"

#include "core/filter/kalman_extended/kalman_extended.h"
#include "core/filter/kalman_unscented/kalman_unscented.h"
#include "core/updater/decay_updater/decay_updater.h"

int main()
{
    int i, iter;
    double residual, residual_min;

    // ------------ MODEL

    std::vector<double> points{0, 1, 2, 5, 6, 7};

    // target model

    std::vector<double> coeffs{1, 2, 3};
    Polynomial_Model target_model(points, coeffs, 2);

    // model to fit
    std::vector<double> start_coeffs(coeffs.size(), 0);
    Polynomial_Model model(points, start_coeffs, 2);

    // ------------ FILTER

    // extended
    Kalman_Extended filter_ext(&model);

    auto R = filter_ext.get_R();
    auto S = filter_ext.get_S();

    double err_par = 0.1;
    for(i=0; i<model.num_pars; i++)
        S[i*model.num_pars + i] = err_par;

    double err_data = 0.1;
    for(i=0; i<model.forward_size; i++)
        R[i*model.forward_size + i] = err_data;

    filter_ext.set_R(R);
    filter_ext.set_S(S);

    // unscented
    Kalman_Unscented filter_uns(&model);

    filter_uns.set_R(R);
    filter_uns.set_S(S);


    // ------------ UPDATER

    Decay_Updater updater(&model, 3, 0.5);

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

    std::cout << "MEASUREMENTS:" << std::endl;
    for(i=0; i<target_model.forward_size; i++)
        std::cout << measurements[i] << " ";
    std::cout << std::endl;

    std::cout << "MODEL PARAMETERS:" << std::endl;
    for(i=0; i<model.num_pars; i++)
        std::cout << model.get_param(i) << " ";
    std::cout << std::endl;

    model.response(response);
    residual_min = model.residual(measurements, response);
    std::cout << "START RESIDUAL:   " << residual_min << std::endl;

    for(iter=0; iter<10; iter++)
    {
        filter_uns.get_update(measurements, upd_vec, upd_cov);
        updater.update(upd_vec, measurements);

        model.response(response);
        std::cout << "+++++   ";
        for(i=0; i<model.forward_size; i++)
            std::cout << response[i] << " ";
        std::cout << std::endl;

        residual = model.residual(measurements, response);
        std::cout << "------   " << residual << std::endl;

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