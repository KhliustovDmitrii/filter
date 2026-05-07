#include <vector>
#include <random>
#include <iostream>
#include <fstream>

#include "model/polynomial_model.h"

#include "core/filter/kalman_extended/kalman_extended.h"
#include "core/filter/kalman_unscented/kalman_unscented.h"
#include "core/updater/decay_updater/decay_updater.h"
#include "types/workspace/filter_workspace.h"

std::ostream& print_vector(std::ostream &where, std::vector<double> &what);
std::ostream& print_matrix_diagonal(std::ostream &where, std::vector<double> &what, size_t msize);
std::ostream& print_model(std::ostream &where, filter::Model& m);

int main()
{
    int i, iter;
    double residual, residual_min;
    double factor;

    std::ofstream mes("polynomial_measurements.XYZ"); // measurements and model response
    std::ofstream coef("polynomial_coefs.XYZ"); // coefficients and variances

    // ------------ MODEL

    std::vector<double> points{0, -0.5, 0.7, 1};

    // target model
    std::vector<double> coeffs{3, 2, 1};
    filter::examples::Polynomial_Model target_model(points, coeffs, 2);

    // model to fit
    std::vector<double> start_coeffs(coeffs.size(), 0);
    filter::examples::Polynomial_Model model(points, start_coeffs, 2);

    // ------------ FILTER

    filter::Kalman_Extended filter(model);

    auto R = filter.get_R();
    auto S = filter.get_S();

    double err_par = 1.;
    for(i=0; i<model.num_pars; i++)
        S[i*model.num_pars + i] = err_par;

    double err_data = 0.5;
    for(i=0; i<model.forward_size; i++)
        R[i*model.forward_size + i] = err_data;

    filter.set_R(R);
    filter.set_S(S);

    auto filter_ws = filter.allocate_workspace();

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
    print_model(std::cout, target_model) << std::endl;
    // print_model(coef, target_model) << std::endl;

    std::cout << "MEASUREMENTS:" << std::endl;
    print_vector(std::cout, measurements) << std::endl;
    print_vector(mes, measurements) << 0 << std::endl;

    std::cout << "MODEL PARAMETERS:" << std::endl;
    print_model(std::cout, model) << std::endl;
    
    print_model(coef, model);
    print_matrix_diagonal(coef, S, model.num_pars) << std::endl;

    model.response(response);
    residual_min = model.residual(measurements, response);
    std::cout << "START RESIDUAL:   " << residual_min << std::endl;

    print_vector(mes, response) << residual_min << std::endl;

    for(iter=0; iter<30; iter++)
    {
        filter.get_update(measurements, upd_vec, upd_cov, *filter_ws);
        factor = updater.update(upd_vec, measurements);
        filter.update_covariance(*filter_ws, upd_cov, factor);


        model.response(response);
        std::cout << "response " << iter << ": ";
        print_vector(std::cout, response) << std::endl;

        residual = model.residual(measurements, response);
        std::cout << "residual " << iter << ": " << residual << std::endl;

        auto P_posterior = filter.get_P();

        print_model(coef, model);
        print_matrix_diagonal(coef, P_posterior, model.num_pars) << std::endl;

        print_vector(mes, response) << residual << std::endl;

        if(residual > residual_min)
        {
            std::cout << "RESIDUAL INCREASED! Terminating..." << std::endl;
            break;
        }

        residual_min = residual;

        if(residual < err_data*0.05) break;
    }

    std::cout << "MODEL PARAMETERS:" << std::endl;
    print_model(std::cout, model) << std::endl;

    std::cout << "RESIDUAL:   " << residual_min << std::endl;

    return 0;
}

std::ostream& print_vector(std::ostream &where, std::vector<double> &what)
{
    for(size_t i=0; i<what.size(); i++)
        where << what[i] << " ";

    return where;
}

std::ostream& print_matrix_diagonal(std::ostream &where, std::vector<double> &what, size_t msize)
{
    for(size_t i=0; i<msize; i++)
        where << what[i*msize + i] << " ";

    return where;
}

std::ostream& print_model(std::ostream &where, filter::Model& m)
{
    for(size_t i=0; i<m.num_pars; i++)
        where << m.get_param(i) << " ";

    return where;
}
