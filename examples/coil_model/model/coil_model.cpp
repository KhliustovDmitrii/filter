#include <cmath>

#include "coil_model.h"

namespace filter::examples
{

void Coil_Model::response(std::vector<double> &resp_arr) const
{
    double r_abs, r_abs_sq;
    double q0, q1, q2, q3;

    // auxillary variables for computations speedup
    std::vector<double> r_1(3, 0);        // pendulum position w r t first coil
    std::vector<double> r_2(3, 0);        // pendulum position w r t second coil
    std::vector<double> O(9, 0);          // orientation matrix for receiver
    std::vector<double> tmp_vec(3, 0);
      
    r_1[0] = params[0] - coil_pos_1[0];
    r_1[1] = params[1] - coil_pos_1[1];
    r_1[2] = params[2] - coil_pos_1[2];
      
    r_2[0] = params[0] - coil_pos_2[0];
    r_2[1] = params[1] - coil_pos_2[1];
    r_2[2] = params[2] - coil_pos_2[2];   
      
    // quaternionic parametrization of rotation matrix
    q1 = params[3];
    q2 = params[4];
    q3 = params[5];
    q0 = sqrt(1 - q1*q1 - q2*q2 - q3*q3);
    if(std::isnan(q0)) q0 = 0;  
      
    O[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    O[1] = 2*(q1*q2 - q0*q3);
    O[2] = 2*(q0*q2 + q1*q3);
      
    O[3] = 2*(q0*q3 + q1*q2);
    O[4] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    O[5] = 2*(q2*q3 - q0*q1);
      
    O[6] = 2*(q1*q3 - q0*q2);
    O[7] = 2*(q0*q1 + q2*q3);
    O[8] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
      
      
    r_abs_sq = r_1[0]*r_1[0] + r_1[1]*r_1[1] + r_1[2]*r_1[2];
    r_abs = sqrt(r_abs_sq);
      
    tmp_vec[0] = (3*r_1[0]*r_1[0]/r_abs_sq - 1)*mom_1[0] + 3*r_1[0]*r_1[1]*mom_1[1]/r_abs_sq + 3*r_1[0]*r_1[2]*mom_1[2]/r_abs_sq;
    tmp_vec[1] = 3*r_1[0]*r_1[1]*mom_1[0]/r_abs_sq + (3*r_1[1]*r_1[1]/r_abs_sq - 1)*mom_1[1] + 3*r_1[1]*r_1[2]*mom_1[2]/r_abs_sq;
    tmp_vec[2] = 3*r_1[0]*r_1[2]*mom_1[0]/r_abs_sq + 3*r_1[1]*r_1[2]*mom_1[1]/r_abs_sq + (3*r_1[2]*r_1[2]/r_abs_sq - 1)*mom_1[2];   
      
    resp_arr[0] = O[0]*tmp_vec[0] + O[1]*tmp_vec[1] + O[2]*tmp_vec[2];
    resp_arr[1] = O[3]*tmp_vec[0] + O[4]*tmp_vec[1] + O[5]*tmp_vec[2];
    resp_arr[2] = O[6]*tmp_vec[0] + O[7]*tmp_vec[1] + O[8]*tmp_vec[2];  
      
    resp_arr[0] = resp_arr[0]/(r_abs*r_abs*r_abs);
    resp_arr[1] = resp_arr[1]/(r_abs*r_abs*r_abs);
    resp_arr[2] = resp_arr[2]/(r_abs*r_abs*r_abs);  
      
    r_abs_sq = r_2[0]*r_2[0] + r_2[1]*r_2[1] + r_2[2]*r_2[2];
    r_abs = sqrt(r_abs_sq);
      
    tmp_vec[0] = (3*r_2[0]*r_2[0]/r_abs_sq - 1)*mom_2[0] + 3*r_2[0]*r_2[1]*mom_2[1]/r_abs_sq + 3*r_2[0]*r_2[2]*mom_2[2]/r_abs_sq;
    tmp_vec[1] = 3*r_2[0]*r_2[1]*mom_2[0]/r_abs_sq + (3*r_2[1]*r_2[1]/r_abs_sq - 1)*mom_2[1] + 3*r_2[1]*r_2[2]*mom_2[2]/r_abs_sq;
    tmp_vec[2] = 3*r_2[0]*r_2[2]*mom_2[0]/r_abs_sq + 3*r_2[1]*r_2[2]*mom_2[1]/r_abs_sq + (3*r_2[2]*r_2[2]/r_abs_sq - 1)*mom_2[2];  
      
    resp_arr[3] = O[0]*tmp_vec[0] + O[1]*tmp_vec[1] + O[2]*tmp_vec[2];
    resp_arr[4] = O[3]*tmp_vec[0] + O[4]*tmp_vec[1] + O[5]*tmp_vec[2];
    resp_arr[5] = O[6]*tmp_vec[0] + O[7]*tmp_vec[1] + O[8]*tmp_vec[2];  
      
    resp_arr[3] = resp_arr[3]/(r_abs*r_abs*r_abs);
    resp_arr[4] = resp_arr[4]/(r_abs*r_abs*r_abs);
    resp_arr[5] = resp_arr[5]/(r_abs*r_abs*r_abs);  
}

double Coil_Model::residual(std::vector<double> &mes, std::vector<double> &resp) const
{
    double res;
    int i;

    res = 0;
    for(i=0; i<forward_size; i++)
        res += (resp[i] - mes[i])*(resp[i] - mes[i])/(weights[i]*weights[i]);
        
    res = res/(forward_size);

    return sqrt(res);
}

void Coil_Model::set_param(int ind, double val)
{
    params[ind] = val;
}

double Coil_Model::get_param(int ind) const
{
    return params[ind];
}

void Coil_Model::quat_normalize()
{
    double q1, q2, q3, q_abs;
   
    q1 = params[3];
    q2 = params[4];
    q3 = params[5];
         
    q_abs = sqrt(q1*q1 + q2*q2 + q3*q3);
         
    if(q_abs > 1)
    {
        params[3] = params[3]/(q_abs+0.001);
        params[4] = params[4]/(q_abs+0.001);
        params[5] = params[5]/(q_abs+0.001);
    }
}
}; // filter::examples