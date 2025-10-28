#include "coil_model.h"

#define RESP_TEST 0

void coil_model::response(double *resp_arr)
{
   double r_abs, r_abs_sq;
   double q0, q1, q2, q3;
      
   r_1[0] = params[0] - coil_pos_1[0];
   r_1[1] = params[1] - coil_pos_1[1];
   r_1[2] = params[2] - coil_pos_1[2];
      
   r_2[0] = params[0] - coil_pos_2[0];
   r_2[1] = params[1] - coil_pos_2[1];
   r_2[2] = params[2] - coil_pos_2[2];
      
   if(RESP_TEST)   
   {
      std::cout << "r vectors\n";
      std::cout << r_1[0] << " " << r_1[1] << " " << r_1[2] << "\n";
      std::cout << r_2[0] << " " << r_2[1] << " " << r_2[2] << "\n";
   }   
      
   // quaternionic parametrization of rotation matrix
   q1 = params[3];
   q2 = params[4];
   q3 = params[5];
   q0 = sqrt(1 - q1*q1 - q2*q2 - q3*q3);
   if(std::isnan(q0)) q0 = 0;
      
   if(RESP_TEST)   
   {   
      std::cout << "sq_sum     " << q1*q1 + q2*q2 + q3*q3 << "\n";
      std::cout << "q0     " << q0 << "\n";
   }   
      
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
      
   if(RESP_TEST)   
   {   
      std::cout << "tmp vector\n";
      std::cout << tmp_vec[0] << " " << tmp_vec[1] << " " << tmp_vec[2] << "\n";
   }   
      
   resp_arr[0] = O[0]*tmp_vec[0] + O[1]*tmp_vec[1] + O[2]*tmp_vec[2];
   resp_arr[1] = O[3]*tmp_vec[0] + O[4]*tmp_vec[1] + O[5]*tmp_vec[2];
   resp_arr[2] = O[6]*tmp_vec[0] + O[7]*tmp_vec[1] + O[8]*tmp_vec[2];
      
   if(RESP_TEST)   
   {   
      std::cout << "resp before\n";
      std::cout << resp_arr[0] << " " << resp_arr[1] << " " << resp_arr[2] << "\n";
   }   
      
   resp_arr[0] = resp_arr[0]/(r_abs*r_abs*r_abs);
   resp_arr[1] = resp_arr[1]/(r_abs*r_abs*r_abs);
   resp_arr[2] = resp_arr[2]/(r_abs*r_abs*r_abs);
      
   if(RESP_TEST)   
   {   
      std::cout << "resp after\n";
      std::cout << resp_arr[0] << " " << resp_arr[1] << " " << resp_arr[2] << "\n";
   }   
      
      
      
      
   r_abs_sq = r_2[0]*r_2[0] + r_2[1]*r_2[1] + r_2[2]*r_2[2];
   r_abs = sqrt(r_abs_sq);
      
   tmp_vec[0] = (3*r_2[0]*r_2[0]/r_abs_sq - 1)*mom_2[0] + 3*r_2[0]*r_2[1]*mom_2[1]/r_abs_sq + 3*r_2[0]*r_2[2]*mom_2[2]/r_abs_sq;
   tmp_vec[1] = 3*r_2[0]*r_2[1]*mom_2[0]/r_abs_sq + (3*r_2[1]*r_2[1]/r_abs_sq - 1)*mom_2[1] + 3*r_2[1]*r_2[2]*mom_2[2]/r_abs_sq;
   tmp_vec[2] = 3*r_2[0]*r_2[2]*mom_2[0]/r_abs_sq + 3*r_2[1]*r_2[2]*mom_2[1]/r_abs_sq + (3*r_2[2]*r_2[2]/r_abs_sq - 1)*mom_2[2];
      
   if(RESP_TEST)   
   {   
      std::cout << "tmp vector\n";
      std::cout << tmp_vec[0] << " " << tmp_vec[1] << " " << tmp_vec[2] << "\n";
   }   
      
   resp_arr[3] = O[0]*tmp_vec[0] + O[1]*tmp_vec[1] + O[2]*tmp_vec[2];
   resp_arr[4] = O[3]*tmp_vec[0] + O[4]*tmp_vec[1] + O[5]*tmp_vec[2];
   resp_arr[5] = O[6]*tmp_vec[0] + O[7]*tmp_vec[1] + O[8]*tmp_vec[2];
      
   if(RESP_TEST)   
   {   
      std::cout << "resp before\n";
      std::cout << resp_arr[3] << " " << resp_arr[4] << " " << resp_arr[5] << "\n";
   }   
      
   resp_arr[3] = resp_arr[3]/(r_abs*r_abs*r_abs);
   resp_arr[4] = resp_arr[4]/(r_abs*r_abs*r_abs);
   resp_arr[5] = resp_arr[5]/(r_abs*r_abs*r_abs);
      
   if(RESP_TEST)   
   {   
      std::cout << "resp after\n";
      std::cout << resp_arr[3] << " " << resp_arr[4] << " " << resp_arr[5] << "\n";
   }   
}

void coil_model::quat_normalize()
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
