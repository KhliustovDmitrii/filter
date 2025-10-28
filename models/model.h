#ifndef _MODEL_H_
#define _MODEL_H_

template <typename E, typename T>
class model
{
public:
   int num_pars;
   int forward_size;
   E *params;
   model() {};
   virtual void response(T *resp_arr)
   {
      
   };
   virtual double update(T *upd_vec, T *resp, T res, T *R)
   {
      return 1;
   };
   virtual void raw_to_conv()
   {
      
   };
   virtual void conv_to_raw()
   {

   };
};

#endif
