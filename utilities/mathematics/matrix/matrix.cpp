#include <cmath>
#include "matrix.h"

void invert_matrix(std::vector<double> &dest, std::vector<double> &source, int dim)
{
   int i, j, k;
   double tmp;
   
   for(i=0; i<dim; i++)
      for(j=0; j<dim; j++)
         dest[i*dim + j] = (i==j)?1:0;
         
   for(i=0; i<dim; i++) // use all rows
   {
      // current diagonal element is zero
      // change with row with non-zero
      j=i;
      while(source[j*dim+i]==0) j++;
      
      for(k=0; k<dim; k++)
      {
         tmp = source[j*dim + k];
         source[j*dim + k] = source[i*dim + k];
         source[i*dim + k] = tmp;
         
         tmp = dest[j*dim + k];
         dest[j*dim + k] = dest[i*dim + k];
         dest[i*dim + k] = tmp;
      }
      
      
      // nonzero diag - main part
      for(j=0; j<dim; j++) // all rows but current
      {
         if(i==j) continue;
         
         if(source[j*dim+i] == 0) continue;
         
         tmp = - source[j*dim+i]/source[i*dim+i];
         
         for(k=i; k<dim; k++)
            source[j*dim+k] = source[j*dim+k] + tmp*source[i*dim+k];
            
         for(k=0; k<dim; k++)
            dest[j*dim+k] = dest[j*dim+k] + tmp*dest[i*dim+k];
         
      }
      
      // normalize current row
      tmp = 1/source[i*dim+i];
      for(k=i; k<dim; k++)
         source[i*dim+k] = source[i*dim+k]*tmp;
         
      for(k=0; k<dim; k++)
         dest[i*dim+k] = dest[i*dim+k]*tmp;
   }
   
}

void matrix_square_root(std::vector<double> &dest, std::vector<double> &source, int dim)
{
   int i, j, k;
   
   for(i=0; i<dim; i++)
   {
      dest[i*dim + i] = source[i*dim + i];
      for(j=0; j<i; j++)
         dest[i*dim + i] = dest[i*dim + i] - dest[i*dim + j]*dest[i*dim + j];
      dest[i*dim + i] = sqrt(dest[i*dim + i]);
      
      for(j=0; j<i; j++)
         dest[j*dim + i] = 0;
         
      for(j=i+1; j<dim; j++)
      {
         dest[j*dim + i] = source[j*dim + i];
         for(k=0; k<i; k++)
            dest[j*dim + i] = dest[j*dim + i] - dest[j*dim + k]*dest[i*dim + k];
            
         dest[j*dim + i] = dest[j*dim + i]/dest[i*dim + i];
      }
   } 
}