#include "text_utils.h"

int parse(std::string buf, double *values)
{
   int is_digit = 0;
   int j = 0;
   int i;
   int read_num = 0;
   char tmp[256];
   int ret_code;
    
   for(i = 0; i < buf.size(); i++)
   {
      if(buf[i]=='*')
      {
         if(read_num==0) ret_code = -1;
         else ret_code = 1;
         continue;
      }
       
      if((buf[i]>='0' && buf[i]<='9')||buf[i]=='.'||buf[i]=='-')
      {
         tmp[j] = buf[i];
         j++;
         is_digit = 1;
      } 
      else
      {
         if(is_digit==1)
         {
            values[read_num] = atof(tmp);
            read_num++;
            is_digit = 0;
            for(j = 0; j < 256; j++)
               tmp[j] = (char) 0;
         }
         j = 0;
      }
   }
   return ret_code;
}
