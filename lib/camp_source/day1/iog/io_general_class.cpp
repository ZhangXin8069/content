
#include "io_general_class.h"
#include <string>
#include <iostream>
#include <vector>
#include <Python.h>

using namespace std;

int general_data_base::get_size()
{   
   int datsize = size / dim[ndim - 1].n_indices;
   return datsize;
}


iog_type* general_data_base::get_dat()
{   
   int* info = (int* )malloc(sizeof(int)*(size / dim[ndim - 1].n_indices * ndim));
   int* info_shape = (int* )malloc(sizeof(int)*2);
   double* Re = (double* )malloc(sizeof(double)*(size / dim[ndim - 1].n_indices));
   double* Im = (double* )malloc(sizeof(double)*(size / dim[ndim - 1].n_indices));

   iog_type* ptr = (iog_type *)malloc(sizeof(iog_type));

   info_shape[0] = size / dim[ndim - 1].n_indices;
   info_shape[1] = ndim-1;

   for (int i = 0; i < size; i += dim[ndim - 1].n_indices)
   {  
      int size0 = size, i0 = i;
      int flag = 0;
      for (int j = 0; j < dim[ndim - 1].n_indices; j++)
         if (fabs(data[i + j]) > -1)
            flag = 1;
      if (flag == 0)
         continue;
      for (int j = 0; j < ndim - 1; j++)
      {
         size0 /= dim[j].n_indices;
   //    printf("%10d->%10d", j, dim[j].indices[i0 / size0]);
         info[i / dim[ndim - 1].n_indices * (ndim-1) + j] = dim[j].indices[i0 / size0] ;
         i0 = i0 % size0;
      }
      
      Re[i/ dim[ndim - 1].n_indices] = data[i + 0];
      Im[i/ dim[ndim - 1].n_indices] = data[i + 1];

      /*
      for (int j = 0; j < dim[ndim - 1].n_indices; j++)
      {    
          dat[i / dim[ndim - 1].n_indices][j] = data[i + j];
      }*/

   }
   ptr->info = info;
   ptr->info_shape = info_shape;
   ptr->Re = Re;
   ptr->Im = Im;

   return ptr;
}

void general_data_base::print_all(bool slim)
{
   double cut = (slim == true) ? 1e-15 : -1.0;
   print_all(cut);
}

void general_data_base::print_all(double cut)
{
   for (int i = 0; i < size; i += dim[ndim - 1].n_indices)
   {
      int size0 = size, i0 = i;
      int flag = 0;
      for (int j = 0; j < dim[ndim - 1].n_indices; j++)
         if (fabs(data[i + j]) > cut)
            flag = 1;
      if (flag == 0)
         continue;
      for (int j = 0; j < ndim - 1; j++)
      {
         size0 /= dim[j].n_indices;
         printf("%10d->%10d", j, dim[j].indices[i0 / size0]);
         i0 = i0 % size0;
      }
      for (int j = 0; j < dim[ndim - 1].n_indices; j++)
         //if(fabs(data[i+j])<1e4&&fabs(data[i+j])>1e-3)
         //printf("Re->%13.5f",data[i+j]);
         //else
         printf("  cmplx_%d->%13.5e", j, data[i + j]);
      printf("\n");
   }
   fflush(stdout);
}

void general_data_base::create_hdf_datset(string *intrprtr)
{
   string pth[size / dim[ndim - 1].n_indices];
   string cmplx[2] = {"Re", "Im"};
   for (int i = 0; i < size; i += dim[ndim - 1].n_indices)
   {

      int size0 = size, i0 = i;
      int flag = 0;
      for (int j = 0; j < dim[ndim - 1].n_indices; j++)
         if (fabs(data[i + j]) > -1)
            flag = 1;
         
      if (flag == 0)
         continue;
      pth[i / dim[ndim - 1].n_indices] = "/";
      for (int j = 0; j < ndim - 1; j++)
      {
         size0 /= dim[j].n_indices;
   //    printf("%10d->%10d", j, dim[j].indices[i0 / size0]);
         pth[i / dim[ndim - 1].n_indices] += intrprtr[j] + '_' + to_string(dim[j].indices[i0 / size0]) + '/';
         i0 = i0 % size0;
      }
      for (int j = 0; j < dim[ndim - 1].n_indices; j++)
      {
 //        pth[i / dim[ndim - 1].n_indices] += cmplx[j];
 //        printf("  cmplx_%d->%13.5e", j, data[i + j]);
      }
   // printf("\n");
   }

   for (int i = 0; i < size / dim[ndim - 1].n_indices; i++)
   {
     cout << pth[i] << endl;
   }
}


void general_data_base::initialize(int flag)
{
   if (data != NULL)
   {
      delete[] data;
      data = NULL;
   }
   if (ndim > 0)
   {
      size = 1;
      for (int i = 0; i < ndim; i++)
         size *= dim[i].n_indices;
      if (flag == 0)
         if (size > 0)
         {
            data = new double[size];
            memset(data, 0, size * sizeof(double));
         }
         else
         {
            printf("the file size is zero\n");
         }
   }
}

void datatype::load_type()
{
   FILE *fp;

   if ((fp = fopen(name, "rb")) == NULL)
   {
      printf("Failed opening file %s for read.\n", name);
   }

   if (fread(&type, sizeof(filetype), 1, fp) != 1)
   {
      printf("Failed reading the head of file %s.\n", name);
   }

   if (fp != NULL)
      fclose(fp);
}

void datatype::save_type()
{
   FILE *fp;

   if ((fp = fopen(name, "wb")) == NULL)
   {
      printf("Failed opening file %s for read.\n", name);
   }

   if (fwrite(&type, sizeof(filetype), 1, fp) != 1)
   {
      printf("Failed reading the head of file %s.\n", name);
   }

   if (fp != NULL)
      fclose(fp);
}

void general_data_base::load()
{
   if (data != NULL)
      delete[] data;
   data = xqcd_file_read_once(name, &type);
   if (ndim > 0)
   {
      size = 1;
      for (int i = 0; i < ndim; i++)
         size *= dim[i].n_indices;
   }
}

void general_data_base::load(int name, std::vector<int> &list)
{
   if (data != NULL)
      delete[] data;
   load_type();
}

