#ifndef _memadd_
#define _memadd_
inline void memadd(double *dest,double *src,int size,double fac=1.0)
{
    for(int i=0;i<size;i++)dest[i]+=src[i]*fac;
}
#endif

#ifndef IO_GENERAL_CL_H
#define IO_GENERAL_CL_H

#define _CPP_
#include	"io_general.h"
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>

using namespace std;

#define INF -65535


struct int_rand
{
    int virtual i_rand(int max) const
    {
      return rand()%max;
    }
};


typedef struct
{
   int* info_shape;
   int* info;
   double* Re;
   double* Im;
   
}iog_type;

class datatype
{
public:
    filetype type;
    
    int &ndim;
    char name[500];
     int size;
    one_dim *dim;
    
    int name_c2i(char *name);
    void name_i2c(int ind,char *name);
    
    datatype():ndim(type.head.n_dimensions)
//    {type.head.n_dimensions=0;ndim=type.head.n_dimensions;dim=type.head.dimensions;}
    {type.head.n_dimensions=0;dim=type.head.dimensions;}
    
    void add_dimension(int name,int size);

    void add_dimension(int name,int size,int *list);

    void add_dimension(int name,int size,double *list);
    
    void add_dimension(one_dim &src);

    void insert_dimension(int name,int id);
    
    void clear_ind()
    {ndim=0;}
    
     void print();
     
     void load_type();
     
     void save_type();
     
     one_dim &fdim(int name)
     {
        return dim[find_dim(name)];
     }
     int nind(int name)
     {
        int ind=find_dim(name);
        if(ind>=0)return dim[ind].n_indices;
        else return 1;
     }
     int find_index(int ind,int name)
     {
        one_dim& dimX=fdim(name);
        for(int i=0;i<dimX.n_indices;i++)
           if(dimX.indices[i]==ind) return i;
        return -1;
     }
     
     void print_void_ind(int name);
     
     int find_dim(int name);
     
     int remove_ind(int name);
     
     int sort_ind(int ind,int *map,bool sort);
     
//     int combine_all_head(*datatype res,int size0);
     
    int common_ind(datatype &res,int name,one_dim &dest);
};

class general_data_base:public datatype
{

     general_data_base* set_pointer(general_data_base &res,general_data_base &tmp);

public:
     
     void set_size(int ind_f,int &size1,int &size2);

//     filetype type;
     double * data;

////////////////////
// initialize
     general_data_base()
     {
         ndim=0;
         data=NULL;
     }
     general_data_base(const char*namex)
     {
         sprintf(name,"%s",namex);
         ndim=0;
         data=NULL;
     }

     void initialize(int flag=0);

     general_data_base(filetype &type0,const char * namex)
     {    sprintf(name,"%s",namex);data=NULL;type=type0;initialize(); }
    
     general_data_base(filetype &type0,const char * namex,double* data0)
     {    sprintf(name,"%s",namex);data=NULL;type=type0;initialize();memcpy(data,data0,size*sizeof(double)); }
    
     ~general_data_base(){if(data!=NULL) {delete [] data;data=NULL;}}
// initialize
////////////////////

///////////////
//io
     void clear(){if(data!=NULL&&size>0)memset(data,0,sizeof(double)*size);}
     
     void free_data(){if(data!=NULL){delete [] data;data=NULL;}}
     
     void load();
     
     void load(int name,std::vector<int> &list);
     
     void save();
     
     void copy(general_data_base &dest);
     
    
     
//io
///////////////

     int get_size();
     
     iog_type* get_dat();

     void print_all(bool slim=false);

     void print_all(double cut);

     void create_hdf_datset(string * intrprtr);


     double * seek(int n,...);
     
     int new_data(general_data_base &dest,int name,int new_size);
     
     int normal();
  
     int jackknife(general_data_base &res);

     int antijack(general_data_base &res);
     
     int sum(general_data_base &res,int name);
     
     int aver(general_data_base &res,int name,int flag=1);
     
     int aver_cfg(general_data_base &res,int aver_flag=1);

     int make_cfg_bin(general_data_base &res,int bin_size);

     int aver_cfg_bin(general_data_base &res,int bin_size);
     
     int nind_cfg();
     
     int pick(general_data_base &res,int name,int ist,int ied=INF);

     int pick(general_data_base &res,int name,int size,int* list,int flag=1);
     
     int pick(general_data_base &res,int name,std::vector<int> list,int flag=0)
     {
        return pick(res,name,list.size(),list.data(),flag);
     }

     int pick(general_data_base &res,int name,std::vector<double> list)
     {
        int *ilist=new int[list.size()];
        for(int i=0;i<list.size();i++) ilist[i]=list[i]*1000000;
        return pick(res,name,list.size(),ilist);
     }
     int pick(general_data_base &res,one_dim &dimX,int type=-1)
     {
        int r_type=(type==-1)?dimX.type:type;
        return pick(res,r_type,dimX.n_indices,dimX.indices,0);
     }

     
     
};

#define	N_MAX_STRING	1024

class	data_iterator
{
	private:
		general_data_base&	rawdata;
		int		nrunningdims;
		int		runningdims[N_MAX_DIMENSIONS][2];
		char		namebuf[N_MAX_STRING];

		void		sort_indices(int *pind, int nind);

		void		copy_data();

	public:
		general_data_base*	pdata;

				data_iterator(general_data_base& rawdata, int nfixdims, ...)
					:rawdata(rawdata)
				{
					filetype	ty;
					va_list		params;
					int		i, j, k;

					ty.head.n_dimensions = nfixdims;
					va_start(params, nfixdims);
					for(i=0; i<nfixdims; i++)
					{
						ty.head.dimensions[i].type = va_arg(params, int);
						for(j=0; j<rawdata.ndim; j++)
							if(ty.head.dimensions[i].type == rawdata.dim[j].type)
							{
								ty.head.dimensions[i].n_indices = rawdata.dim[j].n_indices;
								for(k=0; k<ty.head.dimensions[i].n_indices; k++)
									ty.head.dimensions[i].indices[k] = rawdata.dim[j].indices[k];
								sort_indices(ty.head.dimensions[i].indices, ty.head.dimensions[i].n_indices);
								break;
							}
						if(j==rawdata.ndim)
						{
							printf("Error: dimension #%d does not exist.\n", ty.head.dimensions[i].type);
							return;
						}
					}
					va_end(params);

					pdata = new general_data_base(ty, namebuf);

					nrunningdims = 0;
					for(j=0; j<rawdata.ndim; j++)
					{
						for(i=0; i<ty.head.n_dimensions; i++)
							if(ty.head.dimensions[i].type == rawdata.dim[j].type)
								break;
						if(i==ty.head.n_dimensions)
						{
							runningdims[nrunningdims][0] = rawdata.dim[j].type;
							runningdims[nrunningdims][1] = rawdata.dim[j].indices[0];
							nrunningdims ++;
						}
					}
					for(i=0; i<ty.head.n_dimensions; i++)
						runningdims[nrunningdims+i][0] = ty.head.dimensions[i].type;

					copy_data();
				}

				~data_iterator() { }

		data_iterator	operator++();

};



int check_dim(one_dim &a,one_dim &b);

int check_index(datatype &a, datatype &b);


#endif
