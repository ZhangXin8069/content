
#include 	"io_general_class.h"

double effective_mass(int nt_h,int it,double a,double b,int flag,double x1,double x2)
{
    if(flag==0)return log(a/b);
    int pt1=nt_h-it,pt2=nt_h-it-1;
    double rtn=0.5*(x1+x2);
    double ch1,ch2,sh1,sh2;
    double fun,dfdx,dx;

    for(int iter=0;iter<1000;iter++)
    {
       ch1=cosh(pt1*rtn);sh1=sinh(pt1*rtn);
       ch2=cosh(pt2*rtn);sh2=sinh(pt2*rtn);
       if(flag>0)
       {
          fun=a*ch2-b*ch1;dfdx=a*sh2*pt2-b*sh1*pt1;
       }
       else
       {
          fun=a*sh2-b*sh1;dfdx=a*ch2*pt2-b*ch1*pt1;
       }
       dx=fun/dfdx;
       rtn-=fun/dfdx;
       if(fabs(dx)<1e-7) break;
    }
    return fabs(rtn);
}


void jack(double* res,int nconf,int nd, int nc, double* res_jk)
{
    int i,iconf;
    double *aver=(double*)malloc(sizeof(double)*nc);
    for(i=0;i<nc;i++) {aver[i]=0.0;}
    for(iconf=0;iconf<nconf;iconf++) for(i=0;i<nc;i++) aver[i]+=res[nd*iconf+i];
    for(iconf=0;iconf<nconf;iconf++) for(i=0;i<nc;i++) res_jk[nd*iconf+i]=(aver[i]-res[nd*iconf+i])/(nconf-1);
    if(nconf==1) for(i=0;i<nc;i++) res_jk[i]=aver[i];
    free(aver);
}

void anti_jack(double* res,int nconf,int nd, int nc, double* res_jk)
{
    int i,iconf;
    double *aver=(double*)malloc(sizeof(double)*nc);
    for(i=0;i<nc;i++) {aver[i]=0.0;}
    for(iconf=0;iconf<nconf;iconf++) for(i=0;i<nc;i++) aver[i]+=res[nd*iconf+i];
    for(iconf=0;iconf<nconf;iconf++) for(i=0;i<nc;i++) res_jk[nd*iconf+i]=aver[i]-(nconf-1)*res[nd*iconf+i];
    if(nconf==1) for(i=0;i<nc;i++) res_jk[i]=aver[i];
    free(aver);
}


void jacksig(double* res,int nconf,int nd, int nc, double* aver,double* sig)
{
    int i,iconf;

    for(i=0;i<nc;i++) {aver[i]=0.0;sig[i]=0.0;}
    for(iconf=0;iconf<nconf;iconf++) for(i=0;i<nc;i++) aver[i]+=res[nd*iconf+i];
    for(i=0;i<nc;i++) aver[i]/=nconf;
    for(iconf=0;iconf<nconf;iconf++) for(i=0;i<nc;i++) sig[i]+=pow(aver[i]-res[nd*iconf+i],2);
    for(i=0;i<nc;i++) sig[i]=sqrt(sig[i]*(nconf-1)/(nconf));

}


int datatype::name_c2i(char *name)
{
    if(!strcmp(name,"dim_complex"))return dim_complex;
    if(!strcmp(name,"dim_mass"))return dim_mass;
    if(!strcmp(name,"dim_smear"))return dim_smear;
    if(!strcmp(name,"dim_conf"))return dim_conf;
    if(!strcmp(name,"dim_jackknife"))return dim_jackknife;
    if(!strcmp(name,"dim_operator"))return dim_operator;
    if(!strcmp(name,"dim_t"))return dim_t;
    if(!strcmp(name,"dim_mass2"))return dim_mass2;
    if(!strcmp(name,"dim_t2"))return dim_t2;
    if(!strcmp(name,"dim_channel"))return dim_channel;
    if(!strcmp(name,"dim_momentum"))return dim_momentum;
    if(!strcmp(name,"dim_param"))return dim_momentum;
    if(!strcmp(name,"dim_bootstrap"))return dim_bootstrap;
    return 0;
}

void datatype::name_i2c(int ind,char *name)
{
    if(ind==dim_complex){sprintf(name,"dim_complex");return;}
    if(ind==dim_mass){sprintf(name,"dim_mass");return;}
    if(ind==dim_smear){sprintf(name,"dim_smear");return;}
    if(ind==dim_conf){sprintf(name,"dim_conf");return;}
    if(ind==dim_jackknife){sprintf(name,"dim_jackknife");return;}
    if(ind==dim_operator){sprintf(name,"dim_operator");return;}
    if(ind==dim_t){sprintf(name,"dim_t");return;}
    if(ind==dim_mass2){sprintf(name,"dim_mass2");return;}
    if(ind==dim_t2){sprintf(name,"dim_t2");return;}
    if(ind==dim_channel){sprintf(name,"dim_channel");return;}
    if(ind==dim_momentum){sprintf(name,"dim_momentum");return;}
    if(ind==dim_param){sprintf(name,"dim_param");return;}
    if(ind==dim_temporary){sprintf(name,"dim_temporary");return;}
    if(ind==dim_bootstrap){sprintf(name,"dim_bootstrap");return;}
    sprintf(name,"none");
}


    void datatype::add_dimension(int name,int size)
    {
       int icount=ndim;
        dim[icount].type=name;
        dim[icount].n_indices=size;
        for(int i=0;i<size;i++) dim[icount].indices[i]=i;
        ndim++;
    }

    void datatype::add_dimension(int name,int size,int *list)
    {
       int icount=ndim;
       add_dimension(name,size);
        for(int i=0;i<size;i++) dim[icount].indices[i]=list[i];
    }

    void datatype::add_dimension(int name,int size,double *list)
    {
       int icount=ndim;
       add_dimension(name,size);
        for(int i=0;i<size;i++) dim[icount].indices[i]=(fabs(list[i])<1)?list[i]*1000000:list[i];
    }
    
    void datatype::add_dimension(one_dim &src)
    {
        int icount=ndim;
        dim[icount]=src;
        ndim++;
    }
    
    void datatype::insert_dimension(int name,int id)
    {
        for(int i=ndim-1;i>=0;i--)
           dim[i+1]=dim[i];
        dim[0].type=name;
        dim[0].n_indices=1;
        dim[0].indices[0]=id;
        ndim++;
    }
    
    void datatype::print()
    {
        char name[100];
        for(int i=0;i<ndim;i++)
        {
          
           name_i2c(dim[i].type,name);
           printf("%-20s:(%8d)",name,dim[i].n_indices);
           for(int j=0;j<dim[i].n_indices;j++)
           if(dim[i].indices[i]>10000000)
              printf("%12d",dim[i].indices[j]);
           else
              printf("%8d",dim[i].indices[j]);
           printf("\n");        
        }
        printf("\n");
    }
    
     int datatype::common_ind(datatype &res,int name,one_dim &dest)
     {
         dest.type=name,dest.n_indices=0;
         char dim_name[200];
         int ind1=find_dim(name);
            if(ind1<0){name_i2c(name,dim_name);
             printf("dimension %10s is not exit in data1\n",dim_name);return -1;}
         int ind2=res.find_dim(name);
            if(ind1<0){name_i2c(name,dim_name);
              printf("dimension %10s is not exit in data2\n",dim_name);return -2;}
         for(int i=0;i<dim[ind1].n_indices;i++)
         for(int j=0;j<res.dim[ind2].n_indices;j++)
         if(res.dim[ind2].indices[j]==dim[ind1].indices[i])
            {dest.indices[dest.n_indices]=dim[ind1].indices[i];dest.n_indices++;break;}
         if(dest.n_indices==0){printf("no common index!\n");return -3;}
         return 0;
     }
     
    void general_data_base::print_all(bool slim)
    {
       double cut=(slim==true)?1e-15:-1.0;
       print_all(cut);
    }


    void general_data_base::print_all(double cut)
    {
        print();
        for(int i=0;i<size;i+=dim[ndim-1].n_indices)
        {
          int size0=size,i0=i;
          int flag=0;
          for(int j=0;j<dim[ndim-1].n_indices;j++)
          if(fabs(data[i+j])>cut) flag=1;
          if(flag==0)continue;
          for(int j=0;j<ndim-1;j++)
          {  
             size0/=dim[j].n_indices;
             printf("%10d",dim[j].indices[i0/size0]);
             i0=i0%size0;
          }
          for(int j=0;j<dim[ndim-1].n_indices;j++)
          if(fabs(data[i+j])<1e4&&fabs(data[i+j])>1e-3)
            printf("%13.5f",data[i+j]);
          else
            printf("%13.5e",data[i+j]);
          printf("\n");
        }
        fflush(stdout);
    }

     void general_data_base::initialize(int flag)
     {
        if(data!=NULL){delete [] data;data=NULL;}
        if(ndim>0)
        {  size=1;
           for(int i=0;i<ndim;i++)
             size*=dim[i].n_indices;
          if(flag==0)
          if(size>0)
          { data=new double[size];memset(data,0,size*sizeof(double));}
          else
          {
              printf("the file size is zero\n");
              print();
          }
        }
          
     }
     
     void datatype::load_type()
     {
       FILE *fp;
       
       if((fp=fopen(name, "rb"))==NULL)
        {
                printf("Failed opening file %s for read.\n", name);
        }

        if(fread(&type, sizeof(filetype), 1, fp)!=1)
        {
                printf("Failed reading the head of file %s.\n", name);
        }
        
        if(fp!=NULL) fclose(fp);
     }
     
     void datatype::save_type()
     {
       FILE *fp;
       
       if((fp=fopen(name, "wb"))==NULL)
        {
                printf("Failed opening file %s for read.\n", name);
        }

        if(fwrite(&type, sizeof(filetype), 1, fp)!=1)
        {
                printf("Failed reading the head of file %s.\n", name);
        }
        
        if(fp!=NULL) fclose(fp);
     }
     
/*     int datatype::combine_all_head(*datatype res,int size0)
     {
          if(size0<=0) {printf("size (%6d)<=0!\n",size);return -1;}
          type=res[0].type;
          for(int i=0;i<ndim;i++)
          if(dim[i].n_indices>1)dim[i].n_indices*=-1;//only check the dimension with nind=1
          int flag=0;
          for(int j=1;j<rank;j++)
          { 
             for(int i=0;i<ndim;i++)
             if(dim[i].n_indices<0)
             {
                if(rank[j].dim[i].n_indices+dim[i].n_indices!=0)
                {flag=-1:break;}
             }
             else
             {
               int ioff=-1;
               for(int k=0;k<dim[i].n_indices;k++)
               if(rank[j].dim[i].indices[0]==dim[i].indices[k]){ioff=k;break;}
               if(ioff<0)
               {
                  dim[i].indices[dim[i].n_indices]=rank[j].dim[i].indices[0];
                  dim[i].n_indices++;
               }
             }
             if(flag<0)break;
           }
           if(flag<0){printf("data are not match\n");return -2;}  
           int sizex=1;
           for(int i=0;i<ndim;i++)
           if(dim[i].n_indices>0)sizex*=dim[i].n_indices;
           if(sizex!=size0){printf("the indices can't fulfill the matrix..\n");return -3};
           for(int i=0;i<ndim;i++)
           if(dim[i].n_indices<0)dim[i].n_indices*=-1;
     }*/
     
     void general_data_base::load()
     {
         if(data!=NULL)delete [] data;
         data=xqcd_file_read_once(name,&type);
        if(ndim>0)
        {  size=1;
           for(int i=0;i<ndim;i++)
             size*=dim[i].n_indices;}
     }
     
     void general_data_base::load(int name,std::vector<int> &list)
     {
         if(data!=NULL)delete [] data;
         load_type();
                  
     
     }
     
     void general_data_base::save()
     {
         xqcd_file_write_once(name,&type,data);
     }
     
     void general_data_base::copy(general_data_base &res)
     {
       if(&res!=this)
       {
        res.type=type;
        if(this->data!=NULL)
        {
          res.initialize();
          memcpy(res.data,data,size*sizeof(double));
        }
       }
     }
//io
///////////////

     
     double *general_data_base::seek(int n,...)
     {
        va_list pn;
        va_start(pn,n);
        int *ind=new int[ndim];
        ind[0]=n;
        for(int i=1;i<ndim;i++)
           ind[i]=va_arg(pn,int);
        int off=0;size=1;
        for(int i=ndim-1;i>=0;i--)
        {
            off+=size*(ind[i]<dim[i].n_indices)?ind[i]:0;
            size*=dim[i].n_indices;
        }
        delete [] ind;
        return  data+off;
     }
     
     int datatype::find_dim(int name)
     {   int ind=-1;
         for(int i=0;i<ndim;i++)
            if(dim[i].type==name)ind=i;
         return ind;
     }
     
     void general_data_base::set_size(int ind_f,int &size1,int &size2)
     {
         size1=1,size2=1;
         for(int i=0;i<ndim;i++)
         {  
               if(i<ind_f)size1*=dim[i].n_indices;
               if(i>ind_f)size2*=dim[i].n_indices;
         }
     }
     
     int general_data_base::new_data(general_data_base &res,int name,int new_size)
     {
         general_data_base tmp("");
         tmp.type=type;
         int ind=tmp.find_dim(name);
         res.clear_ind();
         for(int i=0;i<tmp.ndim;i++)
         if(i!=ind) res.add_dimension(tmp.dim[i]);
         else
           if(new_size!=0)res.add_dimension(tmp.dim[i].type,new_size);
         res.initialize();
         return ind;
     }
     
     int general_data_base::move_ind_head(general_data_base &res,int name,int move_flag)
     {
          general_data_base tmp("");tmp.type=type;
          int ind=tmp.find_dim(name);
          if(ind<0)
          {
             res.type=tmp.type;return -1;
          }
          if((ind==0&&move_flag!=0)||(ind==ndim-1&&move_flag==0))
          {
             res.type=tmp.type;return 1;
          }
          res.clear_ind();
          if(move_flag>0)res.add_dimension(tmp.dim[ind]);
          for(int i=0;i<tmp.ndim;i++)if(i!=ind)
             res.add_dimension(tmp.dim[i]);          
          if(move_flag==0)res.add_dimension(tmp.dim[ind]);
          return 0;
     }

     general_data_base* general_data_base::set_pointer(general_data_base &res,general_data_base &tmp)
     {
        if(&res==this)
          return &tmp;
        else
          return this;
     }
     
     int general_data_base::move_ind(general_data_base &res,int name,int move_flag)
     {
        general_data_base tmp("");
        int flag=move_ind_head(tmp,name,move_flag);
        if(data==NULL)return flag;
        if(flag==0)
        {
          general_data_base *tmp1=set_pointer(res,tmp);
          filetype type_tmp=tmp.type;
          copy(*tmp1);
          res.type=type_tmp;
          res.initialize();
          int size1,size2,ic=tmp1->find_dim(name),nset=tmp1->dim[ic].n_indices;
          tmp1->set_size(ic,size1,size2);
            if(move_flag!=0)//move to outside
            for(int i=0;i<size1;i++)
            for(int j=0;j<nset;j++)
              memcpy(res.data+(j*size1+i)*size2,tmp1->data+(i*nset+j)*size2,size2*sizeof(double));
            if(move_flag==0)//move to inside
            for(int i=0;i<size1;i++)
            for(int j=0;j<nset;j++)
            for(int k=0;k<size2;k++)
              res.data[(i*size2+k)*nset+j]=tmp1->data[(i*nset+j)*size2+k];
        }
        else copy(res);
        return flag;
     }
     
     int general_data_base::normal()
     {
         general_data_base tmp("");
         int ind_c=find_dim(dim_conf),ind_j=find_dim(dim_jackknife);
         if(ind_c==0)return 1;
         else if(ind_c>0)
         { move_ind(*this,dim_conf,1);return 0;}
         if(ind_j==0)return 1;
         else if(ind_j>0)
         { move_ind(*this,dim_jackknife,1);return 0;}
         return -1;
//         printf("index is not need to move or not found \n");return 1;   
     } 
  
     int general_data_base::jackknife(general_data_base &res)
     {
        int ic=find_dim(dim_jackknife);
        if(ic>=0)return 0;
        ic=find_dim(dim_conf);
         if(ic<0){print_void_ind(dim_conf);return -1;}
        general_data_base tmp(""),*ptmp=set_pointer(res,tmp);
         move_ind(*ptmp,dim_conf,1);
         res.type=ptmp->type;
         int nconf=ptmp->dim[0].n_indices;
         res.dim[0].type=dim_jackknife;
         res.initialize();
         jack(ptmp->data,nconf,size/nconf,size/nconf,res.data);
         return 0;
     }
     
     int general_data_base::antijack(general_data_base &res)
     {
        int ic=find_dim(dim_conf);
        if(ic>=0)return 0;
        ic=find_dim(dim_jackknife);
         if(ic<0){print_void_ind(dim_jackknife);return -1;}
        general_data_base tmp(""),*ptmp=set_pointer(res,tmp);
         move_ind(*ptmp,dim_jackknife,1);
         res.type=ptmp->type;
         int nconf=ptmp->dim[0].n_indices;
         res.dim[0].type=dim_conf;
         res.initialize();
         anti_jack(ptmp->data,nconf,size/nconf,size/nconf,res.data);
         return 0;
     }
     
     int general_data_base::sum(general_data_base &res,int name)
     {
        int ic=find_dim(name);
        if(ic<0) return -1;
        if(dim[ic].n_indices==1)
        {
           if(&res==this)return 1;
           copy(res);return 0;
        } 
        general_data_base tmp(""),*ptmp=set_pointer(res,tmp);
        copy(*ptmp);
        int size1=1,size2=1,nc=ptmp->dim[ic].n_indices;
        ptmp->set_size(ic,size1,size2);
        res.clear_ind();
        for(int i=0;i<ptmp->ndim;i++)
        if(i!=ic) res.add_dimension(ptmp->dim[i]);
        res.initialize();
         for(int is1=0;is1<size1;is1++)
         for(int ix=0;ix<nc;ix++)
         for(int is2=0;is2<size2;is2++)
            res.data[is1*size2+is2]+=ptmp->data[(is1*nc+ix)*size2+is2];
        return 0;
     }

     void datatype::print_void_ind(int name)
     {
        char ctmp[100];name_i2c(name,ctmp);
        printf("Given dimension %20s is not found\n",ctmp); fflush(stdout);
     }
     
     int general_data_base::aver(general_data_base &res,int name,int flag)
     {
         int ic=find_dim(name);
        if(ic<0){print_void_ind(name); return -1;}
        int size1=1,size2=1,nc=dim[ic].n_indices;
        set_size(ic,size1,size2);
        general_data_base tmp(""),*ptmp=set_pointer(res,tmp);
        copy(*ptmp);
        res.clear_ind();
        int ns=(flag==1)?2:1;
        for(int i=0;i<ptmp->ndim;i++)
        if(i==ic&&flag==1) res.add_dimension(ptmp->dim[i].type,2);
        else if(i!=ic) res.add_dimension(ptmp->dim[i]);
        res.initialize();
         for(int is1=0;is1<size1;is1++)
         for(int ix=0;ix<nc;ix++)
         for(int is2=0;is2<size2;is2++)
            res.data[(is1*ns+0)*size2+is2]+=ptmp->data[(is1*nc+ix)*size2+is2]/nc;
         if(flag!=1) return 0;
         for(int is1=0;is1<size1;is1++)
         for(int ix=0;ix<nc;ix++)
         for(int is2=0;is2<size2;is2++)
            res.data[(is1*2+1)*size2+is2]+=
               pow(ptmp->data[(is1*nc+ix)*size2+is2]-res.data[(is1*2+0)*size2+is2],2);
         int fac=(name==dim_jackknife)?((nc-1)*(nc-1)):1;
         if(name==dim_bootstrap)fac=nc;
         for(int is1=0;is1<size1;is1++)
         for(int is2=0;is2<size2;is2++)
            res.data[(is1*2+1)*size2+is2]=sqrt(res.data[(is1*2+1)*size2+is2]*fac/(nc*(nc-1)));
         return 0;
     }
     
     int general_data_base::nind_cfg()
     {
          int name=-1;
          if(find_dim(dim_conf)>=0) name=dim_conf;
          else
          if(find_dim(dim_jackknife)>=0) name=dim_jackknife;
          else
          if(find_dim(dim_bootstrap)>=0) name=dim_bootstrap;
          if(name==-1)
          {
              printf("The dimension for cfg is not found\n");
              return 1;
          }
          return nind(name);
         
     }
     
     int general_data_base::aver_cfg(general_data_base &res,int aver_flag)
     {
          int name=-1;
          if(find_dim(dim_conf)>=0) name=dim_conf;
          else
          if(find_dim(dim_jackknife)>=0) name=dim_jackknife;
          else
          if(find_dim(dim_bootstrap)>=0) name=dim_bootstrap;
          if(name!=-1)
          {
             if(aver_flag==1)aver(res,name);
             else copy(res);
             res.move_ind(res,name,1);
          }
          else
             return 1;
          return 0;
     }
     
     int general_data_base::make_cfg_bin(general_data_base &res,int bin_size)
     {
         int ic=find_dim(dim_conf);   
         if(ic<0) 
         {
             printf("Dimension conf is not found\n");return 1;
         }
         int size1=1,size2=1,nc=dim[ic].n_indices;
         set_size(ic,size1,size2);
         int nc_bin=nc/bin_size,nc_new=nc_bin*bin_size;
         printf("nc_bin=%4d,nc_new=%4d\n",nc_bin,nc_new);
         fflush(stdout);
         
         
         pick(res,dim_conf,0,nc_new-1);

         general_data_base tmp("");
         tmp.type=res.type;tmp.dim[ic].n_indices=nc_bin;
         tmp.initialize();
         for(int is1=0;is1<size1;is1++)
         for(int ix=0;ix<nc_bin;ix++)
         for(int ix2=0;ix2<bin_size;ix2++)
         for(int is2=0;is2<size2;is2++)
            tmp.data[(is1*nc_bin+ix)*size2+is2]+=
                res.data[((is1*nc_bin+ix)*bin_size+ix2)*size2+is2]/bin_size;
                
         tmp.copy(res);
         
         return 0;
     }
         
     int general_data_base::aver_cfg_bin(general_data_base &res,int bin_size)
     {
         int flag=make_cfg_bin(res,bin_size);
         if(flag==1) return 1;
         res.aver(res,dim_conf);
     }
     

     int datatype::remove_ind(int name)
     {
         int ind_f=find_dim(name);
         if(ind_f<0) {print_void_ind(name);return -1;}
         if(dim[ind_f].n_indices==1)
         {
             for(int i=ind_f;i<ndim;i++)
                 type.head.dimensions[i]=dim[i+1];
             ndim-=1;
             return 0;
         }     
         else return 1;
     }
     
     int general_data_base::pick(general_data_base &res,int name,int ist,int ied)
     {
         int ind_f=find_dim(name);
         if(ind_f<0) {print_void_ind(name);return -1;}
         if(ied==INF)
             pick(res,name,1,&ist,0);
         else
         {
             int it2=ied;
             if(ied>=dim[ind_f].n_indices)it2=dim[ind_f].n_indices-1;
             if(ist<0||ist>it2) return -1;
             int *list=new int[it2-ist+1];
             for(int i=0;i<=it2-ist;i++)
               list[i]=dim[ind_f].indices[i+ist];
             int flag=pick(res,name,it2-ist+1,list,0);
             delete [] list;
             return flag;
         }
         return 0;
     }
     
     int general_data_base::pick(general_data_base &res,int name,int sizet,int* list,int flag)
     {
        char file[500];name_i2c(name,file);
         //printf("picking dim %10s....\n",file);fflush(stdout);
         int *ind=new int[sizet],inf=find_dim(name);
         if(inf<0){print_void_ind(name);return -1;}
         general_data_base tmp(""),*ptmp=set_pointer(res,tmp);copy(*ptmp);
//         printf("set_pointer\n");
         int size1,size2;ptmp->set_size(inf,size1,size2);
         for(int i=0;i<sizet;i++)ind[i]=-1;
         for(int i=0;i<sizet;i++)
          for(int j=0;j<ptmp->dim[inf].n_indices;j++)
           if(list[i]==ptmp->dim[inf].indices[j]){ind[i]=j;break;}
         if(ind[sizet-1]==-1){delete [] ind;
         for(int i=0;i<sizet;i++) printf("%10d",list[i]);
         printf(" list\n");
         for(int i=0;i<ptmp->dim[inf].n_indices;i++) printf("%10d",ptmp->dim[inf].indices[i]);
         printf(" ptmp\n");
         for(int i=0;i<sizet;i++) printf("%10d",ind[i]);
         printf(" ind\n");
           printf("element(s) in the list is not found\n");return -1;}
         res.clear_ind();
         if(flag==1)res.add_dimension(ptmp->dim[inf].type,sizet,list);
         for(int i=0;i<ptmp->ndim;i++)
         if(i!=inf)res.add_dimension(ptmp->dim[i]);
         else if(flag!=1)res.add_dimension(ptmp->dim[inf].type,sizet,list);
//         printf("before\n");
         res.initialize();
//         printf("after\n");
         if(flag==1)
          for(int i=0;i<sizet;i++)
          for(int j=0;j<size1;j++)
            memcpy(res.data+size2*(i*size1+j),ptmp->data+size2*(j*ptmp->dim[inf].n_indices+ind[i])
              ,size2*sizeof(double));
         else
          for(int j=0;j<size1;j++)
          for(int i=0;i<sizet;i++)
            memcpy(res.data+size2*(j*sizet+i),ptmp->data+size2*(j*ptmp->dim[inf].n_indices+ind[i])
              ,size2*sizeof(double));
        delete [] ind;
        return 0;
     }     

     int general_data_base::exclude(general_data_base &res,int name,int index)
     {
        int ind=find_dim(name);
        if(ind<0){print_void_ind(name);return -1;}
        int off=-1;
        for(int i=0;i<dim[ind].n_indices;i++)
           if(dim[ind].indices[i]==index){off=i;break;}
        if(off<0)
        {copy(res);return 1;}
        int size1,size2;
        set_size(ind,size1,size2);
        general_data_base tmp(""),*ptmp=set_pointer(res,tmp);copy(*ptmp);
        res.clear_ind();
        res.type=ptmp->type;
        for(int i=off;i<res.dim[ind].n_indices-1;i++)
          res.dim[ind].indices[i]=res.dim[ind].indices[i+1];
        res.dim[ind].n_indices-=1;
        res.initialize();
        double *p=res.data;
        for(int i=0;i<size1;i++)
        for(int j=0;j<ptmp->dim[ind].n_indices;j++)
        if(j!=off)
        {    memcpy(p,ptmp->data+(i*ptmp->dim[ind].n_indices+j)*size2,size2*sizeof(double));
             p+=size2;}
        return 0;
     }
     
     int general_data_base::folding(general_data_base &res,int* foldfac,int nt2,int its)
     {
          int inds[3],nte=nt2/2;
          general_data_base tmp(""),*ptmp=set_pointer(res,tmp);copy(*ptmp);
          inds[0]=ptmp->new_data(res,dim_t,nte);
          ptmp->set_size(inds[0],inds[1],inds[2]);
          int nt=ptmp->nind(dim_t),ne=nt/nte;
          for(int it=0;it<nt;it++)
          {
             int ipart=((it-its+nt)/nte)%(nt/nte);
             int ieo=1-2*(ipart%2);
             int itp=(nt+ieo*(it-its))%nte;
            for(int i=0;i<inds[1];i++)
            for(int j=0;j<inds[2];j++)
                  res.data[(i*nte+itp)*inds[2]+j]+=ptmp->data[(i*nt+it)*inds[2]+j]/ne*foldfac[it];
          }
          fflush(stdout);
          return 0;
     }


     int general_data_base::eff_mass(general_data_base &res,int fold_type,int nt2,int it_st)
     {
          int inds[3];
          int nt=dim[find_dim(dim_t)].n_indices;
          inds[0]=new_data(res,dim_t,nt);
          set_size(inds[0],inds[1],inds[2]);
          for(int i=0;i<inds[1];i++)
          for(int it=1;it<nt;it++)
          for(int j=0;j<inds[2];j++)
          {
             double min=0.0,max=10.0;

             int it0=(it-it_st-1+nt2*2)%(nt2*2);
             res.data[(i*nt+it)*inds[2]+j]=effective_mass(nt2,it0,data[(i*nt+it-1)*inds[2]+j],data[(i*nt+it)*inds[2]+j]
		,fold_type);
	  }
         return 0;
     }
     
     int general_data_base::combine(general_data_base &file1,general_data_base &file2,bool sort)
     {
         if(file1.data==NULL){file2.copy(*this);return 1;}
         if(file2.data==NULL){file1.copy(*this);return 1;}
     
         int index=check_index(file1,file2);
         if(index<0) return -1;
         int size1,size2;
         if(index==file1.ndim)
         {
            file1.normal();file2.normal();
            for(int i=0;i<file2.dim[0].n_indices;i++)
            file2.dim[0].indices[i]+=100000;
            size1=1;size2=file2.size/file2.dim[0].n_indices;
            index=0;
         }
         else
         {
             file1.set_size(index,size1,size2);
         }
    
         int nd1=file1.dim[index].n_indices;
         int nd2=file2.dim[index].n_indices;
    
         type=file1.type;
         dim[index].n_indices+=nd2;
         for(int i=0;i<nd2;i++)
            dim[index].indices[i+nd1]=
               file2.dim[index].indices[i];
         initialize();
         int *map=new int[nd1+nd2];
         sort_ind(index,map,sort);
    
         for(int i=0;i<size1;i++)
         for(int j=0;j<nd1+nd2;j++)
         if(map[j]<nd1)
            memcpy(data+size2*((nd1+nd2)*i+j),file1.data+size2*(nd1*i+map[j]),sizeof(double)*size2);
         else
            memcpy(data+size2*((nd1+nd2)*i+j),file2.data+size2*(nd2*i+map[j]-nd1),sizeof(double)*size2);
         delete [] map;
         return 0;    
     }
     
     int datatype::sort_ind(int ind,int *map,bool sort) 
     {
        for(int i=0;i<dim[ind].n_indices;i++)map[i]=i;
        if(sort==false)return 1;
        int t;
        for(int i=0;i<dim[ind].n_indices;i++)
        for(int j=i+1;j<dim[ind].n_indices;j++)
        if(dim[ind].indices[i]>dim[ind].indices[j])
        {
           t=map[i];map[i]=map[j];map[j]=t; 
           t=dim[ind].indices[i];dim[ind].indices[i]=dim[ind].indices[j];
           dim[ind].indices[j]=t;
        }
        return 0;
     }


      int general_data_base::bootstrap_sample(general_data_base &res,const int_rand& ran)
      {
        int ic=find_dim(dim_conf);
        if(ic<0){print_void_ind(dim_conf);return -1;}
        if(&res==this){printf("You can't overwrite the original data with bootstrap sample\n");return 1;}
        int ncfg=nind(dim_conf);
        res.type=type;
        res.initialize();
        int size1,size2;
        set_size(ic,size1,size2);

        int r[ncfg];
        for(int i=0;i<ncfg;i++) r[i]=ran.i_rand(ncfg);

        for(int i=0;i<size1;i++)
        for(int j=0;j<ncfg;j++)
        for(int k=0;k<size2;k++)
           res.data[(i*ncfg+j)*size2+k]=data[(i*ncfg+r[j])*size2+k];
        return 0;
      }

      int general_data_base::bootstrap(general_data_base &res,int nboot,const int_rand& r)
      {
        int ic=find_dim(dim_conf);
        if(ic<0){print_void_ind(dim_conf);return -1;}

        general_data_base tmp(""),*ptmp=set_pointer(res,tmp);copy(*ptmp);
        res.type=ptmp->type;
        int nconf=ptmp->dim[ic].n_indices;
        res.dim[ic].type=dim_bootstrap;
        res.dim[ic].n_indices=nboot;
        for(int i=0;i<nboot;i++)res.dim[ic].indices[i]=i;
        res.initialize();
//        res.print();

        int size1,size2;
        set_size(ic,size1,size2);
        
        general_data_base tmp2("");
        for(int i=0;i<nboot;i++)
        {
            ptmp->bootstrap_sample(tmp2,r);
            tmp2.aver(tmp2,dim_conf,0);
            fflush(stdout);
            for(int is=0;is<size1;is++)
            for(int k=0;k<size2;k++)
              res.data[(is*nboot+i)*size2+k]=tmp2.data[is*size2+k];
        }
        return 0;
      }
     
/*     void general_data_base::move_inds(general_data_base &res,int move_flag,int ncount,...)
     {
        va_list pn;
        va_start(pn,n);
        int *ind[]
     
     
        va_end(pn);
     }*/


	void		data_iterator::sort_indices(int *pind, int nind)
	{
		int	i, j, t;

		for(i=0; i<nind; i++)
			for(j=i+1; j<nind; j++)
				if(pind[i]>pind[j])
				{
					t = pind[i];
					pind[i] = pind[j];
					pind[j] = t;
				}
	}

	void		data_iterator::copy_data()
	{
		int	i;
		char	localbuf[N_MAX_STRING];

		for(i=0; i<pdata->size; i++)
		{
			type_disp2index(&pdata->type, runningdims+nrunningdims, i);
			pdata->data[i] = rawdata.data[type_index2disp(&rawdata.type, runningdims)];
		}

		strcpy(namebuf, "iter");
		for(i=0; i<nrunningdims; i++)
		{
			sprintf(localbuf, "_%s%d", xqcd_type_dim_desc[runningdims[i][0]], runningdims[i][1]);
			strcat(namebuf, localbuf);
		}
		strcpy(pdata->name, namebuf);
	}


	data_iterator	data_iterator::operator++()
	{
		int	i, j, k;

		for(i=nrunningdims-1; i>=0; i--)
		{
			for(j=0; j<rawdata.ndim; j++)
				if(rawdata.dim[j].type == runningdims[i][0])
					break;
			if(rawdata.dim[j].indices[rawdata.dim[j].n_indices-1] == runningdims[i][1])
			{
				runningdims[i][1] = rawdata.dim[j].indices[0];
				continue;
			}
			else
			{
				for(k=0; k<rawdata.dim[j].n_indices-1; k++)
					if(rawdata.dim[j].indices[k] == runningdims[i][1])
					{
						runningdims[i][1] = rawdata.dim[j].indices[k+1];
						break;
					}
				break;
			}
		}

		if(i<0)
		{
			delete pdata;
			pdata = NULL;
		}
		else
			copy_data();

		return *this;
	}


int check_dim(one_dim &a,one_dim &b)
{

         if(a.n_indices==b.n_indices)
         {
            int ic=0;
            for(int m=0;m<a.n_indices;m++)
            if(a.indices[m]!=b.indices[m])
                ic++;
            if(ic==0) return 0;
            if(ic==a.n_indices) return 1;//supposing all the index are in ordering.
            //if not those cases
               return -1;
          }
          else
          {
             for(int m=0;m<a.n_indices;m++)
             for(int n=0;n<b.n_indices;n++)
              if(a.indices[m]==b.indices[n]) return -2;
             return 1;
           }
}

int check_index(datatype &a, datatype &b)
{

     if(a.ndim!=b.ndim)
     {printf("Imposible to combine files with different n_dim, exit\n");return -1;}

     int icount=0,index=a.ndim;
     for(int k=0;k<a.ndim;k++)
     if(a.dim[k].type==b.dim[k].type)
     {
         int ind=check_dim(a.dim[k],b.dim[k]);

         if(ind==-1)
         {printf("index of dimension %d is not match, exit\n",k);
          for(int i=0;i<a.dim[k].n_indices;i++)
            printf("%5d",a.dim[k].indices[i]); printf("\n");
          for(int i=0;i<b.dim[k].n_indices;i++)
            printf("%5d",b.dim[k].indices[i]); printf("\n");
               return -1;}
         if(ind==-2)
         {printf("index of dimension %d is not doubleing, exit\n",
                   k);return -1;}
         icount+=ind;
         if(ind==1)index=k;
         if(icount==2){printf("more than one dimension is different, exit\n");return -1;}
     }
     else
     {
//         printf("type of dimension %4d is not match (%6d,%6d), exit\n",k,
//            a.dim[k].type,b.dim[k].type);
          char name1[100],name2[100];
          a.name_i2c(a.dim[k].type,name1);
          a.name_i2c(b.dim[k].type,name2);
         printf("type of dimension %4d is not match (%20s,%20s), exit\n",k,
             name1,name2);
     }

     return index;
     //return the index with different indices. return a.ndim if the datetype is totally the same

}


