
#include	<stdlib.h>
#include	<string.h>
#include	"io_general.h"

#ifdef	MPIIO	// if in the milc MPI environment
#include	"lattice.h"
#include	"io_utils.h"
#include	"make_source.h"
#endif

#ifndef	LITTLE_ENDIAN
#define	LITTLE_ENDIAN	__LITTLE_ENDIAN
#define	BIG_ENDIAN	__BIG_ENDIAN
#endif


char	*xqcd_type_dim_desc[N_MAX_DIMENSION_TYPES]={"other", "x", "y", "z", "t", "d", "c", "d2", "c2", "complex", "mass", "smear", "displacement", "s_01", "s_02", "s_03", "s_11", "s_12", "s_13", "d_01", "d_02", "d_03", "d_11", "d_12", "d_13", "conf", "operator", "momentum", "direction", "t2", "mass2", "column", "row", "temporary", "temporary2", "temporary3", "temporary4", "errorbar", "operator2", "param", "fit_left", "fit_right", "jackknife", "jackknife2", "jackknife3", "jackknife4", "summary", "channel", "channel2", "eigen", "d_row", "d_col", "c_row", "c_col", "parity", "noise", "evenodd", "disp_x", "disp_y", "disp_z", "disp_t", "t3", "t4", "t_source", "t_current", "t_sink","bootstrap", "nothing"};


int get_endian()
{                       
	int iTest;            
	char *pc;     

	iTest = 1;
	pc = (char*) &iTest;
	if ((*(pc)) == 1)
	{     
		return LITTLE_ENDIAN;
	}     
	else  
	{ 
		printf("Warning : The program is only tested on x86-64 environment and may have issues on other machines : the xqcdio data file is assumed to be in LITTLE_ENDIAN.\n");
		return BIG_ENDIAN;
	}                     
} 

void	switch_endian_int(int *pdata, int n)
{
	int	i, j;
	char	*p1, *p2;
	char	cbuf;

	for(i=0; i<n; i++)
		for(j=0; j<sizeof(int)/2; j++)
		{
			cbuf=((char*)(pdata+i))[j];
			((char*)(pdata+i))[j]=((char*)(pdata+i))[sizeof(int)-1-j];
			((char*)(pdata+i))[sizeof(int)-1-j]=cbuf;
		}

	return;
}

void	switch_endian_double(double *pdata, int n)
{
	int	i, j;
	char	*p1, *p2;
	char	cbuf;

	for(i=0; i<n; i++)
		for(j=0; j<sizeof(double)/2; j++)
		{
			cbuf=((char*)(pdata+i))[j];
			((char*)(pdata+i))[j]=((char*)(pdata+i))[sizeof(double)-1-j];
			((char*)(pdata+i))[sizeof(double)-1-j]=cbuf;
		}

	return;
}

int*	make_endian_switch_buffer_int(int *pdata, int n)
{
	int	*buf;
	int	i, j;
	char	*p1, *p2;

	if((buf=(int *)malloc(sizeof(int)*n))==NULL)
	{
		printf("An error occured when allocing memory.\n");
		return NULL;
	}

	for(i=0; i<n; i++)
		for(j=0; j<sizeof(int); j++)
			((char*)(buf+i))[j]=((char*)(pdata+i))[sizeof(int)-1-j];

	return buf;
}

double*	make_endian_switch_buffer_double(double *pdata, int n)
{
	double	*buf;
	int	i, j;
	char	*p1, *p2;

	if((buf=(double *)malloc(sizeof(double)*n))==NULL)
	{
		printf("An error occured when allocing memory.\n");
		return NULL;
	}

	for(i=0; i<n; i++)
		for(j=0; j<sizeof(double); j++)
			((char*)(buf+i))[j]=((char*)(pdata+i))[sizeof(double)-1-j];

	return buf;
}


int	xqcd_file_write_once(char *name, filetype *type, double *pdata)
{
	FILE	*fp;
	int	count;
	int	i;
	int	*pbufint;
	double	*pbufdouble;

	if((fp=fopen(name, "wb"))==NULL)
	{
		printf("Failed opening file %s for write.\n", name);
		return	-1;
	}

	if(get_endian()==LITTLE_ENDIAN)
	{
		if(fwrite(type, sizeof(filetype), 1, fp)!=1)
		{
			printf("Failed writing the head of file %s.\n", name);
			return	-1;
		}
	}
	else
	{
		printf("Converting to little endian\n");

		pbufint=make_endian_switch_buffer_int((int*)type, sizeof(filetype)/sizeof(int));
		if(fwrite(pbufint, sizeof(filetype), 1, fp)!=1)
		{
			printf("Failed writing the head of file %s.\n", name);
			return	-1;
		}
		free(pbufint);
	}

	count=1;
	for(i=0; i<type->head.n_dimensions; i++)
		count*=type->head.dimensions[i].n_indices;
	if(get_endian()==LITTLE_ENDIAN)
	{
		if(fwrite(pdata, sizeof(double), count, fp)!=count)
		{
			printf("Failed writing the data of file %s.\n", name);
			return	-1;
		}
	}
	else
	{
		pbufdouble=make_endian_switch_buffer_double(pdata, count);
		if(fwrite(pbufdouble, sizeof(double), count, fp)!=count)
		{
			printf("Failed writing the data of file %s.\n", name);
			return	-1;
		}
		free(pbufdouble);
	}

	fclose(fp);

	return 0;
}


double*	xqcd_file_read_once(char *name, filetype *type)
{
	FILE	*fp;
	double	*pdata;
	int	count, i;

	if((fp=fopen(name, "rb"))==NULL)
	{
		printf("Failed opening file %s for read.\n", name);
		return	NULL;
	}

	if(fread(type, sizeof(filetype), 1, fp)!=1)
	{
		printf("Failed reading the head of file %s.\n", name);
		return	NULL;
	}

	if(get_endian()!=LITTLE_ENDIAN)
	{
		// debug
		printf("Converting to little endian\n");

		switch_endian_int((int*)type, sizeof(filetype)/sizeof(int));
	}

	count=1;
	for(i=0; i<type->head.n_dimensions; i++)
		count*=type->head.dimensions[i].n_indices;

	if((pdata=(double *)malloc(sizeof(double)*count))==NULL)
	{
		printf("Failed to allocating memory for %d doubles.\n", count);
		return NULL;
	}

	if(fread(pdata, sizeof(double), count, fp)!=count)
	{
		printf("Failed reading the data of file %s.\n", name);
		return	NULL;
	}

	if(get_endian()!=LITTLE_ENDIAN)
		switch_endian_double(pdata, count);

	fclose(fp);

	return pdata;
}


//TODO:
//filehandle	xqcd_file_open(char *name, filetype *type);
//TODO:
//int	xqcd_file_close(filehandle file);
//TODO:
//int	xqcd_file_write(filehandle file, filetype *dataset, double *pdata);
//TODO:
//double*	xqcd_file_read(filehandle file, filetype *dataset);



#ifdef MPIIO


filehandle      xqcd_file_open_write_mpi_milc(char *name, filetype *type, gridtype *pgrid)
{
	filehandle	handle;
	MPI_Offset	offset;
	MPI_Status	status;
//	char		errorbuf[N_MAXERRORBUF];
//	int		errorno;
	int		iret;
	int		i, j;


	node0_printf("Saving file %s ...\n", name);
	if(this_node==0)
		xqcd_type_print(type);

	if((handle=malloc(sizeof(filehandlebase)))==NULL)
	{
		printf("Memory allocation failed when allocating a file handle.\n");
		return NULL;
	}

	strcpy(handle->name, name);
	if(iret=MPI_File_open(MPI_COMM_WORLD, name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &(handle->fp_mpi)))
	{
		printf("Error opening file %s for mpi writing with error code #%d.\n", name, iret);
		return NULL;
	}

	memcpy(&handle->type, type, sizeof(filetype));
	for(handle->rundimensions=0; handle->type.head.dimensions[handle->rundimensions].type != dim_x && handle->type.head.dimensions[handle->rundimensions].type != dim_y && handle->type.head.dimensions[handle->rundimensions].type != dim_z && handle->type.head.dimensions[handle->rundimensions].type != dim_t; handle->rundimensions++);

	//debug
	node0_printf("There are %d run dimensions.\n", handle->rundimensions);

//	handle->blocksize = sizeof(double);
//	for(i=handle->rundimensions; i<handle->type.head.n_dimensions; i++)
//		handle->blocksize *= handle->type.head.dimensions[i].n_indices;

	if(get_endian()==LITTLE_ENDIAN)
		handle->doubletype = MPI_DOUBLE;
	else
	{
		// TODO: make a type for endian switch.
		printf("endian switch not implemented.\n");
		return NULL;
	}

	if(this_node==0)
	{
		offset = 0;
		if(get_endian()==LITTLE_ENDIAN)
			iret = MPI_File_write_at(handle->fp_mpi, offset, &handle->type, sizeof(filetype)/sizeof(int), MPI_INT, &status);
		else
		{
			// TODO: make a type for endian switch.
			printf("endian switch not implemented.\n");
			return NULL;
		}
	}

	int	count;

	count = 1;
	for(i=0; i<handle->type.head.n_dimensions; i++)
		if(handle->type.head.dimensions[i].type == dim_x || handle->type.head.dimensions[i].type == dim_y || handle->type.head.dimensions[i].type == dim_z || handle->type.head.dimensions[i].type == dim_t)
			count = 1;
		else
			count *= handle->type.head.dimensions[i].n_indices;
	MPI_Type_contiguous(count, handle->doubletype, &handle->basetype);
	MPI_Type_commit(&handle->basetype);

	//debug
	node0_printf("The fixed block is %d doubles.\n", count);

	int	*ones, *dispmem, *dispfile;
	int	nmemfile;
	site	*psSite;
	int	flag;

	if((ones=malloc(sizeof(int)*sites_on_node))==NULL || (dispmem=malloc(sizeof(int)*sites_on_node))==NULL || (dispfile=malloc(sizeof(int)*sites_on_node))==NULL)
	{
		printf("Memory allocation failed when allocating variables : ones, dispmem, dispfile.\n");
		return NULL;
	}
	for(i=0; i<sites_on_node; i++)
		ones[i] = 1;
	nmemfile = 0;

	FORALLSITES (dispmem[nmemfile], psSite)
	{
		flag = 0;
		dispfile[nmemfile] = 0;
		for(i=handle->rundimensions; i<handle->rundimensions+4; i++)
		{
			if(i>handle->rundimensions)
				dispfile[nmemfile] *= handle->type.head.dimensions[i].n_indices;
			switch(handle->type.head.dimensions[i].type)
			{
				case dim_x :
					for(j=0; j<handle->type.head.dimensions[i].n_indices && psSite->x != handle->type.head.dimensions[i].indices[j]; j++);
					if(j<handle->type.head.dimensions[i].n_indices)
					{
						flag |= 1;
						dispfile[nmemfile] += j;
					}
					break;
				case dim_y :
					for(j=0; j<handle->type.head.dimensions[i].n_indices && psSite->y != handle->type.head.dimensions[i].indices[j]; j++);
					if(j<handle->type.head.dimensions[i].n_indices)
					{
						flag |= 2;
						dispfile[nmemfile] += j;
					}
					break;
				case dim_z :
					for(j=0; j<handle->type.head.dimensions[i].n_indices && psSite->z != handle->type.head.dimensions[i].indices[j]; j++);
					if(j<handle->type.head.dimensions[i].n_indices)
					{
						flag |= 4;
						dispfile[nmemfile] += j;
					}
					break;
				case dim_t :
					for(j=0; j<handle->type.head.dimensions[i].n_indices && psSite->t != handle->type.head.dimensions[i].indices[j]; j++);
					if(j<handle->type.head.dimensions[i].n_indices)
					{
						flag |= 8;
						dispfile[nmemfile] += j;
					}
					break;
				default :
					printf("Warning! The dimensions of the given type is not correct!\n");
			}
		}
		if(flag == 15 && is_on_grid(pgrid, psSite->x, psSite->y, psSite->z, psSite->t))
		{
			//debug
//			printf("node#%d: nmemfile=%d, (%d,%d,%d,%d), dispfile=%d, dispmem=%d.\n", this_node, nmemfile, psSite->x, psSite->y, psSite->z, psSite->t, dispfile[nmemfile], dispmem[nmemfile]);

			nmemfile++;
			dispmem[nmemfile] = dispmem[nmemfile-1];
		}
	}

	MPI_Type_indexed(nmemfile, ones, dispmem, handle->basetype, &handle->memdatatype);
	MPI_Type_commit(&handle->memdatatype);
	MPI_Type_indexed(nmemfile, ones, dispfile, handle->basetype, &handle->filedatatype);
	MPI_Type_commit(&handle->filedatatype);

	free(ones);
	free(dispfile);
	free(dispmem);


	return handle;
}

void	xqcd_file_close_mpi_milc(filehandle file)
{
	char            errorbuf[N_MAXERRORBUF];
	int             errorno;
	int             iret;

	iret=MPI_File_close(&file->fp_mpi);

	if(get_endian()!=LITTLE_ENDIAN)
		MPI_Type_free(&file->doubletype);
	MPI_Type_free(&file->basetype);
	MPI_Type_free(&file->memdatatype);
	MPI_Type_free(&file->filedatatype);

	free(file);

	return;
}

void	collect_dims_tzyx(filetype *type, gridtype *pgrid, int rundims)
{
	int	t, z, y, x;
	int	i;

	type->head.dimensions[rundims].type=dim_t;
	type->head.dimensions[rundims].n_indices=0;
	type->head.dimensions[rundims+1].type=dim_z;
	type->head.dimensions[rundims+1].n_indices=0;
	type->head.dimensions[rundims+2].type=dim_y;
	type->head.dimensions[rundims+2].n_indices=0;
	type->head.dimensions[rundims+3].type=dim_x;
	type->head.dimensions[rundims+3].n_indices=0;

	for(t=0; t<nt; t++)
	for(z=0; z<nz; z++)
	for(y=0; y<ny; y++)
	for(x=0; x<nx; x++)
		if(is_on_grid(pgrid, x, y, z, t))
		{
			for(i=0; i<type->head.dimensions[rundims].n_indices; i++)
				if(type->head.dimensions[rundims].indices[i]==t)
					break;
			if(i==type->head.dimensions[rundims].n_indices)
			{
				type->head.dimensions[rundims].n_indices++;
				type->head.dimensions[rundims].indices[i] = t;
			}
			for(i=0; i<type->head.dimensions[rundims+1].n_indices; i++)
				if(type->head.dimensions[rundims+1].indices[i]==z)
					break;
			if(i==type->head.dimensions[rundims+1].n_indices)
			{
				type->head.dimensions[rundims+1].n_indices++;
				type->head.dimensions[rundims+1].indices[i] = z;
			}
			for(i=0; i<type->head.dimensions[rundims+2].n_indices; i++)
				if(type->head.dimensions[rundims+2].indices[i]==y)
					break;
			if(i==type->head.dimensions[rundims+2].n_indices)
			{
				type->head.dimensions[rundims+2].n_indices++;
				type->head.dimensions[rundims+2].indices[i] = y;
			}
			for(i=0; i<type->head.dimensions[rundims+3].n_indices; i++)
				if(type->head.dimensions[rundims+3].indices[i]==x)
					break;
			if(i==type->head.dimensions[rundims+3].n_indices)
			{
				type->head.dimensions[rundims+3].n_indices++;
				type->head.dimensions[rundims+3].indices[i] = x;
			}
		}
}

void	make_filetype_propagators(filetype *type, gridtype *pgrid)
{
	int	i;

	type->head.n_dimensions=9;

	type->head.dimensions[0].type=dim_c;
	type->head.dimensions[0].n_indices=SOURCE_COLOR;
	for(i=0; i<SOURCE_COLOR; i++)
		type->head.dimensions[0].indices[i]=i;

	type->head.dimensions[1].type=dim_d;
	type->head.dimensions[1].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[1].indices[i]=i;

	type->head.dimensions[6].type=dim_d2;
	type->head.dimensions[6].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[6].indices[i]=i;

	type->head.dimensions[7].type=dim_c2;
	type->head.dimensions[7].n_indices=SOURCE_COLOR;
	for(i=0; i<SOURCE_COLOR; i++)
		type->head.dimensions[7].indices[i]=i;

	type->head.dimensions[8].type=dim_complex;
	type->head.dimensions[8].n_indices=2;
	for(i=0; i<2; i++)
		type->head.dimensions[8].indices[i]=i;

	collect_dims_tzyx(type, pgrid, 2);
}

void	make_filetype_loop_x(filetype *type, gridtype *pgrid)
{
	int	i;

	type->head.n_dimensions=9;

	type->head.dimensions[0].type=dim_c;
	type->head.dimensions[0].n_indices=SOURCE_COLOR;
	for(i=0; i<SOURCE_COLOR; i++)
		type->head.dimensions[0].indices[i]=i;

	type->head.dimensions[1].type=dim_d;
	type->head.dimensions[1].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[1].indices[i]=i;

	type->head.dimensions[6].type=dim_d_11;
	type->head.dimensions[6].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[6].indices[i]=i;

	type->head.dimensions[7].type=dim_d_12;
	type->head.dimensions[7].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[7].indices[i]=i;

	type->head.dimensions[8].type=dim_complex;
	type->head.dimensions[8].n_indices=2;
	for(i=0; i<2; i++)
		type->head.dimensions[8].indices[i]=i;

	collect_dims_tzyx(type, pgrid, 2);
}

void	make_filetype_loop_x2(filetype *type, gridtype *pgrid)
{
	int	i;

	type->head.n_dimensions=10;

	type->head.dimensions[0].type=dim_c;
	type->head.dimensions[0].n_indices=SOURCE_COLOR;
	for(i=0; i<SOURCE_COLOR; i++)
		type->head.dimensions[0].indices[i]=i;

	type->head.dimensions[1].type=dim_d;
	type->head.dimensions[1].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[1].indices[i]=i;

	type->head.dimensions[6].type=dim_d_11;
	type->head.dimensions[6].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[6].indices[i]=i;

	type->head.dimensions[7].type=dim_d_12;
	type->head.dimensions[7].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[7].indices[i]=i;

	type->head.dimensions[8].type=dim_d_13;
	type->head.dimensions[8].n_indices=SOURCE_SPIN;
	for(i=0; i<SOURCE_SPIN; i++)
		type->head.dimensions[8].indices[i]=i;

	type->head.dimensions[9].type=dim_complex;
	type->head.dimensions[9].n_indices=2;
	for(i=0; i<2; i++)
		type->head.dimensions[9].indices[i]=i;

	collect_dims_tzyx(type, pgrid, 2);
}


int	xqcd_file_write_mpi_milc_begin(filehandle file, double *vec, int *rundims)
{
	int	dims[N_MAX_DIMENSIONS][2];
	int	disp;
	int	i;

	for(i=0; i<file->type.head.n_dimensions; i++)
	{
		dims[i][0] = file->type.head.dimensions[i].type;
		if(i<file->rundimensions)
			dims[i][1] = rundims[i];
		else
			dims[i][1] = file->type.head.dimensions[i].indices[0];
	}
	disp = N_MAX_HEADLENGTH + type_index2disp(&file->type, dims)*sizeof(double);

	MPI_File_set_view(file->fp_mpi, disp, file->basetype, file->filedatatype, "native", MPI_INFO_NULL);
	MPI_File_write_all_begin(file->fp_mpi, vec, 1, file->memdatatype);
}

int	xqcd_file_write_mpi_milc_end(filehandle file, double *vec)
{
	MPI_Status	status;

	MPI_File_write_all_end(file->fp_mpi, vec, &status);
}














#endif



void		xqcd_type_print(filetype *type)
{
	int	i, j;

#ifdef	MPIIO
	printf("The information of the file type (on node# %d): \n", this_node);
#endif
	printf("There are %d dimensions.\n", type->head.n_dimensions);

	for(i=0;i<type->head.n_dimensions;i++)
	{
		printf("Dimension #%d : type %s : %d indices :\n", i, xqcd_type_dim_desc[type->head.dimensions[i].type], type->head.dimensions[i].n_indices);
		for(j=0;j<type->head.dimensions[i].n_indices;j++)
			printf("\t%d", type->head.dimensions[i].indices[j]);
		printf("\n");
	}
	printf("\n");
}



int             xqcd_file_print(char *name, char *tag)
{
	filetype	ty;
	double		*dat;
	int		pdim[N_MAX_DIMENSIONS];
	int		i, j, count;

	if((dat=xqcd_file_read_once(name, &ty))==NULL)
	{
		printf("Error reading file.\n");
		return -1;
	}

	xqcd_type_print(&ty);

	count=1;
	if(*tag)
		printf("LABELS_%s\t", tag);
	for(i=0; i<ty.head.n_dimensions; i++)
	{
		pdim[i]=0;
		count*=ty.head.dimensions[i].n_indices;
		printf("%s\t", xqcd_type_dim_desc[ty.head.dimensions[i].type]);
	}
	printf("DATA\n");

	for(i=0; i<count; i++)
	{
		if(*tag)
			printf("%s\t", tag);
		for(j=0; j<ty.head.n_dimensions; j++)
			printf("%d\t", ty.head.dimensions[j].indices[pdim[j]]);
		printf("%20.15le\n", dat[i]);

		j=ty.head.n_dimensions-1;
		pdim[j]++;
		while(pdim[j]==ty.head.dimensions[j].n_indices)
		{
			pdim[j]=0;
			j--;
			pdim[j]++;
		}
	}

	free(dat);

	return 0;
}



#ifdef	MPIIO

gridtype*	set_grid(gridtype *pgrid, int icol)
{
	switch (par_buf.source_type)
	{
		case	SOURCE_Z3:
		case	SOURCE_Z4:
			pgrid->x0=0;
			pgrid->y0=0;
			pgrid->z0=0;
			if(icol<0)
				pgrid->t0=0;
			else
				pgrid->t0=(int)(icol/(SOURCE_COLOR*SOURCE_SPIN*2));
			pgrid->dx=par_buf.source_schema[0];
			pgrid->dy=par_buf.source_schema[1];
			pgrid->dz=par_buf.source_schema[2];
			if(icol==-1)
				pgrid->dt=1;
			else
				pgrid->dt=par_buf.source_schema[3];
			if (icol<0)
				pgrid->eo=3;
			else
				pgrid->eo=2-(int)(icol/(SOURCE_COLOR*SOURCE_SPIN))%2;
			break;

		case	SOURCE_POINT:
			pgrid->x0=par_buf.source_loc[0];
			pgrid->y0=par_buf.source_loc[1];
			pgrid->z0=par_buf.source_loc[2];
			pgrid->t0=par_buf.source_loc[3];
			pgrid->dx=nx;
			pgrid->dy=ny;
			pgrid->dz=nz;
			pgrid->dt=nt;
			pgrid->eo=3;
			break;

		default:
			printf("Wrong grid parameters. Check the input file.\n");
			return	NULL;
	}

	return	pgrid;
}

#endif


int             is_on_grid(gridtype *pgrid, int x, int y, int z, int t)
{
	if(((x-pgrid->x0)%pgrid->dx) || ((y-pgrid->y0)%pgrid->dy) || ((z-pgrid->z0)%pgrid->dz) || ((t-pgrid->t0)%pgrid->dt))
		return	0;

	if(pgrid->eo==3)
		return	-1;

	if((((x-pgrid->x0)/pgrid->dx) + ((y-pgrid->y0)/pgrid->dy) + ((z-pgrid->z0)/pgrid->dz) + ((t-pgrid->t0)/pgrid->dt) + pgrid->eo)%2)
		return	0;
	else
		return	-1;
}

#ifdef	MPIIO

int             n_site_on_grid(gridtype *pgrid, int timeslice)
{
	int	res;
	int	ix, iy, iz, it;

	res=0;

	for(ix=0;ix<nx;ix++)
	for(iy=0;iy<ny;iy++)
	for(iz=0;iz<nz;iz++)
	{
		if(timeslice<0)
		{
			for(it=0;it<nt;it++)
			{
				if(is_on_grid(pgrid, ix, iy, iz, it))
					res++;
			}
		}
		else
		{
			if(is_on_grid(pgrid, ix, iy, iz, timeslice))
				res++;
		}
	}

	return	res;
}
















#endif



int             type_index2disp(filetype *type, int ind[][2])
{
	int	result;
	int	i, j, k;

	result = 0;
	for(i=0; i<type->head.n_dimensions; i++)
	{
		if(i>0)
			result *= type->head.dimensions[i].n_indices;

		for(j=0; j<N_MAX_DIMENSIONS && ind[j][0]!=type->head.dimensions[i].type; j++);

		if(j<N_MAX_DIMENSIONS)
		{
			for(k=0; k<type->head.dimensions[i].n_indices && type->head.dimensions[i].indices[k]!=ind[j][1]; k++);

			if(k<type->head.dimensions[i].n_indices)
				result += k;
			else
			{
				printf("Error finding index %d for dimension #%d=%s in function type_index2disp() .\n", ind[j][1], i, xqcd_type_dim_desc[type->head.dimensions[i].type]);
				return	-1;
			}
		}
		else
		{
			printf("Error finding dimension #%d=%d in function type_index2disp() .\n", i, type->head.dimensions[i].type);
			return	-1;
		}
	}

	return	result;
}

void            type_disp2index(filetype *type, int ind[][2], int disp)
{
	int	i;

	for(i=type->head.n_dimensions-1; i>=0 ; i--)
	{
		ind[i][0] = type->head.dimensions[i].type;
		ind[i][1] = type->head.dimensions[i].indices[disp % type->head.dimensions[i].n_indices];
		disp = (disp - (disp % type->head.dimensions[i].n_indices)) / type->head.dimensions[i].n_indices;
	}

	if(disp>0)
		printf("Warning: disp in function type_disp2index() is too large!\n");

	return;
}

int	type_str2dimtype(char *str)
{
	int	res;

	for(res=0; res<N_MAX_DIMENSION_TYPES; res++)
		if(xqcd_type_dim_desc[res])
			if(strcmp(xqcd_type_dim_desc[res], str)==0)
				break;

	if(res<N_MAX_DIMENSION_TYPES)
		return	res;
	else
	{
		printf("Warning: unrecogonized dimension string : \"%s\"\n", str);
		return	-1;
	}
}

