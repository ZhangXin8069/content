#ifndef IO_GENERAL_H
#define IO_GENERAL_H

#define	N_MAX_HEADLENGTH 102400
#define	N_MAX_DIMENSIONS 16
#define	N_MAX_INDICES 1024
#define	MAX_N_FILENAMELENGTH 256
#define	N_MAXERRORBUF 1024

#define	N_MAX_DIMENSION_TYPES 1024

#define	IO_SAVE_SERIAL 

#include	<stdio.h>

#ifdef MPIIO
#include	<mpi.h>
#include	"vector-util.h"
#endif

typedef enum{
	dim_other=0,
	dim_x, dim_y, dim_z, dim_t,
	dim_d, dim_c, dim_d2, dim_c2,
	dim_complex,
	dim_mass,
	dim_smear,
	dim_displacement,

	dim_s_01, dim_s_02, dim_s_03, dim_s_11, dim_s_12, dim_s_13,
	dim_d_01, dim_d_02, dim_d_03, dim_d_11, dim_d_12, dim_d_13,

	dim_conf,
	dim_operator,
	dim_momentum,
	dim_direction,
	dim_t2,
	dim_mass2,

	dim_column,
	dim_row,
	dim_temporary,
	dim_temporary2,
	dim_temporary3,
	dim_temporary4,

	dim_errorbar,			// 0 means average, 1 means errorbar, ...

	dim_operator2,

	dim_param,
	dim_fitleft,
	dim_fitright,

	dim_jackknife,
	dim_jackknife2,
	dim_jackknife3,
	dim_jackknife4,

	dim_summary,			// 0 means average, 1 means standard deviation, 2 means minimal value, 3 means maximum value, 4 means standard error, 
					// 5 means median, ...

	dim_channel,
	dim_channel2,

	dim_eigen,

	dim_d_row,			// on matrix multiplication, row is contracted with the left operand, col is contracted with the right operand.
	dim_d_col,
	dim_c_row,
	dim_c_col,

	dim_parity,			// dimension for different parities. we use 1/-1 for +/- parities for baryons.

	dim_noise,
	dim_evenodd,

	dim_disp_x,
	dim_disp_y,
	dim_disp_z,
	dim_disp_t,

	dim_t3,
	dim_t4,
	dim_t_source,
	dim_t_current,
	dim_t_sink,

	dim_nothing,			// do not use this unless for unused data.

        dim_bootstrap,

	// add new dimensions here and add a string name in xqcd_type_dim_desc[] in io_general.c
	// ...

	dim_last
}dimensiontype;

extern	char    *xqcd_type_dim_desc[N_MAX_DIMENSION_TYPES];

typedef struct{
          int	type;
          int	n_indices;
          int	indices[N_MAX_INDICES];
}one_dim;

typedef union{

	struct{

		int	n_dimensions;
		one_dim dimensions[N_MAX_DIMENSIONS];

	}head;

	char	blank[N_MAX_HEADLENGTH];

}filetype;
// The data in files is saved as array[dimensions[0].n_indices][dimensions[1].n_indices][dimensions[2].n_indices]...


typedef struct{

	char		name[MAX_N_FILENAMELENGTH];
	filetype	type;


#ifdef MPIIO

	MPI_File	fp_mpi;
	MPI_Datatype	doubletype;
	MPI_Datatype	basetype;
	MPI_Datatype	memdatatype;
	MPI_Datatype	filedatatype;

	int		rundimensions;
//	int		blocksize;
//	double		*buf;
//	int		*displacements;
//	int		*blocksizes;
//	int		rundimensions, blocksize, nblocks, maxblocks;

#else

	FILE		*fp;

#endif

}filehandlebase;

typedef	filehandlebase*	filehandle;


typedef	struct{
	int	x0, y0, z0, t0;
	int	dx, dy, dz, dt;
	int	eo;			// odd : 1, even : 2, all : 3
}gridtype;


#ifdef  __cplusplus
extern "C" {
#endif



int		type_str2dimtype(char *str);						// return the index of a element of dimensiontype.

void		xqcd_type_print(filetype *type);
int		xqcd_file_print(char *name, char *tag);

int		type_index2disp(filetype *type, int ind[][2]);				// return -1 if any error occurs.
void		type_disp2index(filetype *type, int ind[][2], int disp);		// the indices are stored in ind WITHOUT array boundary check. (N_MAX_DIMENSIONS*2*sizeof(int) for ind is safe)
int		xqcd_file_write_once(char *name, filetype *type, double *pdata);
double*		xqcd_file_read_once(char *name, filetype *type);			// the "type" will be fed in the head of the file.

filehandle	xqcd_file_open(char *name, filetype *type);				// set type=NULL to read. set type!=NULL to write.
int		xqcd_file_close(filehandle file);
int		xqcd_file_write(filehandle file, filetype *dataset, double *pdata);	// set dataset=NULL to write all data.
double*		xqcd_file_read(filehandle file, filetype *dataset);			// set dataset=NULL to read all data.


#ifdef MPIIO

filehandle	xqcd_file_open_write_mpi_milc(char *name, filetype *type, gridtype *pgrid);	// The xyzt dimensions should be together with each other.
int		xqcd_file_write_mpi_milc_begin(filehandle file, double *vec, int *rundims);
int		xqcd_file_write_mpi_milc_end(filehandle file, double *vec);
void		xqcd_file_close_mpi_milc(filehandle file);

void		make_filetype_propagators(filetype *type, gridtype *pgrid);
void		make_filetype_loop_x(filetype *type, gridtype *pgrid);
void		make_filetype_loop_x2(filetype *type, gridtype *pgrid);

/*
filehandle	xqcd_file_open_write_mpi(char *name, filetype *type, int fixdimensions, int maxblocks);
void		xqcd_file_close_mpi(filehandle file);

int		xqcd_file_write_buf_mpi(filehandle file, int *dimensions, double *pdata);
int		xqcd_file_write_commit_mpi(filehandle file);
int		xqcd_file_reduce_buf_mpi(filehandle file);
int		xqcd_file_write_index_node0(filehandle file);

//int		xqcd_file_read_mpi(filehandle file, filetype *dataset, double *pdata);	// set dataset=NULL to read all data.
*/

gridtype*	set_grid(gridtype *pgrid, int icol);	// icol==-1 means all the columns with all time-slices. icol==-2 means all the columns with time grid.
int		is_on_grid(gridtype *pgrid, int x, int y, int z, int t);	// return -1 for true and 0 for false.
int		n_site_on_grid(gridtype *pgrid, int timeslice);	// return the number sites on the grid. timeslice==-1 means to count all timeslices.
#define	FORGRID(grid, x, y, z, t)	for((t)=0;(t)<nt;(t)++) for((z)=0;(z)<nz;(z)++) for((y)=0;(y)<ny;(y)++) for((x)=0;(x)<nx;(x)++) if(is_on_grid(&(grid),(x),(y),(z),(t)))

/*
filehandle      xqcd_file_open_write_mpi_noise(char *name, gridtype *pgrid);
int             xqcd_file_write_buf_mpi_noise(filehandle file, gridtype *pgrid, vector* data, int icol);

filehandle	xqcd_file_open_write_mpi_propagator(char *name, gridtype *pgrid);
int		xqcd_file_write_buf_mpi_vector2propagator(filehandle file, gridtype *pgrid, int color, int dirac, vector data);

filehandle	xqcd_file_open_write_mpi_propagator_set(char *name, gridtype *pgrid, int dimtype, int nset, int *setind);
int		xqcd_file_write_buf_mpi_vector2propagator_set(filehandle file, gridtype *pgrid, int iset, int color, int dirac, vector data);
*/

#endif


int		get_endian();
void		switch_endian_double(double *pdata, int n);
void		switch_endian_int(int *pdata, int n);
int*		make_endian_switch_buffer_int(int *pdata, int n);
double*		make_endian_switch_buffer_double(double *pdata, int n);



#ifdef  __cplusplus
    }
#endif  /* end of __cplusplus */




#endif
