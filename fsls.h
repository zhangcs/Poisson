#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#define GetTime(a) gettimeofday(&a,NULL)
#define mytime(a,b) ((b.tv_sec-a.tv_sec) + (float)(b.tv_usec-a.tv_usec)/1000000.0)

#define PI 3.1415926535897932
#define fsls_max(a,b)  (((a)<(b)) ? (b) : (a))
#define fsls_min(a,b)  (((a)<(b)) ? (a) : (b))
#define fsls_TFree(ptr) ( fsls_Free((char *)ptr), ptr = NULL )
#define fsls_CTAlloc(type, count) ( (type *)fsls_CAlloc((size_t)(count), (size_t)sizeof(type)) )
void  fsls_Free( char *ptr );
char *fsls_CAlloc( size_t count, size_t elt_size );
int   fsls_OutOfMemory( size_t size );

typedef struct
{
   double  *data;
   int     *i;
   int     *j;
   int      num_rows;
   int      num_cols;
   int      num_nonzeros;
   int     *rownnz;
   int      num_rownnz; 
   int      owns_data;
} fsls_CSRMatrix;

#define fsls_CSRMatrixData(matrix)         ((matrix) -> data)
#define fsls_CSRMatrixI(matrix)            ((matrix) -> i)
#define fsls_CSRMatrixJ(matrix)            ((matrix) -> j)
#define fsls_CSRMatrixNumRows(matrix)      ((matrix) -> num_rows)
#define fsls_CSRMatrixNumCols(matrix)      ((matrix) -> num_cols)
#define fsls_CSRMatrixNumNonzeros(matrix)  ((matrix) -> num_nonzeros)
#define fsls_CSRMatrixRownnz(matrix)       ((matrix) -> rownnz)
#define fsls_CSRMatrixNumRownnz(matrix)    ((matrix) -> num_rownnz)
#define fsls_CSRMatrixOwnsData(matrix)     ((matrix) -> owns_data)

typedef struct
{
   double  *data;
   int      size;
   int      owns_data;
   int      num_vectors;
   int      multivec_storage_method;
   int      vecstride, idxstride;

} fsls_Vector;

#define fsls_VectorData(vector)                  ((vector) -> data)
#define fsls_VectorSize(vector)                  ((vector) -> size)
#define fsls_VectorOwnsData(vector)              ((vector) -> owns_data)
#define fsls_VectorNumVectors(vector)            ((vector) -> num_vectors)
#define fsls_VectorMultiVecStorageMethod(vector) ((vector) -> multivec_storage_method)
#define fsls_VectorVectorStride(vector)          ((vector) -> vecstride )
#define fsls_VectorIndexStride(vector)           ((vector) -> idxstride )  

typedef struct
{
   int      n;         // order of the matrix
   int      nx;        // number of nodes along x-direction(excluding boundary nodes)
   int      ny;	       // number of nodes along y-direction(excluding boundary nodes)
   int      nz;        // number of nodes along z-direction(excluding boundary nodes)
   int      nband;     // the number of offdiagonal bands
   
   int     *offsets;   // offsets of the offdiagonal bands (length is nband), offsets are 
                       // ordered in the ascendling manner, the negative and positive values
                       // corresband to lower left bands and upper right bands, respectively  
                     
   double  *diag;      // diagonal entries (length is n)
   double **offdiag;   // off-diagonal entries (dimension is nband X n), 
                       // offdiag[i][j],i=0(1)nband-1,j=0(1)n-1: the j-th entry on the i-th offdiagonal band.
   double  *data_ext;  // data part, including diag_ext and offdiag_ext                    
   
} fsls_BandMatrix;

#define fsls_BandMatrixN(matrix)         ((matrix) -> n)
#define fsls_BandMatrixNx(matrix)        ((matrix) -> nx)
#define fsls_BandMatrixNy(matrix)        ((matrix) -> ny)
#define fsls_BandMatrixNz(matrix)        ((matrix) -> nz)
#define fsls_BandMatrixNband(matrix)     ((matrix) -> nband)
#define fsls_BandMatrixOffsets(matrix)   ((matrix) -> offsets)
#define fsls_BandMatrixDiag(matrix)      ((matrix) -> diag)
#define fsls_BandMatrixOffdiag(matrix)   ((matrix) -> offdiag)
#define fsls_BandMatrixDataExt(matrix)   ((matrix) -> data_ext)

typedef struct
{
   int      size;     // length of the vector	                       
   double  *data;     // data of the vector (length is size)
   double  *data_ext; // data part, including extended data
   
} fsls_XVector;

#define fsls_XVectorSize(vector)     ((vector) -> size)
#define fsls_XVectorData(vector)     ((vector) -> data)
#define fsls_XVectorDataExt(vector)  ((vector) -> data_ext)

int fsls_BandMatrixPrint( fsls_BandMatrix *A, char *file_name );
fsls_BandMatrix *fsls_BandMatrixRead( char *file_name );
void fsls_BandMatrixDestroy( fsls_BandMatrix *matrix );
fsls_BandMatrix *fsls_BandMatrixCreate( int n, int nband );
void fsls_BandMatrixInitialize( fsls_BandMatrix *matrix );
void fsls_TriBand2FullMatrix( fsls_BandMatrix *A, double **full_ptr );
int fsls_Band2FullMatrix( fsls_BandMatrix *A, double **full_ptr );
int fsls_CheckDiagOdd( fsls_BandMatrix *matrix );
int fsls_XVectorPrint( fsls_XVector *vector, char *file_name );
fsls_XVector *fsls_XVectorRead( char *file_name );
fsls_XVector *fsls_XVectorCreate( int size );
int fsls_XVectorInitialize( fsls_XVector *vector );
int fsls_XVectorCopy( fsls_XVector *x, fsls_XVector *y );
int fsls_XVectorSetConstantValues( fsls_XVector *vector, double value );
int fsls_XVectorDestroy( fsls_XVector *vector );

void 
fsls_BuildLinearSystem_5pt2d( int               nx, 
                              int               ny,
                              fsls_BandMatrix **A_ptr, 
                              fsls_XVector    **f_ptr,
                              fsls_XVector    **u_ptr );
void 
fsls_BuildLinearSystem_7pt3d( int               nx, 
                              int               ny,
                              int               nz,
                              fsls_BandMatrix **A_ptr, 
                              fsls_XVector    **f_ptr,
                              fsls_XVector    **u_ptr ); 



int fsls_Band2CSRMatrix( fsls_BandMatrix *B, fsls_CSRMatrix **A_ptr );
int fsls_CSRMatrixPrint( fsls_CSRMatrix *matrix, char *file_name );
fsls_CSRMatrix *fsls_CSRMatrixCreate( int num_rows,int num_cols,int num_nonzeros );
int fsls_CSRMatrixInitialize( fsls_CSRMatrix *matrix ); 
int fsls_CSRMatrixDestroy( fsls_CSRMatrix *matrix );  
fsls_CSRMatrix *fsls_CSRMatrixDeleteZeros( fsls_CSRMatrix *A, double tol );  
int fsls_WriteSAMGData( fsls_CSRMatrix *A, fsls_XVector *b, fsls_XVector *u ); // newly added 2010/08/23


