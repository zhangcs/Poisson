/**
 *  This is a program to generate linear system for 3D Possion Problem.
 *
 *     Consider a three-dimensional Possion equation
 *
 * \f[
 *   \frac{du}{dt}-u_{xx}-u_{yy}-u_{zz} = f(x,y,t)\ \ in\ \Omega = (0,1)\times(0,1)\times(0,1)
 * \f]
 * \f[
 *                 u(x,y,z,0) = 0\ \ \ \ \ \ in\ \Omega
 * \f]
 * \f[
 *                        u = 0\ \ \ \ \ \ \ \ \ on\  \partial\Omega
 * \f]
 *
 *  where f(x,y,z,t) = \f$3*\pi^2*u(x,y,z,t) + sin(\pi*x)*sin(\pi*y)*sin(\pi*z)\f$,
 *  and the solution function can be expressed by
 *
 *             \f$u(x,y,z,t) = sin(\pi*x)*sin(\pi*y)*sin(\pi*z)\f$
 *
 *  Created by peghoty 2010/08/04
 *
 *  Xiangtan University
 *  peghoty@163.com
 *
 */

#include "fsls.h"

int
main( int argc, char *argv[])
{
	struct timeval tStart,tEnd;
	int TTest       = 1;
	int arg_index   = 0;
	int print_usage = 0;

	char *MatFile = NULL;
	char *RhsFile = NULL;
	char *SolFile = NULL;
	char filename[120];

	fsls_BandMatrix *A    = NULL;
	fsls_CSRMatrix  *Acsr = NULL;
	fsls_XVector    *b    = NULL;
	fsls_XVector    *u    = NULL;

	int nx,ny,nz,ngrid,nt,i,j;
    double dt = 0.0;

	nx = 10;
	ny = 10;
	nz = 10;
    nt = 0;

	while (arg_index < argc)
	{
		if ( strcmp(argv[arg_index], "-nx") == 0 )
		{
			arg_index ++;
			nx = atoi(argv[arg_index++]);
		}
        if (arg_index >= argc) break;

		if ( strcmp(argv[arg_index], "-ny") == 0 )
		{
			arg_index ++;
			ny = atoi(argv[arg_index++]);
		}
        if (arg_index >= argc) break;

		if ( strcmp(argv[arg_index], "-nz") == 0 )
		{
			arg_index ++;
			nz = atoi(argv[arg_index++]);
		}
        if (arg_index >= argc) break;

		if ( strcmp(argv[arg_index], "-nt") == 0 )
		{
			arg_index ++;
			nt = atoi(argv[arg_index++]);
		}
        if (arg_index >= argc) break;

		if ( strcmp(argv[arg_index], "-help") == 0 )
		{
			print_usage = 1;
			break;
		}
		else
		{
			arg_index ++;
		}
	}

    if (print_usage)
    {
        printf("\n  Usage: %s [<options>]\n\n", argv[0]);
        printf("  -nx <val> : number of interier nodes in x-direction [default: 10]\n");
        printf("  -ny <val> : number of interier nodes in y-direction [default: 10]\n");
        printf("  -nz <val> : number of interier nodes in z-direction [default: 10]\n");
        printf("  -nt <val> : number of interier nodes in t-direction [default:  0]\n");
        printf("  -help     : using help message\n\n");
        exit(1);
    }

	ngrid = nx*ny*nz;
    if (nt != 0) dt = 1./nt;

	printf("\n +++++++++++++ (nx,ny,nz,nt) = (%d,%d,%d,%d)  ngrid = %d +++++++++++\n\n",nx,ny,nz,nt,ngrid);

	MatFile = "./mat_";
	RhsFile = "./rhs_";
	SolFile = "./sol_";

	/*-----------------------------------------------------
	 * construct a linear system
	 *----------------------------------------------------*/
	if (TTest) GetTime(tStart);

	fsls_BuildLinearSystem_7pt3d(nt, nx, ny, nz, &A, &b, &u);

	if (TTest)
	{
		GetTime(tEnd);
		printf("\n >>> total time: %.3f seconds\n\n",mytime(tStart,tEnd));
	}

	sprintf(filename, "%s%dX%dX%d.dat",MatFile,nx,ny,nz);
	fsls_Band2CSRMatrix(A, &Acsr);

    /**
     * the newly added codes aim to prepare for lapack
     * 2011/12/11 by feiteng
     */
	double *A_full = NULL, *B = NULL;
	int LDA, LDB, NRHS = 1;
	int *IPIV, INFO[1];
	fsls_CSR2FullMatrix(Acsr, &A_full);

	if (nt != 0)
	{
		fsls_dtMatrix(dt, ngrid, ngrid, A_full);
	}

	IPIV = ( int *) malloc ( sizeof (int) * ngrid ) ;
	if( IPIV == NULL ) {
		printf("Allocation Error\n");
		exit(1);
	}
	LDA = LDB = ngrid;

	B = ( double *) malloc ( sizeof(double) * ngrid);
	if( B == NULL ){
		printf("Allocation Error\n");
		exit(1);
	}
	memset(B, 0X0, ngrid*sizeof(double));
	dgetrf_(&ngrid, &ngrid, A_full, &LDA, IPIV, INFO);

    /**
     * the newly added codes aim to use lapack routine to solve the
     * linear system, 2011/12/11 by feiteng
     */
	double err;
	double u_;

	for (i = 0;i < nt; ++i)
	{
		err = u_ = 0.;
		for (j = 0;j < ngrid; ++j)
		{
			B[j] += b->data[i*ngrid+j]*dt;
		}
		dgetrs_("N", &ngrid, &NRHS, A_full, &LDA, IPIV, B, &LDB, INFO);

		for (j = 0;j < ngrid; ++j)
		{
			err += (B[j]-u->data[i*ngrid+j])*(B[j]-u->data[i*ngrid+j]);
			u_ += (u->data[i*ngrid+j])*(u->data[i*ngrid+j]);
		}
        printf("time step %3d: rel err = %e\n", i, err/u_);
	}

	if (nt == 0)
	{
		err = u_ = 0.;
		memcpy(B, fsls_XVectorData(b), sizeof(double) * ngrid);
		dgetrs_("N", &ngrid, &NRHS, A_full, &LDA, IPIV, B, &LDB, INFO);
		for (j = 0;j < ngrid; ++j)
		{
			err += (B[j]-u->data[i*ngrid+j])*(B[j]-u->data[i*ngrid+j]);
			u_ += (u->data[i*ngrid+j])*(u->data[i*ngrid+j]);
		}
        printf("rel err = %e\n", err/u_);
	}

	fsls_CSRMatrixPrint(Acsr,filename);

	sprintf(filename, "%s%dX%dX%d.dat",RhsFile,nx,ny,nz);
	fsls_XVectorPrint(b, filename);

	sprintf(filename, "%s%dX%dX%d.dat",SolFile,nx,ny,nz);
	fsls_XVectorPrint(u, filename);

	/*------------------------------------------------------
	 * added for Chensong's SAMG-testing
	 *-----------------------------------------------------*/
	// fsls_WriteSAMGData(Acsr, b, u);

	/*------------------------------------------------------
	 * free some staff
	 *-----------------------------------------------------*/
	fsls_BandMatrixDestroy(A);
	fsls_CSRMatrixDestroy(Acsr);
	fsls_XVectorDestroy(b);
	fsls_XVectorDestroy(u);

	free(A_full);
	free(B);
	free(IPIV);

	return(0);
}

