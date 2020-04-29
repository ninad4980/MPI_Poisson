/*
 ============================================================================
 Name        : MyMpi.c
 Author      : Ninad
 Version     :
 Copyright   : Your copyright notice
 Description : Program Fairness
 ============================================================================
 */

#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <iostream>
#define nx  5 
#define ny  5 
#define nz  3 
#define maxn nx*ny*nz
#define np  3

int main(int argc, char** argv) {

int rank, value, size, errcnt, toterr, i, j, k;
    MPI_Status status;
    double phi[nx][ny][nz];
    double b[nx][ny][nz];

    double phiNew[nx][ny][nz];
    double     diffNorm, gdiffNorm;
int nz_local= nz/np+2;

    double phiLocal[nz_local][ny][nx];
    double rhoLocal[nz_local][ny][nx];
    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

 //   if (size != 3) MPI_Abort( MPI_COMM_WORLD, 1 );

/* Fill the data as specified */


for (k=0; k<nz_local; k++) {    
for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) {    
	phiLocal[k][j][i] = -1;
	rhoLocal[k][j][i] = -100;
    }}}


for (k=1; k<nz_local-1; k++) {    
for (j=1; j<ny-1; j++) {
	for (i=1; i<nx-1; i++) {    
	phiLocal[k][j][i] = rank;
    }}}


int count = 0;

for (k=0; k<nz_local; k++) {    
for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) {    
	phiNew[k][j][i]= phiLocal[k][j][i] ;
    }}}


do{
	
	/** Send to Right **/
	if (rank < np-1 ) 
		MPI_Send(  &phiLocal[nz_local-2][0][0], nx*ny , MPI_DOUBLE, rank + 1, 0,   MPI_COMM_WORLD );
	/** Rec from left **/
	if (rank > 0 )
		MPI_Recv( &phiLocal[0][0][0], nx*ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status );
	/** Send to left **/
    	if (rank > 0) 
		MPI_Send( &phiLocal[1][0][0], nx*ny, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD );
	/** Rec from Right **/
	if (rank < size - 1) 
		MPI_Recv( &phiLocal[nz_local-1][0][0], nx*ny, MPI_DOUBLE, rank + 1, 1,MPI_COMM_WORLD, &status );
	
	count ++;
	diffNorm = 0.0;
	
	for (k=1; k<nz_local-1; k++) 
		for (j=1; j<ny-1; j++) 
			for (i=1; i<nx-1; i++) {    

				phiNew[k][j][i] = rhoLocal[k][j][i]- (phiLocal[k+1][j][i] + phiLocal[k-1][j][i]
						 + phiLocal[k][j+1][i] + phiLocal[k][j-1][i] +
					      phiLocal[k][j][i+1] + phiLocal[k][j][i-1]) / 6.0;
				diffNorm += (phiNew[k][j][i] - phiLocal[k][j][i]) * 
				            (phiNew[k][j][i] - phiLocal[k][j][i]);
				    }

	/* Only transfer the interior points */
	for (k=1; k<nz_local-1; k++) 
		for (j=1; j<ny-1; j++) 
			for (i=1; i<nx-1; i++) {    
				phiLocal[k][j][i] = phiNew[k][j][i];
				}


	MPI_Allreduce( &diffNorm, &gdiffNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	gdiffNorm = sqrt( gdiffNorm );

	if (rank == 0) printf( "At iteration %d, diff is %e\n", count, gdiffNorm );

    } while (gdiffNorm > 1.0e-2 && count < 100);

if (rank ==1){

	for (k=0; k<nz_local; k++)
		{
	for (j=0; j<ny; j++) 
		{
	for (i=0; i<nx; i++) 
		{		
		std :: cout << phiLocal[k][j][i] << "\t" ;	
			}
			std :: cout << "\t" ;	
		for (i=0; i<nx; i++) 
		{		
		std :: cout << rhoLocal[k][j][i] << "\t" ;	
		}

	std :: cout << "\n";
	}
	std :: cout << "\n";
	}
}


MPI_Finalize( );
    return 0;
}

