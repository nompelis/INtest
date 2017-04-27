#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "intest.h"


//
// Function to create two sets of faces, one of quadrilaterals and one of
// triangles, given i- and j-direction sizes for a quadrilateral subdivision.
//
int make_faces( int im, int jm,
                int **icon, double **x1, 
                int **jcon, double **x2 )
{
   size_t isize;
   int i,j,n;


   isize = (size_t) (im*jm*5);
   *icon = malloc(  isize*sizeof(int));
   *jcon = malloc(2*isize*sizeof(int));
   isize = (size_t) ((im+1)*(jm+1)*3);
   *x1 = malloc(isize*sizeof(double));
   *x2 = malloc(isize*sizeof(double));
   if( *icon == NULL || *jcon == NULL || *x1 == NULL || *x2 == NULL ) {
      if( *icon != NULL ) free( *icon );
      if( *jcon != NULL ) free( *jcon );
      if( *x1   != NULL ) free( *x1 );
      if( *x2   != NULL ) free( *x2 );
      return(1);
   }


   for(j=0;j<jm+1;++j) {
   for(i=0;i<im+1;++i) {
      n = j*(im+1) + i;
   // printf("i= %d  j= %d    n= %d \n",i,j,n);//HACK

      (*x1)[3*n+0] = (double) i;
      (*x1)[3*n+1] = (double) j;
      (*x1)[3*n+2] = 0.0;

      (*x2)[3*n+0] = (double) i + 0.2;
      (*x2)[3*n+1] = (double) j + 0.2;
      (*x2)[3*n+2] = 0.1;
   }}

   n = 0;
   for(j=0;j<jm;++j) {
   for(i=0;i<im;++i) {
      (*icon)[5*n+0] = 4;
      (*icon)[5*n+1] = (j  )*(im+1) + i + 1;
      (*icon)[5*n+2] = (j  )*(im+1) + i+1 + 1;
      (*icon)[5*n+3] = (j+1)*(im+1) + i+1 + 1;
      (*icon)[5*n+4] = (j+1)*(im+1) + i + 1;
      ++n;
   }}

   n = 0;
   for(j=0;j<jm;++j) {
   for(i=0;i<im;++i) {
      (*jcon)[5*n+0] = 3;
      (*jcon)[5*n+1] = (j  )*(im+1) + i + 1;
      (*jcon)[5*n+2] = (j  )*(im+1) + i+1 + 1;
      (*jcon)[5*n+3] = (j+1)*(im+1) + i + 1;
      (*jcon)[5*n+4] = (*jcon)[5*n+3];
      ++n;

      (*jcon)[5*n+0] = 3;
      (*jcon)[5*n+1] = (j+1)*(im+1) + i + 1;
      (*jcon)[5*n+2] = (j  )*(im+1) + i+1 + 1;
      (*jcon)[5*n+3] = (j+1)*(im+1) + i+1 + 1;
      (*jcon)[5*n+4] = (*jcon)[5*n+3];
      ++n;
   }}

   return(0);
}


//
// External (library) pointer to function that will be testing for overlaps
// This variable is a global variable in the face-matching library. It is
// a pointer to a fcuntion with appropriate arguments. The C code outside the
// library will have to set this library's global variable ahead of time. In
// this way the library uses an external function to test for overlaps; it is
// not part of the library.
//
extern double (*incg_FaceOverlapFunction)( int n, const double *xyz1,
                                           int m, const double *xyz2 );

//
// A function to be built by the user in order to perform overlap testing
// (This function presently returns a negative number, which implies no overlap)
//
double my_overlap_function( int n, const double *xyz1,
                            int m, const double *xyz2 )
{
   double area = -99999.0;

//// make every element overlap with every one for testing purposes
   area = 1.0;

   return( area );
}


int code( MPI_Comm *comp )
{
   MPI_Comm comm;
   int nproc,irank;
   int *icon,*jcon,im,jm, nf1,nf2;
   double *x1,*x2;
   int i,j,n;
   int *ifaces,*jfaces;
   double *rpoints,*qpoints;
   size_t isize;
   int ierr=0;


   //
   // Establish MPI stuff
   //
   MPI_Comm_dup( *comp, &comm );
   MPI_Comm_size( comm, &nproc );
   MPI_Comm_rank( comm, &irank );


   //
   // Make face groups
   //
   im = 2;
   jm = 3;
   nf1 = make_faces( im, jm, &icon, &x1, &jcon, &x2 );
   MPI_Allreduce( MPI_IN_PLACE, &nf1, 1, MPI_INT, MPI_SUM, comm );
   if( nf1 != 0 ) {
      printf("Could not make face groups (%d) \n", irank);
      if( icon != NULL ) free( icon );
      if( x1   != NULL ) free( x1 );
      if( jcon != NULL ) free( jcon );
      if( x2   != NULL ) free( x2 );
      return(1);
   }
   nf1 = im*jm;
   nf2 = nf1*2;

   //
   // Make a shift to all faces as if data was sequential when globalized
   //
   for(j=0;j<jm+1;++j) {
   for(i=0;i<im+1;++i) {
      int n = j*(im+1) + i;
      x1[n*3+0] += (double) (irank*im);
      x2[n*3+0] += (double) (irank*im);
   }}

   //
   // Dump the structures
   //
#ifdef _DEBUG_
   if( irank == 0 ) {    // select an MPI process manually (hardwired) here
      FILE *fp = fopen("FACES.dat","w");

      fprintf(fp,"variables = x y z \n");
      fprintf(fp,"zone T=\"quads\", i=%d,j=%d, f=point\n", im+1,jm+1);
      for(j=0;j<jm+1;++j) {
      for(i=0;i<im+1;++i) {
         int n = j*(im+1) + i;
         fprintf(fp, "%lf %lf %lf \n",x1[n*3+0], x1[n*3+1], x1[n*3+2] );
      }}

      fprintf(fp,"zone T=\"triangles\", N=%d,E=%d, F=FEPOINT, ET=QUADRILATERAL\n",
              (im+1)*(jm+1),nf2);
      for(j=0;j<jm+1;++j) {
      for(i=0;i<im+1;++i) {
         n = j*(im+1) + i;
         fprintf(fp, "%lf %lf %lf \n",x2[n*3+0], x2[n*3+1], x2[n*3+2] );
      }}
      n = 0;
      for(j=0;j<jm;++j) {
      for(i=0;i<im;++i) {
         fprintf(fp, "%d %d %d %d \n",
                 jcon[n*5+1], jcon[n*5+2], jcon[n*5+3], jcon[n*5+4] );
         ++n;
         fprintf(fp, "%d %d %d %d \n",
                 jcon[n*5+1], jcon[n*5+2], jcon[n*5+3], jcon[n*5+4] );
         ++n;
      }}

      fclose(fp);
   }
   MPI_Barrier( comm );
#endif


   //
   // Prepare arrays to pass to the face-matching library
   //
   isize = (size_t) (nf1*(5+3));
   ifaces = (int *) malloc(isize*sizeof(int));
   isize = (size_t) (nf1*4*3);
   rpoints = (double *) malloc(isize*sizeof(double));

   isize = (size_t) (nf2*(5+3));
   jfaces = (int *) malloc(isize*sizeof(int));
   isize = (size_t) (nf2*4*3);
   qpoints = (double *) malloc(isize*sizeof(double));

   if( ifaces == NULL || jfaces == NULL ||
       rpoints == NULL || qpoints == NULL ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      printf("Could not allocate memory 2 \n");
      if( ifaces != NULL ) free( ifaces );
      if( jfaces != NULL ) free( jfaces );
      if( rpoints != NULL ) free( rpoints );
      if( qpoints != NULL ) free( qpoints );
      free(icon);
      free(jcon);
      free(x1);
      free(x2);
      return(2);
   }

   // copy element connectivity (surface 1)
   for(i=0;i<nf1;++i) {
      n = i*5;

      ifaces[n + 0] = icon[n + 0];
#ifdef _DEBUG2_
  printf("Should be 4: %d \n", icon[n + 0] );
#endif
      ifaces[n + 1] = i;
      ifaces[n + 2] = -1;
      ifaces[n + 3] = -1;
      ifaces[n + 4] = -1;
   }

   // copy vertex coefficients (surface 1)
   for(i=0;i<nf1;++i) {
      n = i*5;
/////// we need to make a correction here; do not assume 4th index is good
      for(j=0;j<4;++j) {
         rpoints[ (i*4 + j)*3 + 0 ] = x1[ (icon[n+1+j]-1)*3 + 0 ];
         rpoints[ (i*4 + j)*3 + 1 ] = x1[ (icon[n+1+j]-1)*3 + 1 ];
         rpoints[ (i*4 + j)*3 + 2 ] = x1[ (icon[n+1+j]-1)*3 + 2 ];
      }
   }

   // copy element connectivity (surface 2)
   for(i=0;i<nf2;++i) {
      n = i*5;

      jfaces[n + 0] = jcon[n + 0];
#ifdef _DEBUG2_
  printf("Should be 3: %d \n", jcon[n + 0] );
#endif
      jfaces[n + 1] = i;
      jfaces[n + 2] = -1;
      jfaces[n + 3] = -1;
      jfaces[n + 4] = -1;
   }

   // copy vertex coefficients (surface 2)
   for(i=0;i<nf2;++i) {
      n = i*5;
/////// we need to make a correction here; do not assume 4th index is good
      for(j=0;j<4;++j) {
         qpoints[ (i*4 + j)*3 + 0 ] = x2[ (jcon[n+1+j]-1)*3 + 0 ];
         qpoints[ (i*4 + j)*3 + 1 ] = x2[ (jcon[n+1+j]-1)*3 + 1 ];
         qpoints[ (i*4 + j)*3 + 2 ] = x2[ (jcon[n+1+j]-1)*3 + 2 ];
      }
   }

#ifdef _DEBUG2_
if(irank==0){
char filename[20];
sprintf( filename, "TEST.dat");
FILE *fp = fopen( filename, "w" );
fprintf(fp,"variables = x y z\n");
fprintf(fp,"zone T=\"test\", N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",
    4*nf2, nf2 );
for(i=0;i<nf2*4;++i) {
for(j=0;j<3;++j) {
fprintf(fp, " %lf ",qpoints[i*3+j] );
} fprintf(fp, "\n"); }
for(i=0;i<nf2;++i) {
for(j=0;j<4;++j) {
fprintf(fp, " %d ",i*4+j+1 );
} fprintf(fp, "\n"); }
fclose(fp);
}
#endif


   //
   // Assign the global pointer of the library to the user function here
   //
   incg_FaceOverlapFunction = &my_overlap_function;

   //
   // Use the library to resolve dependencies as a test
   //
//while(1) /// LOOP to check for memory leaks
// ierr = incg_PerformFacematch( &comm, nf1, ifaces, rpoints,
//                                      nf2, jfaces, qpoints );

   //
   // We will initialize a face-matching object by requesting that a handle is
   // created for us. The handle will be assocaited internally with an object.
   //
   int ihandle;
//while(1) {
   ierr = incg_Facematch_Init( &ihandle, &comm, nf1, ifaces, rpoints,
                                                nf2, jfaces, qpoints, NULL );

   int num_recv, num_send;
   ierr = incg_Facematch_GetSizes( &ihandle, &num_recv, &num_send );

   int *irdis,*ircnt,*isdis,*iscnt;
   int *irecv,*isend;
   double *recv_area;

   irecv = (int *) malloc(num_recv*sizeof(int));
   isend = (int *) malloc(num_send*sizeof(int));
   recv_area = (double *) malloc(num_send*sizeof(double));
   irdis = (int *) malloc(nproc*sizeof(int));
   ircnt = (int *) malloc(nproc*sizeof(int));
   isdis = (int *) malloc(nproc*sizeof(int));
   iscnt = (int *) malloc(nproc*sizeof(int));
   if( irecv == NULL || 
       isend == NULL || 
       recv_area == NULL || 
       irdis == NULL || 
       ircnt == NULL || 
       isdis == NULL || 
       iscnt == NULL ) { ierr = 1; }
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      free( isend );
      free( irecv );
      free( recv_area );
      free( irdis );
      free( ircnt );
      free( isdis );
      free( iscnt );
   }

   ierr = incg_Facematch_FillArrays( &ihandle,
                              isend, irecv, recv_area,
                              irdis, ircnt, isdis, iscnt );


   ierr = incg_Facematch_Term( &ihandle );
//}



   //
   // Clean up surface structures
   //
   free( ifaces );
   free( jfaces );
   free( rpoints );
   free( qpoints );
   free(icon);
   free(jcon);
   free(x1);
   free(x2);

   //
   // Drop the communicator used in this process
   //
   MPI_Comm_free( &comm );
   if( irank == 0 ) printf("Completed \n");

   return(0);
}

