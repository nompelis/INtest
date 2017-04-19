#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>



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


int code( MPI_Comm *comp )
{
   MPI_Comm comm;
   int nproc,irank;
   int *icon,*jcon,im,jm, nf1,nf2;
   double *x1,*x2;
   int i,j,n;


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
   if( nf1 != 0 ) {
      printf("Could not make face groups (%d) \n", irank);
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
   if( irank == 0 ) {
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


   free(icon);
   free(jcon);
   free(x1);
   free(x2);

   MPI_Comm_free( &comm );
   if( irank == 0 ) printf("Completed \n");

   return(0);
}

