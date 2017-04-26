#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "intest.h"

#ifdef __cplusplus
extern "C" {


//
// Constructor of face-matcher object
//
incg_FaceMatcher::incg_FaceMatcher( MPI_Comm *comp )
{
   MPI_Comm_dup( *comp, &comm );
   MPI_Comm_size( comm, &nproc );
   MPI_Comm_rank( comm, &irank );

   nel1 = 0;
   icon = NULL;
   xi = NULL;
   nel2 = 0;
   jcon = NULL;
   xj = NULL;

   iacc = NULL;
   idis = NULL;
   icnt = NULL;

   nbuf_size = 0;
   ibuf = NULL;
   rbuf = NULL;
   idata = NULL;
   rdata = NULL;

   idst = NULL;
   jdst = NULL;

#ifdef _DEBUG_
   if( irank == 0 ) printf("incg_FaceMatcher object instantiated\n");
#endif
}

//
// Deconstructor of face-matcher object
//
incg_FaceMatcher::~incg_FaceMatcher( )
{
   if( idis != NULL ) free( idis );
   if( icnt != NULL ) free( icnt );
   if( ibuf != NULL ) free( ibuf );
   if( rbuf != NULL ) free( rbuf );
   if( idata != NULL ) free( idata );
   if( rdata != NULL ) free( rdata );

   if( idst != NULL ) free( idst );
   if( jdst != NULL ) free( jdst );

#ifdef _DEBUG_
   if( irank == 0 ) printf("incg_FaceMatcher object deconstructed \n");
#endif
   MPI_Comm_free( &comm );
}

//
// Method to set the internal variables
//
void incg_FaceMatcher::setData( int icnt_, int *icon_, double *xi_,
                                int jcnt_, int *jcon_, double *xj_ )
{
   nel1 = icnt_;
   icon = icon_;
   xi = xi_;

   nel2 = jcnt_;
   jcon = jcon_;
   xj = xj_;
}

//
// Method to enable acceleration of the matching operations by providing
// (pre-calculated) guiding arrays for performing the matching. The input array
// is meant to be sized [nproc] and contain either 0 or 1, with 1 signifying
// that the A surface elements should be checked against the specific process's
// B surface elements; it implies that some bounding box overlap or similar
// operation has been performed externally.
//
void incg_FaceMatcher::setAccel( int *isb_ )
{
   if( isb_ != NULL ) {
#ifdef _DEBUG_
      if( irank == 0 ) printf("Search acceleration array provided \n");
#endif
      iacc = isb_;
   }
}


//
// Method to prepare the object
//
int incg_FaceMatcher::prepare( void )
{
   size_t isize;
   int n;
   int ierr = 0;


   // sanity checks
   if( nel1 <= 0 || nel2 <= 0 ||
       icon == NULL || jcon == NULL ||
       xi == NULL || xj == NULL ) {
      if( irank == 0 ) printf("The object data have not been set\n");
      return(1);
   }

   // surface element distribution arrays
   isize = (size_t) nproc;
   idis = (int *) malloc(isize*sizeof(int));
   icnt = (int *) malloc(isize*sizeof(int));
   isize += 1;
   idst = (int *) malloc(isize*sizeof(int));
   jdst = (int *) malloc(isize*sizeof(int));
   if( idst == NULL || jdst == NULL || idis == NULL || icnt == NULL ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( idst != NULL ) free( idst ); idst = NULL;
      if( jdst != NULL ) free( jdst ); jdst = NULL;
      if( idis != NULL ) free( idis ); idis = NULL;
      if( icnt != NULL ) free( icnt ); icnt = NULL;
      return(-1);
   }

   for(n=0;n<nproc+1;++n) idst[n] = 0;
   idst[irank+1] = nel1;
   MPI_Allreduce( MPI_IN_PLACE, idst, nproc+1, MPI_INT, MPI_SUM, comm );

   for(n=0;n<nproc+1;++n) jdst[n] = 0;
   jdst[irank+1] = nel2;
   MPI_Allreduce( MPI_IN_PLACE, jdst, nproc+1, MPI_INT, MPI_SUM, comm );

   nbuf_size = jdst[0];       // buffer size for surface2 broadcasts
   for(n=1;n<nproc+1;++n) {
#ifdef _DEBUG2_
if(irank==0) printf("idst[%d]= %d  ",n,idst[n]);
#endif
      idst[n] += idst[n-1];
#ifdef _DEBUG2_
if(irank==0) printf("idst[%d]= %d\n",n,idst[n]);
#endif

      if( jdst[n] > nbuf_size ) nbuf_size = jdst[n];   // buffer size
      jdst[n] += jdst[n-1];
   }

   for(n=0;n<nproc;++n) {
      icnt[n] = jdst[n+1] - jdst[n];       // original count
      if( iacc != NULL ) {
         if( iacc[n] == 0 ) icnt[n] = 0;   // count is zero
      }
   }
   idis[0] = 0;
   for(n=1;n<nproc;++n) idis[n] = idis[n-1] + icnt[n-1];
#ifdef _DEBUG2_
if(irank==0) for(n=0;n<nproc;++n)
   printf("idis[%d]= %d   icnt[%d]= %d \n",n,idis[n], n,icnt[n]);
#endif


   // buffers for surface2 broadcasts and...
   isize = (size_t) (nbuf_size);
   ibuf = (int *)    malloc( isize*5   * sizeof(int) );
   rbuf = (double *) malloc( isize*4*3 * sizeof(double) );
   // ...buffers for entire foreign surface2 element set
   isize = 0;
   for(n=0;n<nproc;++n) isize += (size_t) (icnt[n]);
//printf("IRANK= %d   SIZE= %ld \n",irank, (long) isize );//HACK
   idata = (int *)    malloc( isize*5   * sizeof(int) );
   rdata = (double *) malloc( isize*4*3 * sizeof(double) );
   if( ibuf == NULL || rbuf == NULL || idata == NULL || rdata == NULL ) ierr= 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( ibuf != NULL ) free( ibuf ); ibuf = NULL;
      if( rbuf != NULL ) free( rbuf ); rbuf = NULL;
      if( idata != NULL ) free( idata ); idata = NULL;
      if( rdata != NULL ) free( rdata ); rdata = NULL;

      free( idst ); idst = NULL;
      free( jdst ); jdst = NULL;
      free( idis ); idis = NULL;
      free( icnt ); icnt = NULL;
      return(-2);
   }

#ifdef _DEBUG2_
n = irank;
if(irank==n){    // CHANGE THE PROCESS NUMBER TO TEST OTHER PROCESSES
char filename[20];
sprintf( filename, "TEST_%.5d.dat", n );
FILE *fp = fopen( filename, "w" );
int knt = nel2;
fprintf(fp,"variables = x y z\n");
fprintf(fp,"zone T=\"test_%d\", N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",
    n, 4*knt, knt );
for(int i=0;i<knt*4;++i) {
for(int k=0;k<3;++k) {
fprintf(fp, " %lf ",xj[i*3+k] );
} fprintf(fp, "\n"); }
for(int i=0;i<knt;++i) {
for(int j=0;j<4;++j) {
fprintf(fp, " %d ",i*4+j+1 );
} fprintf(fp, "\n"); }
fclose(fp);
}
#endif

   // linear broadcast of surface2 elements
   for(n=0;n<nproc;++n) {
      int *idum;
      double *rdum;
      int knt = jdst[n+1] - jdst[n];     // count of surface2 elements of "n"

      if( icnt[n] > 0 ) {                // have overlap with this process
         idum = &( idata[ 5*idis[n] ] ); // point straight into main arrays
         rdum = &( rdata[ 4*3*idis[n] ] );
      } else {
         idum = ibuf;                    // point into the buffers
         rdum = rbuf;
      }

      if( irank == n ) {                 // this is the broadcasting process
         for(int i=0;i<nel2;++i) {       // sender fills buffers from origin
            for(int j=0;j<4;++j) {       // always four nodes (tri. or quad.)
               int ii = i*4 + j;         // index of node in buffer
               rdum[ ii*3 + 0 ] = xj[ ii*3 + 0 ];
               rdum[ ii*3 + 1 ] = xj[ ii*3 + 1 ];
               rdum[ ii*3 + 2 ] = xj[ ii*3 + 2 ];
            }

            idum[ i*5 + 0 ] = i*1000 + 111;  // HACK
            idum[ i*5 + 1 ] = i*1000 + 222;  // HACK
            idum[ i*5 + 2 ] = i*1000 + 333;  // HACK
            idum[ i*5 + 3 ] = i*1000 + 444;  // HACK
            idum[ i*5 + 4 ] = i*1000 + 555;  // HACK
         }
      }

      MPI_Bcast( idum, 5*knt, MPI_INT, n, comm );
      MPI_Bcast( rdum, 4*3*knt, MPI_DOUBLE, n, comm );

      if( irank == 0 ) printf("Broadcast completed for process %d \n", n);
#ifdef _DEBUG_
if(irank==0){    // CHANGE THE PROCESS NUMBER TO TEST OTHER PROCESSES
char filename[20];
sprintf( filename, "TEST_%.5d.dat", irank );
FILE *fp = fopen( filename, "w" );
fprintf(fp,"variables = x y z\n");
fprintf(fp,"zone T=\"test_%d\", N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",
    n, 4*knt, knt );
for(int i=0;i<knt*4;++i) {
for(int k=0;k<3;++k) {
fprintf(fp, " %lf ",rdum[i*3+k] );
} fprintf(fp, "\n"); }
for(int i=0;i<knt;++i) {
for(int j=0;j<4;++j) {
fprintf(fp, " %d ",i*4+j+1 );
} fprintf(fp, "\n"); }
fclose(fp);
}
#endif
   }

   // drop the buffers used for transfers
   free( ibuf ); ibuf = NULL;
   free( rbuf ); rbuf = NULL;
   nbuf_size  = 0;

   // prepare the lists of overlaps for each owned element of surface 1
   for(n=0;n<nel1;++n) {
      std::list< overlap_s > ltmp;

      lists.push_back( ltmp );

#ifdef _DEBUG_
      if( irank == 0 ) {
         printf("Overlap list for process %d at %p, size %ld \n",
                irank, &( lists[n] ), lists[n].size() );
      }
#endif
   }

   return(0);
}

//
// Method to perform the face-matching
//
int incg_FaceMatcher::perform( void )
{
   int n;
   double xyz1[4*3], xyz2[4*3];
   int ierr = 0;





#ifdef _DEBUG_
   MPI_Barrier( comm );
   if( irank == 0 ) {
      for(n=0;n<nel1;++n)
         printf("Overlap list after filling for process %d at %p, size %ld \n",
                irank, &( lists[n] ), lists[n].size() );
   }
#endif

#ifdef _DEBUG_
if(irank==irank){    // CHANGE THE PROCESS NUMBER TO TEST OTHER PROCESSES
char filename[20];
sprintf( filename, "LIST_%.5d.dat", irank );
FILE *fp = fopen( filename, "w" );
   fprintf( fp, "###### Process %d element overlaps ###### \n", irank );

   for(n=0;n<nel1;++n) {
      fprintf( fp, "### Element %d ### \n", n );

      std::list< overlap_s > *lp = &( lists[n] );
      fprintf( fp, " %ld \n", lp->size() );

      std::list< overlap_s > :: iterator il;
      for( il = lp->begin(); il != lp->end(); ++ il ) {
      // struct overlap_s *odum = (struct overlap_s *) il;

         fprintf( fp, " %d %d %lf \n", (*il).iproc, (*il).ielem, (*il).a );
      }
   }
fclose(fp);
}
#endif

   return(0);
}



//
// API function to invoke the functionality from the C side
//
int incg_PerformFacematch( MPI_Comm *comm, int icnt1, int *ilist, double *x1,
                                           int icnt2, int *jlist, double *x2 )
{
   int nproc,irank;
   MPI_Comm_size( *comm, &nproc );
   MPI_Comm_rank( *comm, &irank );
   printf("The function was invoked by MPI process %d \n", irank);
   int ierr=0;


   // constructing the object
   incg_FaceMatcher fm( comm );

   // provide variables to the object
   fm.setData( icnt1, ilist, x1,  icnt2, jlist, x2 );

   // call this if we are using some pre-search-based acceleration technique
   (void) fm.setAccel( NULL );
   // (this is a test)
int ijunk[] = {
 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
   (void) fm.setAccel( ijunk );

   // method to create internal structures
   ierr = fm.prepare();
   if( ierr != 0 ) return(1);

   // method to perform face-matching
   ierr = fm.perform();
   if( ierr != 0 ) return(1);




   return(0);
}


}
#endif

