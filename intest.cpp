#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <indsde.h>

#include "intest.h"
//#include <unistd.h> // HACK remove after done with sleep()

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

   overlap_func = NULL;

   isdis = NULL;
   iscnt = NULL;

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

   if( isdis != NULL ) free( isdis );
   if( iscnt != NULL ) free( iscnt );

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
fprintf(fp,"### surface 1 ### \n");
fprintf(fp,"variables = x y z\n");
fprintf(fp,"zone T=\"own_%d\", N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",
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

            idum[ i*5 + 0 ] = jcon[i*5 + 0]; // number of nodes of element
            idum[ i*5 + 1 ] = irank;         // owner/sender process
            idum[ i*5 + 2 ] = i;             // index in owner process
            idum[ i*5 + 3 ] = i*1000 + 444;  // HACK
            idum[ i*5 + 4 ] = i*1000 + 555;  // HACK
         }
      }

      MPI_Bcast( idum, 5*knt, MPI_INT, n, comm );
      MPI_Bcast( rdum, 4*3*knt, MPI_DOUBLE, n, comm );

      if( irank == 0 ) printf("Broadcast completed for process %d \n", n);
#ifdef _DEBUG2_
if(irank==0){    // CHANGE THE PROCESS NUMBER TO TEST OTHER PROCESSES
char filename[20];
sprintf( filename, "TEST_%.5d_%.5d.dat", irank, n );
FILE *fp = fopen( filename, "w" );
fprintf(fp,"### incoming surface 2 ### \n");
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
         printf("Overlap list (elem %d) for process %d at %p, size %ld \n",
                n, irank, &( lists[n] ), lists[n].size() );
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

   //
   // sanity check
   //
   if( overlap_func == NULL ) {
      if( irank == 0 ) printf("No overlap-test function pointer provided \n");
      return(1);
   }

   //
   // Loop over all surface 1 elems and match against accumulated per-proc elems
   //
   for(int i=0;i<nel1;++i) {
      int ic1 = icon[i*5 + 0];

      for(int mm=0;mm<4;++mm) {
         for(int nn=0;nn<3;++nn) {
            xyz1[mm*3 + nn] = xi[i*4*3 + mm*3 + nn];
      }}
#ifdef _DEBUG2_
  if(irank==0) {
     printf("ELEMENT: %d \n", i );
     for(int mm=0;mm<ic1;++mm) {
        printf(" %lf %lf %lf \n",
           xyz1[mm*3 + 0],
           xyz1[mm*3 + 2],
           xyz1[mm*3 + 1] );
     }
  }
#endif

      for(n=0;n<nproc;++n) {
         for(int j=idis[n];j<idis[n]+icnt[n];++j) {
//if(irank==0) printf("QUERYING: %d   %d,%d \n", i, n, j );

            int ic2 = idata[i*5 + 0];
            for(int mm=0;mm<4;++mm) {
               for(int nn=0;nn<3;++nn) {
                  xyz2[mm*3 + nn] = rdata[j*4*3 + mm*3 + nn];
            }}

#ifdef _DEBUG2_
  if(irank==0) {
     printf("OTHER ELEMENT: (proc %d) %d \n", n,j );
     for(int mm=0;mm<ic2;++mm) {
        printf(" %lf %lf %lf \n",
           xyz2[mm*3 + 0],
           xyz2[mm*3 + 2],
           xyz2[mm*3 + 1] );
     }
  }
//sleep(1);
#endif

            // invoke an overlap function
            double area = (*overlap_func)( ic1, xyz1, ic2, xyz2 );
         // printf("AREA= %lf \n", area );

            if( area > 0.0 ) {
               // perform sanity check
               if( n != idata[j*5 + 1] ) { ierr = 1; }

               struct overlap_s otmp;
               // assign structure members
               otmp.iproc = n;
               otmp.ielem = idata[j*5 + 2];
               otmp.a = area;
               // add overlap to element
               lists[i].push_back( otmp );
            }

         }
      }
   }

   //
   // deal with any sanity check errors
   //
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) printf("SANITY CHECK FAILED IN OVERLAP TESTS!\n");

      return(2);
   }

#ifdef _DEBUG_
   MPI_Barrier( comm );
   if( irank == 0 ) {
      for(n=0;n<nel1;++n)
         printf("Overlap list (elem %d) post filling for process %d at %p, size %ld\n",
                n, irank, &( lists[n] ), lists[n].size() );
   }
#endif

#ifdef _DEBUG2_
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


   //
   // drop arrays (that now live in expanded form in the list containers)
   // clean-up counters and displacements for re-using
   //
   free( ibuf ); ibuf = NULL;
   free( rbuf ); rbuf = NULL;
   for(n=0;n<nproc;++n) {
      icnt[n] = 0;
      idis[n] = 0;
   }

   //
   // Count per-process data depndencies
   // (The counts correspond to elements of surface2 in each process that
   // overlap with local surface1 elements.)
   //
   for(n=0;n<nel1;++n) {
      std::list< overlap_s > *lp = &( lists[n] );
      std::list< overlap_s > :: iterator il;
      for( il = lp->begin(); il != lp->end(); ++ il ) {
         int ip = (*il).iproc;
         icnt[ip] += 1;
      }
   }
   // build displacement array
   for(n=1;n<nproc;++n) idis[n] = idis[n-1] + icnt[n-1];
   nelem_recv = idis[nproc-1] + icnt[nproc-1];

#ifdef _DEBUG_
   MPI_Barrier( comm );
   if( irank == 0 ) {
      printf("Process %d reporting \n", irank );
      for(n=0;n<nproc;++n)
         printf(" Rank: %d   Displacement: %d   Count: %d \n",
                n, idis[n], icnt[n] );
      printf("Total number of items to receive: %d \n", nelem_recv );
   }
#endif

   //
   // create arrays to use for equivalent sending counts and displacements
   //
   size_t isize = (size_t) nproc;
   isdis = (int *) malloc(isize*sizeof(int));
   iscnt = (int *) malloc(isize*sizeof(int));
   if( isdis == NULL || iscnt == NULL ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) printf("Could not allocate arrays in perform()\n");
      if( isdis != NULL) free( isdis ); isdis = NULL;
      if( iscnt != NULL) free( iscnt ); iscnt = NULL;
      return(-1);
   }

   //
   // negotiate sending counts and form displacement
   //
   ierr = inMPI_DSDE_Global_WindowExchange( nproc, irank, comm, icnt, iscnt );
   if( ierr != 0 ) return(3);

   isdis[0] = 0;
   for(n=0;n<nproc;++n) isdis[n] = isdis[n-1] + iscnt[n-1];
   nelem_send = isdis[nproc-1] + iscnt[nproc-1];

#ifdef _DEBUG_
   MPI_Barrier( comm );
   // in the present test diagnostic irank=1 should have all zeros
   if( irank == 0 ) {
      printf("Process %d reporting \n", irank );
      for(n=0;n<nproc;++n)
         printf(" Rank: %d   Displacement: %d   Count: %d \n",
                n, isdis[n], iscnt[n] );
      printf("Total number of items to send: %d \n", nelem_send );
   }
#endif

   return(0);
}


//
// Setter method to assign the internal function pointer to the external
// overlap-test function.
//
void incg_FaceMatcher::setOverlapFunction( 
                          double (*func)( int n, const double *xyz1,
                                          int m, const double *xyz2 ) )
{
   overlap_func = func;
}


//
// Global pointer variable to an external function to perform the overlap
//
double (*incg_FaceOverlapFunction)( int n, const double *xyz1,
                                    int m, const double *xyz2 ) = NULL;

//
// Getter method to return the receiving count size
//
int incg_FaceMatcher::getSizeRecv( void ) const
{
   return( nelem_recv );
}

//
// Getter method to return the sending count size
//
int incg_FaceMatcher::getSizeSend( void ) const
{
   return( nelem_send );
}


//
// Getter method to return the sending count size
//
int incg_FaceMatcher::formArrays(
                   int *isend_, int *irecv_, double *area_, 
                   int *irdis_, int *ircnt_, int *isdis_, int *iscnt_ )
{
#ifdef _DEBUG_
   if( irank == 0 ) printf("Call to deliver arrays \n");
#endif
   int n,ierr=0;

   if( isdis == NULL || iscnt == NULL ) {
      if( irank == 0 ) {
         printf("Object not in the right state; API calls out of order? \n");
      }
      return(1);
   }

   // allocate an array of requests
   MPI_Request *ireq = (MPI_Request *) malloc( nproc*2*sizeof(MPI_Request) );
   if( ireq == NULL ) ierr = 1;
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) {
         printf("Could not allocate memory for MPI requests array (strange)\n");
      }
      if( ireq != NULL ) free( ireq );
      return(-1);
   }

   // copy receive index array and area/weights array
   for(n=0;n<nproc;++n) icnt[n] = 0;   // use as counters
   for(n=0;n<nel1;++n) {
      std::list< overlap_s > *lp = &( lists[n] );
      std::list< overlap_s > :: iterator il;
      for( il = lp->begin(); il != lp->end(); ++il ) {
         int ip = (*il).iproc;

         irecv_[ idis[ip] + icnt[ip] ] = (*il).ielem;
         area_[ idis[ip] + icnt[ip] ] = (*il).a;

         icnt[ip] += 1;
      }
   }
   // sanity check
   for(n=1;n<nproc;++n) {
      if( idis[n-1] + icnt[n-1] != idis[n] ) ierr = 1;
   }
   MPI_Allreduce( MPI_IN_PLACE, &ierr, 1, MPI_INT, MPI_SUM, comm );
   if( ierr != 0 ) {
      if( irank == 0 ) {
         printf("A sanity check did not pass; contact IN immediately!\n");
      }
      return(2);
   }

/*
   // negotiate sending arrays
   for(n=0;n<nproc;++n) {
      MPI_Irecv( &( isend_[ isdis[n] ] ), iscnt[n],
                 MPI_INT, n, 1000+n, comm, &( ireq[nproc+n] ) );
   }
   for(n=0;n<nproc;++n) {
      MPI_Issend( &( irecv_[ idis[n] ] ), icnt[n],
                 MPI_INT, n, 1000+irank, comm, &( ireq[n] ) );
   }
   MPI_Waitall( 2*nproc, ireq, MPI_STATUS_IGNORE );
*/

   // drop the requests array
   free( ireq );

   // copy displacement and count arrays
   for(n=0;n<nproc;++n) {
      irdis_[n] = idis[n];
      ircnt_[n] = icnt[n];
      isdis_[n] = isdis[n];
      iscnt_[n] = iscnt[n];
   }

#ifdef _DEBUG_
for(n=0;n<nproc;++n) { if( irank == n ) {
char filename[20];
sprintf( filename, "TEST_%.5d.dat", n );
FILE *fp = fopen( filename, "w" );
fprintf(fp,"### Sending groups ### \n");
for(int m=0;m<nproc;++m) {
fprintf(fp,"Receiving procees: %d \n",m);
for(int i=isdis[m];i<isdis[m]+iscnt[m];++i) fprintf(fp," %d: %d \n",m,i);
}
fclose(fp);
} MPI_Barrier( comm ); }
#endif

   return(0);
}


//
// API function to invoke the functionality from the C side
// (This function is meant to test the first stage of the object's operation
// and is not going to be used when this functionality is made into a library.
// It is supposed to perform the basic matching operations and then clean up.)
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
   incg_FaceMatcher *fm = new incg_FaceMatcher( comm );

   // provide variables to the object
   fm->setData( icnt1, ilist, x1,  icnt2, jlist, x2 );

   // call this if we are using some pre-search-based acceleration technique
   (void) fm->setAccel( NULL );
   // (this is a test)
int ijunk[] = {
 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
   (void) fm->setAccel( ijunk );

   // method to create internal structures
   ierr = fm->prepare();
   if( ierr != 0 ) {
      delete fm;
      return(1);
   }

   // set pointer to external function to perform face-matching
   // (See notes on how to set externally the variable that is the argument!)
   fm->setOverlapFunction( incg_FaceOverlapFunction );

   // method to perform face-matching
   ierr = fm->perform();
   if( ierr != 0 ) {
      delete fm;
      return(2);
   }

   // drop object
   delete fm;

   return(0);
}


//
// Global vector of face-matcher objects
//
std::vector< incg_FaceMatcher * > global_fm_objects;


//
// API function prototype to instantiate an object and return a handle
//
int incg_Facematch_Init( int *handle, MPI_Comm *comm,
                                           int icnt1, int *ilist, double *x1,
                                           int icnt2, int *jlist, double *x2,
                                           int *iacc )
{
   int nproc,irank;
   MPI_Comm_size( *comm, &nproc );
   MPI_Comm_rank( *comm, &irank );
   printf("The function was invoked by MPI process %d \n", irank);
   int ierr=0;


   // constructing the object
   incg_FaceMatcher *fm = new incg_FaceMatcher( comm );

   // provide variables to the object
   fm->setData( icnt1, ilist, x1,  icnt2, jlist, x2 );

   // call this if we are using some pre-search-based acceleration technique
   (void) fm->setAccel( iacc );
   // (this is a test)
int ijunk[] = {
 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
   (void) fm->setAccel( ijunk );

   // method to create internal structures
   ierr = fm->prepare();
   if( ierr != 0 ) {
      delete fm;
      return(1);
   }

   // set pointer to external function to perform face-matching
   // (See notes on how to set externally the variable that is the argument!)
   fm->setOverlapFunction( incg_FaceOverlapFunction );

   // method to perform face-matching
   ierr = fm->perform();
   if( ierr != 0 ) {
      delete fm;
      return(2);
   }


   //
   // place the object in the vector
   //
   int iplace=-1,k=0;
   while( k < (int) (global_fm_objects.size()) ) {
      if( global_fm_objects[k] == NULL ) iplace = k;
      ++k;
   }
   if( iplace == -1 ) {
      global_fm_objects.push_back( fm );
      *handle = (int) (global_fm_objects.size()) - 1;
   } else {
      global_fm_objects[ iplace ] = fm;
      *handle = iplace;
   }

   return(0);
}


//
// API function to deconstruct a handle's object
//
int incg_Facematch_Term( int *handle )
{
   if( *handle >= (int) (global_fm_objects.size()) ) return(1);

   if( global_fm_objects[ *handle ] == NULL ) return(2);

   incg_FaceMatcher *fm = global_fm_objects[ *handle ];
   delete fm;

   global_fm_objects[ *handle ] = NULL;

   return(0);
}

//
// API function to retreive sizes from the handle's object
//
int incg_Facematch_GetSizes( int *handle, int *num_recv, int *num_send )
{
   if( *handle >= (int) (global_fm_objects.size()) ) return(1);

   if( global_fm_objects[ *handle ] == NULL ) return(2);

   incg_FaceMatcher *fm = global_fm_objects[ *handle ];

   *num_recv = fm->getSizeRecv();
   *num_send = fm->getSizeSend();

   return(0);
}

//
// API function prototype to retreive array data from the handle's object
//
int incg_Facematch_FillArrays( int *handle,
                               int *isend, int *irecv, double *recv_area,
                               int *irdis, int *ircnt, int *isdis, int *iscnt )
{
   if( *handle >= (int) (global_fm_objects.size()) ) return(1);

   if( global_fm_objects[ *handle ] == NULL ) return(2);

   incg_FaceMatcher *fm = global_fm_objects[ *handle ];


   // get the array data
   fm->formArrays( isend, irecv, recv_area, irdis, ircnt, isdis, iscnt );


   return(0);
}





}  // extern C
#endif

