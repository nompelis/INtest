
#ifndef _INCG_FACEMATCH_H_
#define _INCG_FACEMATCH_H_

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#ifdef __cplusplus
#include <vector>
#include <list>
#endif

#ifdef __cplusplus
extern "C" {
#endif

//
// API function prototype visible to the C side to invoke the functionality
//
int incg_PerformFacematch( MPI_Comm *comm, int icnt1, int *ilist, double *x1,
                                           int icnt2, int *jlist, double *x2 );

//
// API function prototype to instantiate an object and return a handle
//
int incg_Facematch_Init( int *handle, MPI_Comm *comm,
                                           int icnt1, int *ilist, double *x1,
                                           int icnt2, int *jlist, double *x2,
                                           int *iacc );

//
// API function prototype to deconstruct a handle's object
//
int incg_Facematch_Term( int *handle );

//
// API function prototype to retreive sizes from the object
//
int incg_Facematch_GetSizes( int *handle, int *num_recv, int *num_send );

//
// API function prototype to retreive array data from the object
//
int incg_Facematch_FillArrays( int *handle,
                               int *isend, int *irecv, double *recv_area,
                               int *irdis, int *ircnt, int *isdis, int *iscnt );

#ifdef __cplusplus
}
#endif

//
// a structure to hold overlap data
//
struct overlap_s {
   int iproc;
   int ielem;
   double a;
};


#ifdef __cplusplus
extern "C" {


//
// The object that performs face-matching operations and creates appropriate
// arrays
//
class incg_FaceMatcher {
 public:
   incg_FaceMatcher( MPI_Comm *comp );
   ~incg_FaceMatcher();

   void setData( int icnt_, int *icon_, double *xi_,
                 int jcnt_, int *jcon_, double *xj_ );
   void setAccel( int *ibs_ );

   int prepare( void );
   int perform( void );
   void setOverlapFunction( double (*func)( int n, const double *xyz1,
                                            int m, const double *xyz2 ) );

   int getSizeRecv( void ) const;
   int getSizeSend( void ) const;

   int formArrays( int *isend_, int *irecv_, double *area_,
                   int *irdis_, int *ircnt_, int *isdis_, int *iscnt_ );

 protected:

 private:
   // communication variables and handles
   MPI_Comm comm;
   int irank;
   int nproc;

   // incoming data and pointers
   int nel1,nel2;
   int *icon,*jcon;
   double *xi,*xj;

   // internal arrays, etc
   int *idst,*jdst;
   int *iacc,*idis,*icnt;
   int nbuf_size, *ibuf;
   double *rbuf;
   int *idata;
   double *rdata;
   std::vector< std::list< overlap_s > > lists;
   double (*overlap_func)( int n, const double *xyz1,
                           int m, const double *xyz2 );

   int nelem_recv, nelem_send;
   int *isdis,*iscnt;
};



}
#endif

#endif

