
#ifndef _INTEST_H_
#define _INTEST_H_

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

#ifdef __cplusplus
}
#endif

//
// a structure to holde overlap data
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
};



}
#endif

#endif

