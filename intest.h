
#ifndef _INTEST_H_
#define _INTEST_H_

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


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


#ifdef __cplusplus
extern "C" {

class incg_FaceMatcher {
 public:
   incg_FaceMatcher( MPI_Comm *comp );
   ~incg_FaceMatcher();

   void setData( int icnt_, int *icon_, double *xi_,
                 int jcnt_, int *jcon_, double *xj_ );

 protected:

 private:
   MPI_Comm comm;
   int irank;
   int nproc;

   int nel1,nel2;
   int *icon,*jcon;
   double *xi,*xj;

};



}
#endif

#endif

