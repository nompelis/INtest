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


#ifdef _DEBUG_
   if( irank == 0 ) printf("incg_FaceMatcher object instantiated\n");
#endif
}

//
// Deconstructor of face-matcher object
//
incg_FaceMatcher::~incg_FaceMatcher( )
{


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
// API function to invoke the functionality from the C side
//
int incg_PerformFacematch( MPI_Comm *comm, int icnt1, int *ilist, double *x1,
                                           int icnt2, int *jlist, double *x2 )
{
   int nproc,irank;
   MPI_Comm_size( *comm, &nproc );
   MPI_Comm_rank( *comm, &irank );
   printf("The function was invoked by MPI process %d \n", irank);


   // constructing the object
   incg_FaceMatcher fm( comm );

   // provide variables to the object
   fm.setData( icnt1, ilist, x1,  icnt2, jlist, x2 );




   return(0);
}


}
#endif

