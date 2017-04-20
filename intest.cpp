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
#ifdef _DEBUG_
   printf("incg_FaceMatcher object instantiated\n");
#endif

   MPI_Comm_dup( *comp, &comm );
   MPI_Comm_size( comm, &nproc );
   MPI_Comm_rank( comm, &irank );


}

//
// Deconstructor of face-matcher object
//
incg_FaceMatcher::~incg_FaceMatcher( )
{


   MPI_Comm_free( &comm );
}



//
// API function to invoke the functionality from the C side
//
int incg_PerformFacematch( MPI_Comm *comm, int icnt1, double *x1,
                                           int icnt2, double *x2 )
{
   int nproc,irank;
   MPI_Comm_size( *comm, &nproc );
   MPI_Comm_rank( *comm, &irank );
   printf("The function was invoked by MPI process %d \n", irank);


   incg_FaceMatcher fm( comm );





   return(0);
}


}
#endif

