#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>



int main( int argc, char *argv[] )
{
   int irequired,iprovided;
   int nproc,irank;
   MPI_Comm comm;
   MPI_Info mpi_info;
   char hname[1000];
   int ilen;


   irequired = 0;
   comm = MPI_COMM_WORLD;
   MPI_Init_thread( &argc, &argv, irequired, &iprovided );
   MPI_Comm_size( comm, &nproc );
   MPI_Comm_rank( comm, &irank );
   MPI_Info_create( &mpi_info );
   MPI_Get_processor_name( hname, &ilen );


   if( irank == 0 ) printf("Test running on %d proccesses on \"%s\" \n", nproc,hname);



   MPI_Finalize();

   return(0);
}

