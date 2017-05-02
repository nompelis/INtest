#include <stdio.h>
#include <stdlib.h>

// This file is going to be used by the library API for creating fortran
// bindings. It is needed in order to introduce the overlap function to the
// API form the Fortran side. In particular, it is used to set the global
// pointer variable to point to the overlap function. To avoid conflicts with
// the functionality on the pure C side of the API, perhaps the user overlap
// function can be #include-ed in this file and on the C side file. A permanent
// solution would be to include the overlap function in this library.

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
// Function to be called from fortran in order to set the overlap function
//
void incg_facematch_setoverlapfunction_()
{
   incg_FaceOverlapFunction = &my_overlap_function;
}

