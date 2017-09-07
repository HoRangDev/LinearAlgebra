#include "GaussianElimination.h"

int main( )
{
   LA::Matrix<3, 3> mat{ };
   mat[ 0 ] = { 2, -1, 1 }; 
   mat[ 1 ] = { 1, 1, 1 };
   mat[ 2 ] = { 1, 1, 2 };

   LA::Vector<3> vec{ 3,6,9 };

   auto res = LA::ForwardElimination( mat, vec );
   res = res;

   auto sol = LA::GaussianElimination( mat, vec );
   sol = sol;

   return 0;
}