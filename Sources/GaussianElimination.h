#include "Common.h"

namespace LA
{
   template <size_t M = 3, size_t N = 3>
   Matrix<M, N + 1>&& ForwardElimination( const Matrix<M, N>& mat, const Vector<M>& constant )
   {
      auto augumentedMat = MakeAugumentedMatrix<M, N>( mat, constant );
      for ( size_t unknownIdx = 0; unknownIdx < ( M - 1 ); ++unknownIdx )
      {
         // Compute new rows for forward Elimination
         for ( size_t row = unknownIdx + 1; row < M; ++row )
         {
            double proportion 
               = ( augumentedMat[ unknownIdx ][ unknownIdx ] / augumentedMat[ row ][ unknownIdx ]);
            for ( size_t col = 0; col <= N; ++col )
            {
               augumentedMat[ row ][ col ] *= proportion;
            }
         }

         // Compute Forward Elimination
         for ( size_t row = unknownIdx + 1; row < M; ++row )
         {
            for ( size_t col = 0; col <= N; ++col )
            {
               augumentedMat[ row ][ col ] -= augumentedMat[ unknownIdx ][ col ];
            }
         }
      }

      return std::move( augumentedMat );
   }
};