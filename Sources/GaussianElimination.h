#include "Common.h"

namespace LA
{
   template <size_t M = 3, size_t N = 3>
   Matrix<M, N + 1>&& ForwardElimination( const Matrix<M, N>& mat, const Vector<M>& constant )
   {
      auto augumentedMat = MakeAugumentedMatrix<M, N>( mat, constant );
      for ( size_t unknownIdx = 0; unknownIdx < ( N - 1 ); ++unknownIdx )
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

   template <size_t M = 3, size_t N = 3>
   Vector<N>&& GaussianElimination( const Matrix<M, N>& mat, const Vector<M>& constant )
   {
      // mat and constant no more using.
      auto resultMat = ForwardElimination( std::move( mat ), std::move( constant ) );
      auto solution = Vector<N>{ };

      for ( long long unknownIdx = N-1; unknownIdx >= 0; --unknownIdx )
      {
         for ( long long col = N-1; col >= unknownIdx; --col )
         {
            if ( col == unknownIdx )
            {
               solution[ unknownIdx ] 
                  = resultMat[ unknownIdx ][ N ] / resultMat[ unknownIdx ][ unknownIdx ];
            }
            else
            {
               resultMat[ unknownIdx ][ N ] -= resultMat[ unknownIdx ][ col ] * solution[ col ];
            }
         }
      }

      return std::move( solution );
   }
};