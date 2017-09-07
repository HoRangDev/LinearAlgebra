#pragma once

#include <vector>
#include <array>
#include <cmath>

namespace LA
{
   template <size_t M, size_t N>
   using Matrix = std::array<std::array<double, N>, M>;

   template <size_t N>
   using Vector = std::array<double, N>;

   template < size_t M, size_t N >
   Matrix<M, N + 1>&& MakeAugumentedMatrix( const Matrix<M, N>& mat, const Vector<M>& constant )
   {
      Matrix<M, N + 1> augumentedMat{ };
      for ( size_t row = 0; row < M; ++row )
      {
         // Given matrix to augumented matrix.
         for ( size_t column = 0; column < N; ++column )
         {
            augumentedMat[ row ][ column ] = mat[ row ][ column ];
         }

         // Constant vector to augumented matrix.
         augumentedMat[ row ][ N ] = constant[ row ];
      }

      return std::move( augumentedMat );
   }
}