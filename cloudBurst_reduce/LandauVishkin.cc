#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <cstdlib>
#include <assert.h>

#include "LandauVishkin.h"
#include "core/MemoryUtils.h"

//  For a number of differences e and diagonal d
//  L[d][e] denote largest row i such that D[i][l] =e  and D[i][l]
//  on diagonal D
//  for more details refer Landau vishkin 1986,88 paper

void LandauVishkin::configure(int32_t _maxAlignDiff) {
  // initialize 2d matrix
  maxAlignDiff = _maxAlignDiff;
  matrix  = new (themis::memcheck) int32_t* [maxAlignDiff*2+1];
  diagMatrix  = new (themis::memcheck) int32_t* [maxAlignDiff*2+1];
  for (int32_t i = 0 ; i < 2*maxAlignDiff+1; i++) {
    matrix[i] = new (themis::memcheck) int32_t[maxAlignDiff+1];
    diagMatrix[i] = new (themis::memcheck) int32_t[maxAlignDiff+1];
  }
  dist = new (themis::memcheck) int32_t[maxAlignDiff+1];
  what = new (themis::memcheck) int32_t[maxAlignDiff+1];
  noAlignment.setVals(0, 0, NULL, NULL, 0);
  badAlignment.setVals(-1, -1, NULL, NULL, 0);
  goodAlignment.setVals(0, 0, NULL, NULL, 0);
}

LandauVishkin::~LandauVishkin() {
  for (int32_t i = 0 ; i < 2*maxAlignDiff+1; i++) {
    delete[] matrix[i];
    delete[] diagMatrix[i];
  }
  delete[] matrix;
  delete[] diagMatrix;
  delete[] dist;
  delete[] what;
}

//  kdifference
//  Landau-Vishkin k-difference algorithm to align strings
//  dynamic programming with matrix
AlignInfo& LandauVishkin::kdifference(
  const byte* text, int32_t textLength, const byte* pattern,
  int32_t patternLength, int k) {
  if (patternLength == 0 || textLength == 0) {
    return noAlignment;
  }
  // Compute the dynamic programming to see how the strings align
  for (int32_t numDiffE = 0; numDiffE <= k; numDiffE++) {
    for (int32_t diagonalD = -numDiffE; diagonalD <= numDiffE; diagonalD++) {
      int32_t row = -1;

      if (numDiffE > 0) {
        if (abs(diagonalD) < numDiffE) {
          int32_t up = matrix[k+diagonalD][numDiffE-1] + 1;
          if (up > row) {
            row = up;
            diagMatrix[k+diagonalD][numDiffE] = 0;
          }
        }

        if (diagonalD > -(numDiffE-1)) {
          int32_t left = matrix[k+diagonalD-1][numDiffE-1];
          if (left > row) {
            row = left;
            diagMatrix[k+diagonalD][numDiffE] = -1;
          }
        }

        if (diagonalD < numDiffE-1) {
          int32_t right = matrix[k+diagonalD+1][numDiffE-1]+1;
          if (right > row) {
            row = right;
            diagMatrix[k+diagonalD][numDiffE] = +1;
          }
        }
      } else {
        row = 0;
      }

      while ((row < patternLength) && (row+diagonalD < textLength)
        && (pattern[row] == text[row+diagonalD])) {
        row++;
      }

      matrix[k+diagonalD][numDiffE] = row;
      if ((row+diagonalD == textLength) || (row == patternLength)) {
        // reached the end of the pattern or text
        int32_t distlen = numDiffE+1;

        int32_t E = numDiffE;
        int32_t D = diagonalD;

        what[numDiffE] = 2;  // always end at end-of-string

        while (numDiffE >= 0) {
          int32_t b = diagMatrix[k+diagonalD][numDiffE];
          if (numDiffE > 0) {
            what[numDiffE-1] = b;
          }

          dist[numDiffE] = matrix[k+diagonalD][numDiffE];
          if (numDiffE < E) {
            dist[numDiffE+1] -= dist[numDiffE];
          }

          diagonalD += b;
          numDiffE--;
        }
        goodAlignment.setVals(row+D, E, dist, what, distlen);
        // say how far we reached in the text (reference)
        return goodAlignment;
      }
    }
  }
  return badAlignment;
}
//  extend
//  align the strings either for either k-mismatch or k-difference
AlignInfo& LandauVishkin::extend(
  const byte* refbin, int32_t refLength, const byte* qrybin,
  int32_t qryLength, int32_t k, bool allowDiff) {
  if (allowDiff) {
    if(refLength < 1 || qryLength < 1) {
      return badAlignment;
    }
    byte* ref = dnaStringObj.dnaToArr(refbin, refLength);
    byte* qry = dnaStringObj.dnaToArr(refbin, qryLength);
    AlignInfo& a = kdifference(ref, refLength, qry, qryLength, k);
    delete[] ref;
    delete[] qry;
    return a;
  } else {
    return kmismatch_bin(refbin, refLength, qrybin, qryLength, k);
  }
}
