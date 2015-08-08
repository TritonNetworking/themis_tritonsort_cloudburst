#include<stdint.h>
#include<string.h>
#include "AlignInfo.h"

AlignInfo::AlignInfo() {
  setVals(0, 0, NULL, NULL, 0);
}

//  constructor
AlignInfo::AlignInfo(
  int32_t len, int32_t k, int32_t* pdist, int32_t* pweightMatrix,
  int32_t dlen) {
  setVals(len, k, pdist, pweightMatrix, dlen);
}

//  setvals
void AlignInfo::setVals(
  int32_t len, int32_t alignDiff, int32_t* pdist, int32_t* pweightMatrix,
  int32_t dlen) {
  alignlen = len;
  differences = alignDiff;
  dist = pdist;
  weightMatrix = pweightMatrix;
  distlen = dlen;
}

// isBazeaYatesSeed 
// Since an alignment may be recompute k+1 times for each of the k+1 seeds,
// see if the current alignment is the leftmost alignment by checking for
// differences in the proceeding chunks of the query

bool AlignInfo::isBazeaYatesSeed(int32_t qlen, int32_t kmerlen) {
  int32_t numBuckets = qlen / kmerlen;
  int32_t lastbucket = -1;
  int32_t distdelta = 0;
  int32_t pos = 0;
  for (int32_t i = 0; i < distlen; i++) {
    pos += dist[i] + distdelta;
    distdelta = 0;
    if (weightMatrix[i] == 2) {
      // end of string
      continue;
    } else if (weightMatrix[i] == -1) {
      // gap character occurs between pos and pos+1
      if (pos % kmerlen == 0) {
        // occurs right between buckets, skip
        continue;
      }
    }
    int32_t bucket = pos / kmerlen;
    if (bucket - lastbucket > 1) {
      return false;
    }
    lastbucket = bucket;
  }
  return (lastbucket == numBuckets-1);
}

