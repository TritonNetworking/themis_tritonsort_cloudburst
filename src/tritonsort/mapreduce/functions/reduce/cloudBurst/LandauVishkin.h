#ifndef _LANDAU_VISHKIN_H
#define _LANDAU_VISHKIN_H

#include "AlignInfo.h"
#include "AlignmentRecord.h"
#include "mapreduce/functions/map/cloudBurst/DNAString.h"

class LandauVishkin {
public:
  AlignInfo noAlignment;
  AlignInfo badAlignment;
  AlignInfo goodAlignment;
  DNAString dnaStringObj;
  int32_t** matrix;  // 2d array [2K+1][K+1]size
  int32_t** diagMatrix;  // 2d array [2K+1][K+1]size
  int32_t* dist;  // 2d array [2K+1][K+1]size
  int32_t* what;  // 2d array [2K+1][K+1]size
  int32_t maxAlignDiff;
  // initialize runtime buffers
  virtual ~LandauVishkin();
  void configure(int32_t k);

  AlignInfo& kdifference(
    const byte* text, int32_t textLength, const byte* pattern,
    int32_t patternLength, int32_t k);

  AlignInfo& extend(
    const byte* refbin, int32_t refLength, const byte* qrybin,
    int32_t qryLength, int32_t K, bool ALLOW_DIFFERENCES);

  // Ever so slightly optimized version of kmismatch_bin, inlined for
  // performance
  // text corresponds to reference and pattern corresponds to query
  inline AlignInfo& kmismatch_bin(
    const byte* text, int32_t textLength, const byte* pattern,
    int32_t patternLength, int32_t k) {
    if (patternLength == 0) {
      // 0-length query patterns should be considered as if they did not align.
      return noAlignment;
    } else if (textLength < patternLength) {
      // The reference (text) must be long enough to fully contain the query.
      return badAlignment;
    }

    int32_t alignmentLength = std::min<int32_t>(textLength, patternLength);

    int32_t distanceFromLastDifference = 0;
    int32_t differences = 0;

    const byte* lastTextByte = text + alignmentLength - 1;
    // Convert alignment length to base pairs.
    alignmentLength *= 2;

    // Check the last byte separately since it might contain an odd number of
    // base pairs.
    while (text < lastTextByte) {
      // Each byte encodes 2 sets of base pairs, so we need to compare 4 bit
      // nibbles. The most efficient way to do this is to XOR the two bytes and
      // check if a left or right shift is 0.
      uint8_t byteXOR = *text ^ *pattern;

      // Check high bits
      if ((byteXOR >> 4) != 0) {
        if (differences >= k) {
          return badAlignment;
        }
        dist[differences] = distanceFromLastDifference;
        distanceFromLastDifference = 1;
        ++differences;
      } else {
        ++distanceFromLastDifference;
      }

      // Check low bits
      if ((byteXOR & 0x0F) != 0) {
        if (differences >= k) {
          return badAlignment;
        }
        dist[differences] = distanceFromLastDifference;
        distanceFromLastDifference = 1;
        ++differences;
      } else {
        ++distanceFromLastDifference;
      }

      // Advance to the next byte.
      ++text;
      ++pattern;
    }

    // Check the last byte, which might have a space in the second nibble.
    uint8_t byteXOR = *text ^ *pattern;
    if ((*text & 0x0F) == dnaStringObj.space ||
        ((*pattern & 0x0F) == dnaStringObj.space)) {
      // One of the strings has a space in the second nibble, so only compare the
      // first set of base pairs.

      // Check high bits
      if ((byteXOR >> 4) != 0) {
        if (differences >= k) {
          return badAlignment;
        }
        dist[differences] = distanceFromLastDifference;
        distanceFromLastDifference = 1;
        ++differences;
      } else {
        ++distanceFromLastDifference;
      }

      // Since the last base pair set is ignored, reduce the alignment by 1.
      --alignmentLength;
    } else {
      // Check high bits
      if ((byteXOR >> 4) != 0) {
        if (differences >= k) {
          return badAlignment;
        }
        dist[differences] = distanceFromLastDifference;
        distanceFromLastDifference = 1;
        ++differences;
      } else {
        ++distanceFromLastDifference;
      }

      // Check low bits
      if ((byteXOR & 0x0F) != 0) {
        if (differences >= k) {
          return badAlignment;
        }
        dist[differences] = distanceFromLastDifference;
        distanceFromLastDifference = 1;
        ++differences;
      } else {
        distanceFromLastDifference++;
      }
    }

    // Reaching this point means there were at most k mismatches.
    memset(what, 0, differences * sizeof(int32_t));

    // Record the number of matches until the end of the alignment.
    dist[differences] = distanceFromLastDifference;
    what[differences] = 2;

    goodAlignment.setVals(
      alignmentLength, differences, dist, what, differences + 1);
    return goodAlignment;
  }
};
#endif
