#ifndef CLOUD_BURST_REDUCE_FUNCTION_H
#define CLOUD_BURST_REDUCE_FUNCTION_H

#include <vector>

#include "mapreduce/functions/reduce/ReduceFunction.h"
#include "mapreduce/functions/map/cloudBurst/DNAString.h"
#include "mapreduce/functions/map/cloudBurst/MerRecord.h"
#include "mapreduce/functions/reduce/cloudBurst/AlignmentRecord.h"
#include "mapreduce/functions/reduce/cloudBurst/LandauVishkin.h"

/**
   CloudBurst reduce function and associated helper classes based on the
   open-source Java implementation of CloudBurst by Michael Schatz.
   See http://sourceforge.net/p/cloudburst-bio/wiki/CloudBurst/
   ----

   The input to the CloudBurstReduceFunction is two sets of tuples: a reference
   set and a query set. The last byte of each key is used to distinguish between
   reference and query, so that the input to a single invocation of reduce() is
   either a set of reference records or a set of query records.

   The goal of this function is to output alignment information if a query
   record has a small enough number of differences from a reference record.
   Query records are compared to reference records with the same seed, which is
   a subsequence that is guaranteed to be the same given the constraints on
   number of differences. Seeds are stored in the first |key| - 2 bytes of each
   key.

   A tuple's value contains information about the sequences on either side of
   the seed, which are referred to as flanks. By definition, flanks are the
   parts of the sequences that are allowed to differ. Seeds are extended with
   flanks in the extend() function and aligned in the alignBatch() function.
 */
class CloudBurstReduceFunction : public ReduceFunction {
public:
  /// Constructor
  /**
     \param maxAlignDiff the maximum number of mismatches to allow

     \param seedLength the length of the seed, which is the portion of the DNA
     sequence that must match exactly

     \param allowDifferences whether or not indels (insertions/deletions) should
     be allowed, or simply mismatches

     \param blockSize the number of query tuples to batch up before aligning

     \param redundancy the number of copies of low complexity seeds to use
   */
  CloudBurstReduceFunction(
    uint32_t maxAlignDiff, uint32_t seedLength, uint32_t allowDifferences,
    uint32_t blockSize, uint32_t redundancy);

  /// \sa ReduceFunction::reduce
  void reduce(
    const uint8_t* key, uint64_t keyLength,
    KVPairIterator& iterator, KVPairWriterInterface& writer);

  /**
     Clear state when reducing a new buffer since old pointers may be invalid.
   */
  void configure();

private:
  /**
     Clear state so a new seed can be aligned.
   */
  void clearState();

  /**
     Align stored query tuples to stored reference tuples with the same seed,
     and write out any matches that are within the maximum number of
     differences.

     \param writer the KVPairWriterInterface associated with the Reducer
   */
  void alignBatch(KVPairWriterInterface& writer);

  /**
     Extend seeds using the flank information in the tuple value, and compare
     the extended query and reference records.

     \param queryTuple the MerRecord representation of a query tuple to compare

     \param referenceTuple the MerRecord representation of a reference tuple to
     compare

     \return an AlignmentRecord that has the differences field set to -1 if the
     records do NOT align
   */
  AlignmentRecord* extend(
    const MerRecord& queryTuple, const MerRecord& referenceTuple);

  AlignmentRecord noalignment;
  AlignmentRecord fullalignment;
  LandauVishkin  landauVishkinObj;
  std::vector<MerRecord> referenceTuples;
  std::vector<MerRecord> queryTuples;
  DNAString   dnaStringObj;
  uint32_t maxAlignDiff;
  uint32_t seedLength;
  uint32_t blockSize;
  uint32_t redundancy;
  bool allowDifferences;
  bool filterAlignments;
  bool readingReferenceTuples;
  const uint8_t* referenceKey;
  uint32_t referenceKeyLength;
};

#endif // CLOUD_BURST_REDUCE_FUNCTION_H
