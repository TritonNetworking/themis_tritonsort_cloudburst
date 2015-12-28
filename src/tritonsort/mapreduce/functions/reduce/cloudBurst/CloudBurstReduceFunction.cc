#include "CloudBurstReduceFunction.h"
#include "core/MemoryUtils.h"
#include "mapreduce/common/KeyValuePair.h"
#include "mapreduce/functions/reduce/cloudBurst/AlignInfo.h"

CloudBurstReduceFunction::CloudBurstReduceFunction(
  uint32_t _maxAlignDiff, uint32_t _seedLength, uint32_t _allowDifferences,
  uint32_t _blockSize, uint32_t _redundancy)
  : noalignment(-1, -1, -1, -1, true),
    maxAlignDiff(_maxAlignDiff),
    seedLength(_seedLength),
    blockSize(_blockSize),
    redundancy(_redundancy),
    allowDifferences(_allowDifferences),
    readingReferenceTuples(false),
    referenceKey(NULL),
    referenceKeyLength(0) {
  landauVishkinObj.configure(maxAlignDiff);
}

void CloudBurstReduceFunction::reduce(
  const uint8_t* key, uint64_t keyLength,
  KVPairIterator& iterator, KVPairWriterInterface& writer) {

  uint64_t tuplesRead = 0;
  MerRecord merIn;
  // Reduce::
  // Pair up Reference and Query tuples which has same key (-minus last 2 byte)
  // Last byte differentiates whether tuple is refTuple or queryTuple
  // Group these tuples together and send to alignBatch function
  // which extends these strings with max K differences

  // Reference mers are first, save them away

  uint8_t lastKeyByte = key[keyLength - 1];
  if (lastKeyByte == 0) {
    // Reference tuples should have a 0 in the last byte of the key.
    readingReferenceTuples = true;
    // Since we expect reference tuples before query tuples, we know this is a
    // new seed, so we can clear out any old records.
    referenceTuples.clear();
    ASSERT(queryTuples.empty(),
           "Got a new reference seed but the set of query tuples is non-empty. "
           "Query tuples should be written and cleared at the end of reduce()");
    // Set the reference key so we can check if future query records have the
    // same seed.
    referenceKey = key;
    referenceKeyLength = keyLength;
  } else if (lastKeyByte == 1) {
    // Query tuples should have a 1 in the last byte of the key.
    readingReferenceTuples = false;
    if (referenceTuples.empty()) {
      // There are no reference tuples for this seed, so just return.
      return;
    } else {
      // Verify that the query seed matches the previous reference seed, which
      // is the key minus the last 2 bytes.
      ASSERT(referenceKey != NULL,
             "Reference key should not be NULL if reference records exist.");
      if (keyLength != referenceKeyLength ||
          memcmp(key, referenceKey, keyLength - 2) != 0) {
        // These query records don't correspond to the previous reference
        // records, so clear the set of reference records and return.
        referenceTuples.clear();
        return;
      }
    }
  } else {
    ABORT("Last byte of the key should be 0 (reference) or 1 (query). Got %u",
          lastKeyByte);
  }

  KeyValuePair tuple;

  while (iterator.next(tuple)) {
    merIn.fromBytes(tuple.getValue(), tuple.getValueLength());

    if (merIn.isReference) {
      ABORT_IF(!readingReferenceTuples,
               "Got a reference tuple (tuple %llu) but expected only query "
               "tuples", tuplesRead);

      // Store reference tuples in a vector because we don't expect to see any
      // query tuples until the next invocation of reduce()
      referenceTuples.push_back(merIn);
    } else {
      ABORT_IF(readingReferenceTuples,
               "Got a query tuple (tuple %llu) but expected only reference "
               "tuples", tuplesRead);

      // Store query tuples until we have a batch of size 'blockSize'
      queryTuples.push_back(merIn);
      if (queryTuples.size() == blockSize) {
        // Perform DNA alignment and write out a batch of records.
        alignBatch(writer);
        queryTuples.clear();
      }
    }

    ++tuplesRead;
  }

  // Write out any left over records.
  if (!queryTuples.empty()) {
    alignBatch(writer);
  }

  // If we just finished aligning query tuples to reference tuples, clear out
  // everything since we are expecting a new set of reference tuples.
  if (!readingReferenceTuples) {
    clearState();
  }
}

void CloudBurstReduceFunction::configure() {
  // Clear state between buffers.
  clearState();
}

void CloudBurstReduceFunction::clearState() {
  referenceTuples.clear();
  queryTuples.clear();
  referenceKey = NULL;
  referenceKeyLength = 0;
}

void CloudBurstReduceFunction::alignBatch(KVPairWriterInterface& writer) {
  int32_t numRefTuples = referenceTuples.size();
  int32_t numQueryTuples = queryTuples.size();

  // join together the query-ref shared mers
  if ((numRefTuples != 0) && (numQueryTuples != 0)) {
    // Align reads to the references in blocks of blockSize x BLOCK_SIZE
    // to improve cache locality
    // define a qry block between [queryTuplesIndex, lastQueryTupleIndex)
    for (int32_t queryTuplesIndex = 0; queryTuplesIndex < numQueryTuples;
      queryTuplesIndex += blockSize) {
      int32_t lastQueryTupleIndex = queryTuplesIndex + blockSize;
      if (lastQueryTupleIndex > numQueryTuples) {
        lastQueryTupleIndex = numQueryTuples;
      }
      // define a ref block between [startRefTupleIndex, lastRefTupleIndex)
      for (int32_t startRefTupleIndex = 0; startRefTupleIndex < numRefTuples;
        startRefTupleIndex += blockSize) {
        int32_t lastRefTupleIndex = startRefTupleIndex + blockSize;
        if (lastRefTupleIndex > numRefTuples) {
          lastRefTupleIndex = numRefTuples;
        }
        // for each element in [queryTuplesIndex, lastQueryTupleIndex)
        for (int32_t curq = queryTuplesIndex; curq < lastQueryTupleIndex;
          curq++) {
          MerRecord& queryRecord = queryTuples[curq];
          int32_t queryID = queryRecord.id;
          // for each element in [startRefTupleIndex, lastRefTupleIndex)
          for (int32_t curr = startRefTupleIndex; curr < lastRefTupleIndex;
            curr++) {
            AlignmentRecord* rec =
              extend(queryRecord, referenceTuples[curr]);
            if (rec->differences == -1) {
              continue;
            }
            KeyValuePair outputKVPair;
            outputKVPair.setKey(
              reinterpret_cast<uint8_t*>(&queryID), sizeof(queryID));
            byte* value = fullalignment.toBytes();
            outputKVPair.setValue(
              static_cast<uint8_t*>(value), fullalignment.outputSize);
            writer.write(outputKVPair);
          }
        }
      }
    }
  }
}

AlignmentRecord* CloudBurstReduceFunction::extend(
  const MerRecord& qrytuple, const MerRecord& reftuple) {
  int32_t refStart    = reftuple.offset;
  int32_t refEnd      = reftuple.offset + seedLength;
  int32_t differences = 0;

  if (qrytuple.leftFlankLength != 0) {
    // at least 1 read base on the left needs to be aligned
    int32_t realleftflanklen =
      dnaStringObj.dnaArrLen(qrytuple.leftFlank, qrytuple.leftFlankLength);
    // aligned the pre-reversed strings!
    AlignInfo& a = landauVishkinObj.extend(
      reftuple.leftFlank, reftuple.leftFlankLength, qrytuple.leftFlank,
      qrytuple.leftFlankLength, maxAlignDiff, allowDifferences);

    if (a.alignlen == -1) {
      return &noalignment;
    }  // alignment failed
    if (!a.isBazeaYatesSeed(realleftflanklen, seedLength)) {
      return &noalignment;
    }
    refStart -= a.alignlen;
    differences = a.differences;
  }
  if (qrytuple.rightFlankLength != 0) {
    AlignInfo& b = landauVishkinObj.extend(
      reftuple.rightFlank, reftuple.rightFlankLength, qrytuple.rightFlank,
      qrytuple.rightFlankLength, maxAlignDiff - differences, allowDifferences);

    if (b.alignlen == -1) {
      return &noalignment;
    }  // alignment failed
    refEnd += b.alignlen;
    differences += b.differences;
  }
  fullalignment.refID = reftuple.id;
  fullalignment.refStart = refStart;
  fullalignment.refEnd = refEnd;
  fullalignment.differences = differences;
  fullalignment.isRC = qrytuple.isRC;
  return &fullalignment;
}
