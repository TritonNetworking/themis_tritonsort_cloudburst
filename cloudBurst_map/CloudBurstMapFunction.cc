#include <string>
#include "CloudBurstMapFunction.h"

typedef uint8_t byte;

//constructor
CloudBurstMapFunction::CloudBurstMapFunction(
  uint32_t _maxAlignDiff, uint32_t _redundancy, int32_t _minReadLen,
  int32_t _maxReadLen)
  : chunkOverlap(1024),
    maxAlignDiff(_maxAlignDiff),
    maxReadLen(_maxReadLen),
    minReadLen(_minReadLen),
    redundancy(_redundancy) {
    //calculate seedLength based on min read len and K
  seedLen = minReadLen/(maxAlignDiff + 1);
  flankLen = maxReadLen - seedLen + maxAlignDiff;
}

// Get source of buffer to figure out whether
// KV Pair is from reference set or Query Set
void CloudBurstMapFunction::configure(KVPairBuffer* buffer) {
  const std::string &fileName = buffer->getSourceName();
  ABORT_IF(fileName.empty(),
           "CloudBurst requires a FileByteStreamConverter to set filenames on "
           "map input buffers");
  isRef = fileName.find("ref", 0) != std::string::npos;
}

void CloudBurstMapFunction::map(
  KeyValuePair& kvPair, KVPairWriterInterface& writer) {

  MerRecord seedInfo;
  DNAString dnaStringObj;
  byte* merInfo;
  if (redundancy > 1) {
    seedBuffer = new byte[(seedLen + 3)/4 + 2];
  } else {
    seedBuffer = new byte[(seedLen + 3)/4 + 1];
  }
  // DNA sequences stores as KV Pair (id, Seqinfo)
  // Seqinfo is tuple (sequence, start_offset)
  // Map function emits KV pairs(seed, MerInfo)
  // seed is a sequence of length S
  // Merinfo is tuple (id, position, isRef, isRc, left_flank, right_flank)
  FastaRecord fastaRecord(kvPair.getValue(), kvPair.getValueLength());
  byte* seq = fastaRecord.sequence;

  int32_t realOffsetStart = fastaRecord.offset;
  bool isLast = fastaRecord.lastChunk;
  int32_t seqLen = fastaRecord.sequenceLength;
  memcpy(&seedInfo.id, kvPair.getKey(), sizeof(seedInfo.id));

  seedInfo.isReference = isRef;
  seedInfo.isRC = false;

  if (isRef) {
    // Sequence is chunk of the reference
    int32_t startOffset = 0;

    // If I'm not the first chunk, shift over so there is room for left flank
    if (realOffsetStart != 0) {
      startOffset = chunkOverlap + 1 - flankLen - seedLen;
      realOffsetStart += startOffset;
    }
    // stop so the last mer will just fit
    int32_t end = seqLen - seedLen + 1;
    // If I'm not the last chunk, stop so the right flank will fit as well
    if (!isLast) {
      end -= flankLen;
    }
    // emit the mers starting at every position in the range
    for (int32_t start = startOffset, realOffset = realOffsetStart;
      start < end; start++, realOffset++) {
      // don't bother with seeds with N's
      // DNA is expressed as combination of A,D,G,C
      if (dnaStringObj.arrHasN(seq, start, seedLen)) {
        continue;
      }
      seedInfo.offset = realOffset;
      // figure out the ranges for the flanking sequence
      int32_t leftStart = start - flankLen;
      if (leftStart < 0) {
        leftStart = 0;
      }
      int32_t leftLen = start - leftStart;
      int32_t rightStart = start + seedLen;
      int32_t rightEnd = rightStart + flankLen;
      if (rightEnd > seqLen) {
        rightEnd = seqLen;
      }
      int32_t rightLen = rightEnd - rightStart;
      KeyValuePair outputKVPair;
      merInfo = seedInfo.toBytes(seq, leftStart, leftLen, rightStart, rightLen);
      int32_t outputLen = seedInfo.toBytesLen(leftLen, rightLen);
      outputKVPair.setValue(static_cast<byte*>(merInfo), outputLen);
      if ((redundancy >1) && (dnaStringObj.repSeed(seq, start, seedLen))) {
        for (uint32_t r = 0; r < redundancy; r++) {
          int32_t length = dnaStringObj.arrToSeed
            (seq, start, seedLen, seedBuffer, 0 , r , redundancy, 0);
          outputKVPair.setKey(seedBuffer, length);
          writer.write(outputKVPair);
        }
      } else {
        int32_t length = dnaStringObj.arrToSeed(seq,
          start, seedLen, seedBuffer, 0, 0, redundancy, 0);
        outputKVPair.setKey(seedBuffer, length);
        writer.write(outputKVPair);
      }
      delete[] merInfo;
    }  // END OF FOR
    delete[] seq;
    delete[] seedBuffer;
  } else {
    // Skip reads that can't possibly align end-to-end with<= K differences
    uint32_t numN = 0;
    for (int32_t i = 0; i < seqLen; i++) {
      if (seq[i] == 'N') {
        numN++;
      }
    }
    if (numN > maxAlignDiff) {
      return;
    }
    for (int32_t rc = 0; rc < 2 ; rc++) {
      if (rc == 1) {
        // reverse complement the sequence
        dnaStringObj.reverseComplementSequenceInPlace(seq, seqLen);
        seedInfo.isRC = true;
      }
      // only emit the non-overlapping mers
      for (int32_t i = 0; i + seedLen <= seqLen; i += seedLen) {
        int32_t len;
        if (dnaStringObj.arrHasN(seq, i, seedLen)) {
          continue;
        }
        if ((redundancy > 1) && (dnaStringObj.repSeed(seq, i, seedLen))) {
          len = dnaStringObj.arrToSeed(seq, i,
            seedLen, seedBuffer, 0, seedInfo.id, redundancy, 1);
        } else {
          len = dnaStringObj.arrToSeed(seq, i,
            seedLen, seedBuffer, 0, 0, redundancy, 1);
        }
        seedInfo.offset = i;
        // figure out the ranges for the flanking sequence
        int32_t leftStart = 0;
        int32_t leftLen = i;
        int32_t rightStart = i + seedLen;
        int32_t rightLen = seqLen - rightStart;
        KeyValuePair outputKVPair;
        merInfo = seedInfo.toBytes(
          seq, leftStart, leftLen, rightStart, rightLen);
        int32_t outputLen = seedInfo.toBytesLen(leftLen, rightLen);
        outputKVPair.setKey(seedBuffer, len);
        outputKVPair.setValue(static_cast<uint8_t*>(merInfo), outputLen);
        writer.write(outputKVPair);
        delete[] merInfo;
      }
    }
    delete[] seq;
    delete[] seedBuffer;
  }
}
