#ifndef _FASTA_RECORD_H_
#define _FASTA_RECORD_H_

#include"DNAString.h"
typedef uint8_t byte;

class FastaRecord {
public:
  FastaRecord();

  FastaRecord(const uint8_t* rawByteString, int32_t valueLength);

  virtual ~FastaRecord();

  byte* toBytes(int32_t arrLen);

  ///\todo(AR) All these non-const public members make my brain bleed
  const int dnaStartPos;
  byte* sequence;
  bool lastChunk;
  int32_t offset;
  int32_t sequenceLength;
  DNAString dnaStringObj;
};
#endif  // _FASTA_RECORD_H
