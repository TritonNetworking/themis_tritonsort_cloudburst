#include<arpa/inet.h>
#include<cstdlib>
#include<stdint.h>
#include<string.h>

#include "FastaRecord.h"
#include "core/MemoryUtils.h"

typedef uint8_t byte;

FastaRecord::FastaRecord()
  : dnaStartPos(5),
    sequence(NULL),
    lastChunk(false),
    offset(0),
    sequenceLength(0) {
}

FastaRecord::~FastaRecord() {
}

//  Fasta file reader.
//  DNA sequences are available in fasta format
//  Input file is compressed in byte format and read back in ascii
//  chars
FastaRecord::FastaRecord(const uint8_t* rawByteString, int32_t valueLength)
  : dnaStartPos(5) {
  lastChunk = rawByteString[0];
  memcpy(&offset, &rawByteString[1], sizeof(offset));
  //  we need htonl here because of ordering assumption made in input file
  offset = htonl(offset);
  // sequence is double the string size
  // pos 5 marks begining of dna string
  // offset (4) + lastChunk (1)
  sequence = dnaStringObj.dnaToArr(
    rawByteString, dnaStartPos, valueLength - dnaStartPos);
  sequenceLength = dnaStringObj.dnaToArrLen(
    rawByteString, dnaStartPos, valueLength - dnaStartPos);
}

// Convert FastaFormat in byte format
byte* FastaRecord::toBytes(int32_t arrLen) {
  DNAString dnaStringObj;
  byte* dna = dnaStringObj.arrToDNA(sequence, arrLen);
  int32_t dnaLen = (arrLen + 1)/2;
  int32_t len = dnaStartPos +  dnaLen;
  byte* buffer = new (themis::memcheck) byte[len];
  buffer[0] = (byte) (lastChunk ? 1 : 0);
  memcpy(buffer, &offset, sizeof(offset));
  memcpy(buffer + dnaStartPos, dna, dnaLen);
  delete[] dna;
  return buffer;
}
