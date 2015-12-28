#include<assert.h>
#include<stdint.h>
#include<string.h>

#include "MerRecord.h"
#include "core/MemoryUtils.h"

typedef uint8_t byte;
// constructor
MerRecord::MerRecord()
  : offsetIndex(sizeof(isReference)),
    idIndex(offsetIndex+sizeof(offset)),
    isReference(false),
    isRC(false),
    offset(0),
    id(0),
    leftFlank(NULL),
    rightFlank(NULL),
    sbuffer(NULL) {
}
//  constructor
MerRecord::MerRecord(byte* t, int32_t len)
  : offsetIndex(sizeof(isReference)),
    idIndex(offsetIndex+sizeof(offset)) {
  fromBytes(t, len);
}

MerRecord::~MerRecord() {
}

// Convert MerRecord to Bytes
int32_t MerRecord::toBytesLen(int32_t leftlen, int32_t rightlen) {
  DNAString dnaStringObj;
  // +1 for hardstop between left and right flank
  int32_t len = sizeof(isReference) + sizeof(offset) + sizeof(id) + 1;
  if (leftlen > 0)  {
    len += dnaStringObj.arrToDNALen(leftlen);
  }
  if (rightlen > 0) {
    len += dnaStringObj.arrToDNALen(rightlen);
  }
  return len;
}

// Pack the MerRecord information into a BytesWritable
// extract the flanking sequence on-the-fly do avoid copying as much as possible
byte* MerRecord::toBytes(
  byte* seq, int32_t leftstart, int32_t leftlen, int32_t rightstart,
  int32_t rightlen) {
  DNAString dnaStringObj;
  // +1 for hardstop between left and right flank
  int32_t len = sizeof(isReference) + sizeof(offset) + sizeof(id) + 1;
  if (leftlen > 0)  {
    len += dnaStringObj.arrToDNALen(leftlen);
  }
  if (rightlen > 0) {
    len += dnaStringObj.arrToDNALen(rightlen);
  }
  byte* sbuffer = new (themis::memcheck) byte[len];
  // byte is used to indicate whether the given sequence is from
  // reference string and whether the sequence is replicated
  // 0x11 - indicate  reference and rc set
  // 0x00 - indicate rederence and rc not set
  sbuffer[0] = (byte) ((isReference ? 0x01 : 0x00) |
    (isRC ? 0x10 : 0x00));
  memcpy(&sbuffer[offsetIndex], &offset, sizeof(offset));
  memcpy(&sbuffer[idIndex], &id, sizeof(id));
  int32_t pos = 9;
  if (leftlen > 0) {
    pos += dnaStringObj.arrToDNAStrRev(seq, leftstart, leftlen, sbuffer, pos);
  }
  sbuffer[pos] = dnaStringObj.hardstop;
  pos++;
  if (rightlen > 0) {
    pos += dnaStringObj.arrToDNAStr(seq, rightstart, rightlen, sbuffer, pos);
  }
  return sbuffer;
}

//  change MerRecord tobytes
byte* MerRecord::toBytes(int32_t id) {
  byte* buffer = new (themis::memcheck) byte[4];
  memcpy(buffer, &id, sizeof(id));
  return buffer;
}

// Unpack the raw bytes and set the MerRecord fields
void MerRecord::fromBytes(const byte* bytes, int32_t length) {
  DNAString dnaStringObj;

  // byte is used to indicate whether the given sequence is from
  // reference string and whether the sequence is replicated
  // 0x11 - indicate  reference and rc set
  // 0x00 - indicate reference and rc not set
  isReference = (bytes[0] & 0x01) == 0x01;
  isRC        = (bytes[0] & 0x10) == 0x10;
  offset = *reinterpret_cast<const int32_t*>(bytes + offsetIndex);
  id = *reinterpret_cast<const int32_t*>(bytes + idIndex);

  // Flanks are formatted as left;right, after all of the header information.
  uint32_t flankStart = sizeof(isReference) + sizeof(offset) + sizeof(id);

  // Left flank begins immediately after headers.
  leftFlank = bytes + flankStart;
  leftFlankLength = 0;
  for (int32_t i = flankStart; i < length; ++i) {
    if (bytes[i] == dnaStringObj.hardstop) {
      leftFlankLength = i - flankStart;
      // The right flank will start after the hardstop.
      flankStart = i + 1;
      break;
    }
  }
  // Right flank begins after the hardstop, or is the entire portion after the
  // header if no hardstop was found.
  rightFlank = bytes + flankStart;
  rightFlankLength = length - flankStart;
}
