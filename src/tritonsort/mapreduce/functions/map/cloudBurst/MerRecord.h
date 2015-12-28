#ifndef _MER_RECORD_H_
#define _MER_RECORD_H_

#include"DNAString.h"

typedef uint8_t byte;
class MerRecord {
public:
  MerRecord();
  ~MerRecord();
  MerRecord(byte* t, int32_t len);
  int32_t toBytesLen(int32_t leftlen, int32_t rightlen);
  byte* toBytes(byte* seq, int32_t leftstart, int32_t leftlen,
      int32_t rightstart, int32_t rightlen);
  void fromBytes(const byte* bytes, int32_t length);
  byte* toBytes(int32_t id);

  /// \todo(AR) These fields should be private or const
  int32_t offsetIndex;
  int32_t idIndex;
  bool isReference;
  bool isRC;
  int32_t offset;
  int32_t id;
  const uint8_t* leftFlank;
  uint32_t leftFlankLength;
  const uint8_t* rightFlank;
  uint32_t rightFlankLength;

private:
  byte* sbuffer;
};
#endif  // _MER_RECORD_H
