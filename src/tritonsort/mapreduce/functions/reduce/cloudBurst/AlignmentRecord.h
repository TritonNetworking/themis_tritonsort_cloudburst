#ifndef _ALIGNMENT_RECORD_H_
#define _ALIGNMENT_RECORD_H_
#include<stdint.h>
#include<string.h>

typedef uint8_t byte;

class AlignmentRecord {
public:
  byte* sbuffer;
  const int32_t outputSize;
  const int32_t refIDIndex;
  const int32_t refStartIndex;
  const int32_t refEndIndex;
  const int32_t refDifferencesIndex;
  int32_t refID;
  int32_t refStart;
  int32_t refEnd;
  int32_t differences;
  bool isRC;
  AlignmentRecord();
  AlignmentRecord(AlignmentRecord& other);
  AlignmentRecord(byte* b);
  virtual ~AlignmentRecord();
  AlignmentRecord& operator= (const AlignmentRecord& that);
  AlignmentRecord(
    int32_t refid, int32_t refstart, int32_t refend,
    int32_t differences, bool rc);
  void set(AlignmentRecord other);
  byte* toBytes();
  void fromBytes(byte* raw);
private:
};
#endif  //  _ALIGNMENT_RECORD_H
