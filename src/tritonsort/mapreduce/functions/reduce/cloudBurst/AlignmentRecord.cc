#include "AlignmentRecord.h"
#include "core/MemoryUtils.h"

// size 17 comes from the total size of the output format
// which if of form number, number, number, number, flag
AlignmentRecord::AlignmentRecord()
  : outputSize(17),
    refIDIndex(1),
    refStartIndex(5),
    refEndIndex(9),
    refDifferencesIndex(13) {
  sbuffer = new (themis::memcheck) byte[outputSize];
}

AlignmentRecord::~AlignmentRecord() {
  delete[] sbuffer;
}

AlignmentRecord& AlignmentRecord::operator= (const AlignmentRecord& that) {
  if (this != &that) {
    if (that.sbuffer != NULL) {
      memcpy(sbuffer, that.sbuffer, outputSize);
    }
    refID = that.refID;
    refStart = that.refStart;
    refEnd = that.refEnd;
    differences = that.differences;
    isRC = that.isRC;
  }
  return *this;
}

// size 17 comes from the total size of the output format
AlignmentRecord::AlignmentRecord(
  int32_t _refID, int32_t _refStart, int32_t _refEnd, int32_t _differences,
  bool _isRC)
  : outputSize(17),
    refIDIndex(1),
    refStartIndex(5),
    refEndIndex(9),
    refDifferencesIndex(13),
    refID(_refID),
    refStart(_refStart),
    refEnd(_refEnd),
    differences(_differences),
    isRC(_isRC) {
  sbuffer = new (themis::memcheck) byte[outputSize];
}


// size 17 comes from the total size of the output format
AlignmentRecord::AlignmentRecord(AlignmentRecord& other)
  : outputSize(17),
    refIDIndex(1),
    refStartIndex(5),
    refEndIndex(9),
    refDifferencesIndex(13) {
  refID = other.refID;
  refStart = other.refStart;
  refEnd = other.refEnd;
  differences = other.differences;
  isRC = other.isRC;
}


// size 17 comes from the total size of the output format
AlignmentRecord::AlignmentRecord(byte* b)
  : outputSize(17),
    refIDIndex(1),
    refStartIndex(5),
    refEndIndex(9),
    refDifferencesIndex(13) {
  fromBytes(b);
}

void AlignmentRecord::set(AlignmentRecord other) {
  refID = other.refID;
  refStart = other.refStart;
  refEnd = other.refEnd;
  differences = other.differences;
  isRC = other.isRC;
}

byte* AlignmentRecord:: toBytes() {
  sbuffer[0] = (byte) (isRC ? 1 : 0);
  *(reinterpret_cast<int32_t*>(sbuffer + refIDIndex)) = refID;
  *(reinterpret_cast<int32_t*>(sbuffer + refStartIndex)) = refStart;
  *(reinterpret_cast<int32_t*>(sbuffer + refEndIndex)) = refEnd;
  *(reinterpret_cast<int32_t*>(sbuffer + refDifferencesIndex)) = differences;
  return sbuffer;
}

void AlignmentRecord:: fromBytes(byte* raw) {
  isRC = (raw[0] == 1);
  memcpy(&refID, &raw[refIDIndex], sizeof(refID));
  memcpy(&refStart, &raw[refStartIndex], sizeof(refStart));
  memcpy(&refEnd, &raw[refEndIndex], sizeof(refEnd));
  memcpy(&differences, &raw[refDifferencesIndex], sizeof(differences));
}
