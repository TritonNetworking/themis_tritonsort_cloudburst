#include<assert.h>
#include<stdio.h>
#include<stdint.h>
#include<string.h>

#include "DNAString.h"
#include "core/MemoryUtils.h"

typedef uint8_t byte;
// cconstructor
DNAString::DNAString()
  : dnaA(0x00),
    dnaC(0x01),
    dnaG(0x02),
    dnaT(0x04),
    dnaN(0x08),
    space(0x0F),
    hardstop(0xFF) {
  initializeLetterToDNA(letterToDNA);
  initializeLetterToSeed(letterToSeed);
  initializeDNAToLetter(dnaToLetter);
  initializeSeedToLetter(seedToLetter);
  initializeRC(rcLetter);
}

DNAString::~DNAString() {
}

// initialize
void DNAString::initializeRC(byte* retval) {
  memset(retval, 'N' , 256);
  retval['A'] = 'T';
  retval['a'] = 't';
  retval['T'] = 'A';
  retval['t'] = 'a';
  retval['C'] = 'G';
  retval['c'] = 'g';
  retval['G'] = 'C';
  retval['g'] = 'c';
}

//initialize
byte DNAString::rc(byte letter) {
  return rcLetter[letter];
}

//reverseInPlace
void DNAString::reverseComplementSequenceInPlace(byte* arr, int32_t len) {
  int32_t i = 0, j = len - 1;
  for (; i < j; i++, j--) {
    byte t = arr[i];
    arr[i] = rcLetter[arr[j]];
    arr[j] = rcLetter[t];
  }
  if (i == j) {
    arr[i] = rcLetter[arr[i]];
  }
}

void DNAString::initializeLetterToDNA(byte* retval) {
  memset(retval, dnaN, 256);
  retval['A'] = dnaA;
  retval['a'] = dnaA;
  retval['C'] = dnaC;
  retval['c'] = dnaC;
  retval['G'] = dnaG;
  retval['g'] = dnaG;
  retval['T'] = dnaT;
  retval['t'] = dnaT;
  retval[' '] = space;
  retval[';'] = hardstop;
}

// Seeds don't have N's so use 2 bits / bp
void DNAString::initializeLetterToSeed(byte* retval) {
  memset(retval,  0, 256);
  retval['A'] = 0;
  retval['a'] = 0;
  retval['C'] = 1;
  retval['c'] = 1;
  retval['G'] = 2;
  retval['g'] = 2;
  retval['T'] = 3;
  retval['t'] = 3;
}

// Seeds don't have N's so use 2 bits / bp
void DNAString::initializeSeedToLetter(byte* retval) {
  memset(retval, 'N', 256);
  retval[0] = 'A';
  retval[1] = 'C';
  retval[2] = 'G';
  retval[3] = 'T';
}

void DNAString::initializeDNAToLetter(byte* retval) {
  memset(retval, 'N', 256);
  retval[dnaA]    = 'A';
  retval[dnaC]    = 'C';
  retval[dnaG]    = 'G';
  retval[dnaT]    = 'T';
  retval[dnaN]    = 'N';
  retval[space]    = ' ';
  retval[hardstop] = ';';
}

int32_t DNAString::dnaArrLen(const byte* arr, int32_t len) {
  int32_t retval = len * 2;
  if (len > 0) {
    if ((arr[len-1] & 0x0F) == space) {
      retval--;
    }
  }
  return retval;
}

byte DNAString::byteToDNA(byte letter) {
  return letterToDNA[letter];
}

byte DNAString::byteToSeed(byte letter) {
  return letterToSeed[letter];
}

byte DNAString::seedToByte(int32_t seed) {
  return seedToLetter[seed & 0x03];
}

byte DNAString::dnaToByte(byte dna) {
  return dnaToLetter[dna];
}

int32_t DNAString::arrToDNALen(int32_t len) {
  return (len+1)/2;
}

//  arrToSeedLen
//  Calculate seed Length
int32_t DNAString::arrToSeedLen(int32_t len, int32_t REDUNDANCY) {
  // 4 bytes per basepair + reserve 1 byte for ref/qry flag
  if (REDUNDANCY > 1) {
    // use an extra byte for redundancy
    return ((len+3)/4) + 1 + 1;
  }
  return (len+3)/4 + 1;
}

//  arrTOseed
//  change sequence to seed type for next phase
int32_t DNAString::arrToSeed(
  byte* arr, int32_t arrPos, int32_t len, byte* seed, int32_t seedPos,
  int32_t id, int32_t REDUNDANCY, int32_t ISQRY) {

  int32_t seedlen = (len+3)/4+1;
  int32_t arrend = arrPos + len;
  while (arrPos+3 < arrend) {
    seed[seedPos] = (byte) (((byteToSeed(arr[arrPos])   << 6) |
          (byteToSeed(arr[arrPos+1]) << 4) |
          (byteToSeed(arr[arrPos+2]) << 2) |
          (byteToSeed(arr[arrPos+3]))));

    seedPos++;
    arrPos+=4;
  }

  int32_t remaining = arrend-arrPos;
  if (remaining == 3) {
    seed[seedPos] = (byte) (((byteToSeed(arr[arrPos])   << 6) |
          (byteToSeed(arr[arrPos+1]) << 4) |
          (byteToSeed(arr[arrPos+2]) << 2)));
    seedPos++;
  } else if (remaining == 2) {
    seed[seedPos] = (byte) (((byteToSeed(arr[arrPos])  << 6) |
          (byteToSeed(arr[arrPos+1]) << 4)));
    seedPos++;
  } else if (remaining == 1) {
    seed[seedPos] = (byte) (((byteToSeed(arr[arrPos])   << 6)));
    seedPos++;
  }

  if (REDUNDANCY > 1) {
    seed[seedPos] = (byte) ((id % REDUNDANCY) & 0xFF);
    seedPos++;
    seedlen++;
  }
  seed[seedPos] = (byte) ISQRY;
  return seedlen;
}

//  check if string is bunch of repeats
bool DNAString::repSeed(byte* sequence, int32_t start, int32_t SEED_LEN) {
  byte first = sequence[start];

  for (int32_t i = 1; i < SEED_LEN; i++) {
    if (sequence[i+start] != first) {
      return false;
    }
  }
  return true;
}

//  seedToDNAStr
byte* DNAString::seedToArr(byte* seed, int32_t SEEDLEN, int32_t REDUNDANCY) {
  byte* retval = new (themis::memcheck) byte[SEEDLEN];
  int32_t outPos = 0;
  int32_t seedPos = 0;
  int32_t slen = arrToSeedLen(SEEDLEN, 0);

  // the very last byte has the isref flag
  // otherwise all but the last byte will have 4 bp
  while (seedPos < slen - 2) {
    retval[outPos]   = seedToByte(seed[seedPos] >> 6);
    retval[outPos+1] = seedToByte(seed[seedPos] >> 4);
    retval[outPos+2] = seedToByte(seed[seedPos] >> 2);
    retval[outPos+3] = seedToByte(seed[seedPos]);
    outPos += 4;
    seedPos++;
  }
  int32_t diff = SEEDLEN - outPos;
  if (diff == 4) {
    retval[outPos]   = seedToByte(seed[seedPos] >> 6);
    retval[outPos + 1] = seedToByte(seed[seedPos] >> 4);
    retval[outPos + 2] = seedToByte(seed[seedPos] >> 2);
    retval[outPos + 3] = seedToByte(seed[seedPos]);
  } else if (diff == 3) {
    retval[outPos]   = seedToByte(seed[seedPos] >> 6);
    retval[outPos + 1] = seedToByte(seed[seedPos] >> 4);
    retval[outPos + 2] = seedToByte(seed[seedPos] >> 2);
  } else if (diff == 2) {
    retval[outPos]   = seedToByte(seed[seedPos] >> 6);
    retval[outPos + 1] = seedToByte(seed[seedPos] >> 4);
  } else if (diff == 1) {
    retval[outPos]   = seedToByte(seed[seedPos] >> 6);
  }

  return retval;
}

//  arrToDNAStr
int32_t DNAString::arrToDNAStr(
  byte* arr, int32_t arrPos, int32_t len, byte* out, int32_t outPos) {
  int32_t dnaLen = (len+1)/2;
  int32_t arrend = arrPos + len;
  while (arrPos+1 < arrend) {
    out[outPos] = (byte) ((byteToDNA(arr[arrPos]) << 4)
        | (byteToDNA(arr[arrPos+1])));
    outPos++;
    arrPos+=2;
  }

  if (arrPos < arrend) {
    out[outPos] = (byte) ((byteToDNA(arr[arrPos]) << 4) | (space));
  }

  return dnaLen;
}

//  arrToDNAStrRev
//  Sequence length to dna and reverse
int32_t DNAString::arrToDNAStrRev(
  byte* arr, int32_t arrStart, int32_t len, byte* out, int32_t outPos) {
  //  calculate dnaLen from the sequence length (1 byte for every 2bytes)
  int32_t dnaLen = (len+1)/2;
  int32_t arrPos = arrStart + len - 1;

  while (arrPos > arrStart) {
    out[outPos] = (byte) ((byteToDNA(arr[arrPos]) << 4)
        | (byteToDNA(arr[arrPos - 1])));
    outPos++;
    arrPos -= 2;
  }
  if (arrPos == arrStart) {
    out[outPos] = (byte) ((byteToDNA(arr[arrPos]) << 4) | (space));
  }
  return dnaLen;
}


byte* DNAString::arrToDNA(byte* arr, int32_t start, int32_t len) {
  int32_t dnaLen = (len + 1)/2;
  byte* dna = new (themis::memcheck) byte[dnaLen];
  arrToDNAStr(arr, start, len, dna, 0);
  return dna;
}

byte* DNAString::arrToDNA(byte* arr, int32_t len) {
  return arrToDNA(arr, 0, len);
}

byte* DNAString::dnaToArr(const byte* dna, int32_t dnaPos, int32_t dnaLen) {
  if (dnaLen == 0) {
    return new (themis::memcheck) byte[1];
  }
  int32_t arrlen = dnaLen*2;
  int32_t arrPos = 0;
  if ((dna[dnaPos + dnaLen - 1] & 0x0F) == space) {
    arrlen--;
  }
  byte* arr = new (themis::memcheck) byte[arrlen];
  int32_t dnaend = dnaPos + dnaLen - 1;
  // don't check the last byteacter, may be a 'space'
  while (dnaPos < dnaend) {
    arr[arrPos]   = dnaToByte((byte)((dna[dnaPos] & 0xF0) >> 4));
    arr[arrPos+1] = dnaToByte((byte)((dna[dnaPos] & 0x0F)));
    dnaPos++;
    arrPos+=2;
  }
  arr[arrPos] = dnaToByte((byte)((dna[dnaPos] & 0xF0) >> 4));
  arrPos++;
  if ((dna[dnaPos] & 0x0F) != space) {
    assert(arrPos < arrlen);
    arr[arrPos] = dnaToByte((byte)((dna[dnaPos] & 0x0F)));
  }
  return arr;
}

byte* DNAString::dnaToArr(const byte* dna, int32_t len) {
  return dnaToArr(dna, 0, len);
}

int32_t DNAString::dnaToArrLen(const byte* dna, int32_t dnaPos,
    int32_t dnaLen) {
  int32_t arrlen = dnaLen*2;
  if ((dna[dnaPos + dnaLen - 1] & 0x0F) == space) {
    arrlen--;
  }
  return arrlen;
}

byte* DNAString::stringToBytes(std::string src) {
  int32_t srcLen = src.length();
  byte* ret = new (themis::memcheck) byte[srcLen];
  for (int32_t i = 0; i < srcLen; i++) {
    ret[i] = (byte) src[i];
  }
  return ret;
}


bool DNAString::arrHasN(byte* sequence, int32_t start, int32_t len) {
  for (int32_t n = start; n < start+len; n++) {
    if (letterToDNA[sequence[n]] == dnaN) {
      return true;
    }
  }
  return false;
}
