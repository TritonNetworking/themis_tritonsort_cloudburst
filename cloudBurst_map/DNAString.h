#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include<stdio.h>
#include<stdint.h>
#include<string>

typedef uint8_t byte;
class DNAString {
public:
  const byte dnaA;
  const byte dnaC;
  const byte dnaG;
  const byte dnaT;
  const byte dnaN;
  const byte space;
  const byte hardstop;
  byte letterToDNA[256];
  byte letterToSeed[256];
  byte dnaToLetter[256];
  byte seedToLetter[256];
  byte rcLetter[256];
  byte* nostr;
  DNAString();
  virtual ~DNAString();
  void initializeRC(byte* retval);
  byte rc(byte letter);
  void reverseComplementSequenceInPlace(byte* arr, int32_t len);
  void initializeLetterToDNA(byte* retval);
  void initializeLetterToSeed(byte* retval);
  void initializeSeedToLetter(byte* retval);
  void initializeDNAToLetter(byte* retval);
  int32_t dnaArrLen(const byte* arr, int32_t len);
  byte byteToDNA(byte letter);
  byte byteToSeed(byte letter);
  byte seedToByte(int32_t seed);
  byte dnaToByte(byte dna);
  int32_t arrToDNALen(int32_t len);
  int32_t arrToSeedLen(int32_t len, int32_t REDUNDANCY);
  int32_t arrToSeed(
    byte* arr, int32_t arrpos, int32_t len, byte* seed, int32_t seedpos,
    int32_t id, int32_t REDUNDANCY, int32_t ISQRY);
  bool repSeed(byte* seq, int32_t start, int32_t SEED_LEN);
  byte* seedToArr(byte* seed, int32_t SEEDLEN, int32_t REDUNDANCY);
  int32_t arrToDNAStr(
    byte* arr, int32_t arrpos, int32_t len, byte* out, int32_t outpos);
  int32_t arrToDNAStrRev(
    byte* arr, int32_t arrstart, int32_t len, byte* out, int32_t outpos);
  byte* arrToDNA(byte* arr, int32_t start, int32_t len);
  byte* arrToDNA(byte* arr, int32_t len);
  byte* dnaToArr(const byte* dna, int32_t dnapos, int32_t dnalen);
  byte* dnaToArr(const byte* dna, int32_t len);
  int32_t dnaToArrLen(const byte* dna, int32_t dnapos, int32_t dnalen);
  byte* stringToBytes(std::string src);
  bool arrHasN(byte* seq, int32_t start, int32_t len);
private:
};
#endif  //  DNA_STRING_H
