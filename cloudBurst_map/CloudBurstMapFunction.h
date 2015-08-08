#ifndef MAPRED_CLOUD_BURST_MAP_FUNCTION_H
#define MAPRED_CLOUD_BURST_MAP_FUNCTION_H

#include "DNAString.h"
#include "FastaRecord.h"
#include "MerRecord.h"
#include "mapreduce/functions/map/MapFunction.h"

/**
   CloudBurst map function and associated helper classes based on the
   open-source Java implementation of CloudBurst by Michael Schatz.
   See http://sourceforge.net/p/cloudburst-bio/wiki/CloudBurst/
 */
class CloudBurstMapFunction : public MapFunction {
public:
  CloudBurstMapFunction(
    uint32_t _maxAlignDiff, uint32_t _redundancy, int32_t _minReadLen,
    int32_t _maxReadLen);
private:
  uint32_t chunkOverlap;
  uint32_t flankLen;
  uint32_t maxAlignDiff;
  uint32_t maxReadLen;
  uint32_t minReadLen;
  uint32_t  redundancy;
  int32_t seedLen;
  unsigned char* seedBuffer;
  bool isRef;
  std::string refPath;
  void map(KeyValuePair& kvPair, KVPairWriterInterface& writer);
  void configure(KVPairBuffer* buffer);
};

#endif  // MAPRED_CLOUD_BURST_MAP_FUNCTION_H
