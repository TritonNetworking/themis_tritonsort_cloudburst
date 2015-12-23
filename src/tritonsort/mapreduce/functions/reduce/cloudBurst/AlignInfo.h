#ifndef _ALIGN_INFO_H
#define _ALIGN_INFO_H

class AlignInfo {
public:
  int32_t alignlen;
  int32_t differences;
  AlignInfo();
  AlignInfo(
    int32_t len, int32_t k, int32_t* pdist, int32_t* pweightMatrix,
    int32_t dlen);
  void setVals(
    int32_t len, int32_t k, int32_t* pdist, int32_t* pweightMatrix,
    int32_t dlen);
  // ------------------------- isBazeaYatesSeed --------------------------
  // Since an alignment may be recompute k+1 times for each of the k+1 seeds,
  // see if the current alignment is the leftmost alignment by checking for
  // differences in the proceeding chunks of the query
  bool isBazeaYatesSeed(int32_t qlen, int32_t kmerlen);
private:
  int32_t* dist;
  int32_t* weightMatrix;
  int32_t distlen;
};
#endif  //  _ALIGN_INFO_H
