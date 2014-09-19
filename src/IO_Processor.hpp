#ifndef __IO_Processor__
#define __IO_Processor__

#include <iostream>
#include <vector>
#include "sequenceUtil.hpp"

class IO_Processor {
  
public:
  
  IO_Processor(){};
  void initialize(char *mapping, long length, unsigned long long bound, unsigned long long &log);
  vector <kmer_pair> kmer_buffer;
  
private:
  char *map;
  long length;
  unsigned long long bound_kmer;
  void generate_kmer_buffer(unsigned long long &log);   
};

#endif