#include "pairhmm_common.h"

#ifdef __APPLE__
#define WEAK __attribute__((weak_import))
#else
#define WEAK __attribute__((weak))
#endif

#define MAX_READ_LENGTH 32*1024
#define MAX_HAP_LENGTH 32*1024
#define MAX_NUM_RESULTS 64*1024*1024

#define ROWS 26
#define COLS 8

typedef struct {
  int length;
  const char* bases;
  const char* q;
  const char* i;
  const char* d;
  const char* c;
} Read;

typedef struct {
  int length;
  const char* bases;
} Haplotype;

typedef struct {
  int num_reads;
  int num_haps;
  long num_cells;
  Read* reads;
  Haplotype* haps;
  float* results;
} Batch;

typedef struct {
  char base;
  char position;
  short hap_num;
  float y_init;
} HapData;

typedef struct {
  char base;
  char position;
  short read_num;
  float mx;
  float my;
  float gg;
  float mm_1m_qual;
  float mm_qual_div3;
  float gm_1m_qual;
  float gm_qual_div3;
} ReadData;

typedef struct {
  int result_read_num;
  int result_hap_num;
  float result;
} ResultData;

float fpga_pairhmm(testcase testcase);
