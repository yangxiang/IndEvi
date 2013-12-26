#include <assert.h>
#include <limits>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "IndEviStatic.hh"

std::string IndEviStatic::ds_fpath;
Itemset::item_t IndEviStatic::max_item = 0;
Itemset::item_t IndEviStatic::min_item = 0;
Transactionset::transaction_t IndEviStatic::num_transactions = 0;
double IndEviStatic::beta = -1.0;
unsigned int IndEviStatic::merge_k = std::numeric_limits<unsigned int>::max();
int64_t IndEviStatic::start_usec = 0;

int64_t t_usec()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (int64_t)tv.tv_sec * pow(10, 6) + tv.tv_usec;
}

double IndEviStatic::elapsed_time_sec()
{
  return (double)(t_usec() - IndEviStatic::start_usec) / pow(10, 6);
}

void IndEviStatic::init_timer()
{
  IndEviStatic::start_usec = t_usec();
}
