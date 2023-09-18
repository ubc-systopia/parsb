
#ifndef GAPBS_SB_UTIL_H
#define GAPBS_SB_UTIL_H

#include <cinttypes>
#include <stdio.h>
#include <string>

#include "timer.h"


/*
GAP Benchmark Suite
Author: Scott Beamer

Miscellaneous helpers that don't fit into classes

Adapted by atrostan on 13/10/22.

*/


static const int64_t kRandSeed = 27491095;


inline void PrintLabel(const std::string &label, const std::string &val) {
  printf("%-21s%7s\n", (label + ":").c_str(), val.c_str());
}

inline void PrintTime(const std::string &s, double seconds) {
  printf("%-21s%3.5lf\n", (s + ":").c_str(), seconds);
}

inline void PrintStep(const std::string &s, int64_t count) {
  printf("%-14s%14" PRId64 "\n", (s + ":").c_str(), count);
}

inline void PrintStep(int step, double seconds, int64_t count = -1) {
  if (count != -1)
    printf("%5d%11" PRId64 "  %10.5lf\n", step, count, seconds);
  else
    printf("%5d%23.5lf\n", step, seconds);
}

inline void PrintStep(const std::string &s, double seconds, int64_t count = -1) {
  if (count != -1)
    printf("%5s%11" PRId64 "  %10.5lf\n", s.c_str(), count, seconds);
  else
    printf("%5s%23.5lf\n", s.c_str(), seconds);
}

// Runs op and prints the time it took to execute labelled by label
#define TIME_PRINT(label, op)       \
  {                                 \
    Timer t_;                       \
    t_.Start();                     \
    (op);                           \
    t_.Stop();                      \
    PrintTime(label, t_.Seconds()); \
  }


template<typename T_>
class RangeIter {
  T_ x_;

  public:
  explicit RangeIter(T_ x) : x_(x) {}
  bool operator!=(RangeIter const &other) const { return x_ != other.x_; }
  T_ const &operator*() const { return x_; }
  RangeIter &operator++() {
    ++x_;
    return *this;
  }
};

template<typename T_>
class Range {
  T_ from_;
  T_ to_;

  public:
  explicit Range(T_ to) : from_(0), to_(to) {}
  Range(T_ from, T_ to) : from_(from), to_(to) {}
  RangeIter<T_> begin() const { return RangeIter<T_>(from_); }
  RangeIter<T_> end() const { return RangeIter<T_>(to_); }
};
#endif//GAPBS_SB_UTIL_H
