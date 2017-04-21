#include "timestamp.hpp"
#include <stdio.h>

void PrintTimeStamp(std::string label, time_t *t) {
  time_t tt;
  time(&tt);
  if(t != NULL)
    *t = tt;
  printf("[%10s] %s", label.c_str(), ctime(&tt));
}
