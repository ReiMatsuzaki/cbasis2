#ifndef TIMER_HPP
#define TIMER_HPP

#include <map>
#include <string>

namespace
{
  using std::map;
  using std::string;
  using std::cout;
  using std::endl;
}

class Timer {

  enum Status {
    kNone,
    kStart,
    kEnd,
  };  

  struct TimeData {
    int id;
    clock_t t0;
    clock_t t1;
    Status status;
  };

private:
  /* ====== Member Field ====== */
  int current_id;
  map<string, TimeData> label_data_list;

public:
  /* ====== Constructor ======= */
  Timer() : current_id(0) {}

public:
  /* ===== Pure Functional ======= */
  Status GetStatus(const string& label);
  double GetTime(const string& label);

  /* ====== Side Effect ========== */
  void Reset();
  void Start(const string& label) ;
  void End(const string& label);

  /* ====== Input / Output ======= */
  void Display();

};
#endif



