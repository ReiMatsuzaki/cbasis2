#include <stdexcept>
#include <iostream>
#include <time.h>
#include "timer.hpp"
#include "../utils/macros.hpp"

using namespace std;

Timer::Status Timer::GetStatus(const string& label) {

  Status status =
    (label_data_list.find(label) == 
     label_data_list.end()) ?
    kNone :
    label_data_list[label].status;

  return status;
}
double Timer::GetTime(const string& label) {
  
  if(this->GetStatus(label) != kEnd) {
    string msg;
    SUB_LOCATION(msg);
    msg += "status must be End\n";
    msg += "label : ";
    msg += label;
    throw runtime_error(msg);
  }

  TimeData data = label_data_list.find(label)->second;
  double dt = (1.0*(data.t1 - data.t0)) / CLOCKS_PER_SEC;
  return dt;
}
void Timer::Reset()  { 
  label_data_list.clear(); 
}
void Timer::Start(const string& label) {

  if(GetStatus(label) != kNone) {
    string msg = "In Timer::Start, ";
    msg += label + " is already started or ended.";
    throw std::runtime_error(msg);
  }
  TimeData data;
  data.id = current_id;
  ++current_id;
  data.t0 = clock();
  data.status = kStart;
  
  label_data_list.insert(make_pair(label, data));
}
void Timer::End(const string& label) {

  if(GetStatus(label) == kNone) {
    string msg = "In Timer::End, ";
    msg += label + " is not started.";
    throw std::runtime_error(msg);
  }
  if(GetStatus(label) == kEnd) {
    string msg = "In Timer::End, ";
    msg += label  + " is already ended.";
    throw std::runtime_error(msg);
  }

  label_data_list[label].t1 = clock();
  label_data_list[label].status = kEnd;

}
void Timer::Display() {

  cout << "label : time / s" << endl;

  for(map<string, TimeData>::iterator it = 
	label_data_list.begin();
      it != label_data_list.end();
      ++it)
    {
      string label = it->first;
      double dt = this->GetTime(label);
      
      cout << label << " : " << dt << endl;	
    }
}
