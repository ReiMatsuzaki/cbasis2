
#include <iostream>
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/symgtos_io.hpp"
using namespace std;

int main () {

    ifstream f;
    f.open(filename.c_str(), ios::in);
    string json;
    string line;
    while(getline(f, line)) {
      json += line + "\n";
    }
    SymGTOsReadJson(gtos, json);

}
