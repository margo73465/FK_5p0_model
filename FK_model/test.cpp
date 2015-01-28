#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main() {
  ifstream objective("scaled_MAP.dat", ios::in);
  // ifstream objective("scaled_MAP.txt", ios::in);
  int i;
  double obj_u;

  for (i = 0; i < 499; i++) {
    objective >> obj_u;
    cout << obj_u << endl;
  }

  // vector<int> test(10);
  // int i;

  // for (i = 0; i < 10; i++) {
  //   test[i] = i;
  //   cout << test[i] << endl;
  // }

  // double stim_time = 1000.03;
  // int n = 1000;

  // if (n == stim_time) {
  //   cout << "HELLO" << endl;
  // }


  return 0;
}
