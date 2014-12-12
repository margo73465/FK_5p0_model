#include <stdlib.h>
#include <iostream>
#include <vector>

using namespace std;

int main() {

  vector<int> test(10);
  int i;

  for (i = 0; i < 10; i++) {
    test[i] = i;
    cout << test[i] << endl;
  }

  double stim_time = 1000.03;
  int n = 1000;

  if (n == stim_time) {
    cout << "HELLO" << endl;
  }


  return 0;
}
