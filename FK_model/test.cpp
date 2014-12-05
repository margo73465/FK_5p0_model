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


  return 0;
}
