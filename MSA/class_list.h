#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sys/types.h>
//#include <dirent.h>
#include <vector>
#include <errno.h>
#include <list>
#include <algorithm>
#include <functional>
#include <string>
#include <sstream>
//#include <cstddef>
#include <cmath>
#include <numeric>

using namespace std;
/*
class Report {
  string *path;
  ofstream output;
public:
  Report (string);
  void word (string);
  void close();
};

Report::Report (string a) {
  output.open(a.c_str());
  if(output.is_open())
    path=&a;
  else
    cout << "fails to open " << a << " !\n";
}

void Report::word (string a) {
  output << a << endl;
}

void Report::close () {
  output.close();
}


class Numtostr{
public:
  string convert (unsigned int);
  string convert_float (double);
  //	string convert_fix (unsigned int, int);
};

string Numtostr::convert (unsigned int a) {
  string num;
  num.clear();
  unsigned int k=0, record[8];
  do{
    record[k++]=a%10;
    a/=10;
  }while(a);
  while(k--) num+=record[k]+48;
  return num;
}

string Numtostr::convert_fix (unsigned int a, int l) {
unsigned int b;
string num;
num.clear();	
while(l--) {
b=a/pow(10.0,l);
num+=b%10+48;
}
return num;
}

string Numtostr::convert_float (double a) {
  unsigned int b, k, record[8];
  string num;
  num.clear();
  b = a*1000;
  k = b%10;
  b/=10;
  if (k>4)
    b++;
  k = 0;
  do{
    record[k++]=b%10;
    b/=10;
  }while(b);
  while(k-->2) num+=record[k]+48;
  num+='.';
  num+=record[k]+48;
  num+=record[0]+48;
  return num;
}
*/
