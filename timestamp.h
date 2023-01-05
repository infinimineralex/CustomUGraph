#include<string>
#include<iostream>
using namespace std;

#ifndef TIMESTAMP_H
#define TIMESTAMP_H

class TimeStamp{
public:
	TimeStamp(): d(0), f(0){};
	int d; //discovery time
	int f; //finish time

};

#endif
