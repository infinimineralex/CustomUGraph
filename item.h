#include<iostream>
#include<vector>
#include<string>
using namespace std;

#ifndef ITEM_H
#define ITEM_H


class Item{
        public:
                Item(): key(0), priority(0) {};
				Item(int akey, int aval){
					key = akey;
					priority = aval;
				};
               int key;
			   int priority;
};

#endif
