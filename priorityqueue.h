#include<iostream>
using namespace std;

#include "item.h"
#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

class PriorityQueue{
        public:
        PriorityQueue(){
			aheap = new Item[2];
			keys = new int[2];
			size = 0;
			capacity = 2;
			totalKeys = 0;
		};

		PriorityQueue(Item *array, int length); //takes a pointer to an input array with Item objects and constructs a heap from it
		PriorityQueue(const PriorityQueue &another); //a copy constructor: copies data members from another 
		~PriorityQueue(){
			if(aheap){
				delete [] aheap;
				aheap = NULL;
			}
			if(keys){
				delete [] keys;
				keys = NULL;
			}
		}; //destructor of the class: deallocates space taken by the heap

		int getCapacity(){ return capacity; };
		int getSize(){ return size; };
		void reheapifyUp(int i);
		void reheapifyDown(int i);
		void pop();
		Item getMin();
		bool push(int akey, int apr);
		bool updatePriority(int akey, int apr);  
		void print();
		bool isKey(int akey){
			if(keys[akey] == -1)
				return false;
			else
				return true;
		}

		private:
		Item* aheap = NULL;
		int *keys = NULL;
		int size;
		int capacity;
		int totalKeys;
};

#endif
