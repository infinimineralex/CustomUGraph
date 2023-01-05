#include<iostream>
#include<string>

using namespace std;

#include "priorityqueue.h"

void PriorityQueue::print(){
	for(int i = 0; i < size; i++){
		cout << aheap[i].priority << " " ;

	}
	cout << endl;
	for(int i = 0; i < size; i++){
		cout << aheap[i].key << " " ;

	}
	cout << endl;
	for(int i = 0; i < totalKeys; i++){
		cout << keys[i] << " " ;

	}
	cout << endl;
}
PriorityQueue::PriorityQueue(Item *array, int length){ //takes a pointer to an input array with Item objects and constructs a heap from it
	if(array == NULL)
		return;
	int firstParent = (length - 2)/2;
	if(aheap){
		delete [] aheap;
		aheap = NULL;
	}
	if(keys)
		delete [] keys;
	aheap = new Item[length];
	keys = new int[length];
	size = length;
	capacity = length;
	totalKeys = length;

	for(int i = 0; i < length; i++){
		aheap[i] = array[i];
		int ind = aheap[i].key;
		keys[ind] = i;
	}
	for(int j = firstParent; j >= 0; j--)
		reheapifyDown(j);
	
}
		


void PriorityQueue::reheapifyUp(int i){
	if(i <= 0)
		return;
	int parent = (i - 1)/2;
	if(parent >= 0){
		int curKey = aheap[i].key;
		int parKey = aheap[parent].key;
		if(aheap[parent].priority > aheap[i].priority)
		{
			Item temp = aheap[parent];
			aheap[parent] = aheap[i];
			aheap[i] = temp;
			keys[curKey] = parent;
			keys[parKey] = i;
		}
	}
	reheapifyUp(parent);
}//reheapifyUp
void PriorityQueue::reheapifyDown(int i){
	if(i >= size)
		return;
	int left = 2*i + 1;
	int right = 2*i + 2;
	int minChild = left;
	if(left < size){
		if(right < size){
			if(aheap[right].priority < aheap[left].priority)
				minChild = right;

		}
		if(aheap[minChild].priority < aheap[i].priority){
			int keyCur = aheap[i].key;
			int keyChild = aheap[minChild].key;
			Item temp = aheap[minChild];
			aheap[minChild] = aheap[i];
			aheap[i] = temp;
			keys[keyCur] = minChild;
			keys[keyChild] = i;

			reheapifyDown(minChild);
		}
	}//if left < size
	
}//PriorityQueue
		
void PriorityQueue::pop(){
	if(size == 0)
		return;
	keys[aheap[0].key] = -1;
	int smallestKey = aheap[size - 1].key;
	keys[smallestKey] = 0;
	aheap[0] = aheap[size - 1];
	reheapifyDown(0);
	size--;
}
		
Item PriorityQueue::getMin(){
	return aheap[0];
}
		
bool PriorityQueue::push(int akey, int apr){
	if(akey != size)
		return false;
	if(size >= capacity){
		capacity = capacity * 2;
		Item *temp = new Item[capacity];
		int *tempKeys = new int[capacity];
		for(int i = 0; i < size; i++){
			temp[i] = aheap[i];
			tempKeys[i] = keys[i];
		}//for
		delete [] aheap;
		delete [] keys;
		aheap = temp;
		keys = tempKeys;
		temp = NULL;
		tempKeys = NULL;
	}//if heap is full

	aheap[size] = Item(akey, apr);
	keys[size] = size;
	reheapifyUp(size);
	size++;
	totalKeys++;
	return true;
	
}//push
		
bool PriorityQueue::updatePriority(int akey, int apr){
	int index = keys[akey];
	if(aheap[index].priority < apr)
		return false;
	aheap[index].priority = apr;
	reheapifyUp(index);
	return true;
}     	