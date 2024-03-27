#ifndef _BINARY_HEAP_HPP
#define _BINARY_HEAP_HPP

#include<iostream>

#include"tree_node.hpp"

using namespace std;

class binary_heap
{
	public:
		binary_heap(void) {};
		int init(const double * values,const unsigned int size);
		void heap_sort(void);
		void heap_updt(const double old_value,const double new_value,const unsigned int position);
		bool heap_ness(void);
	private:
		void sift_down(const int l,const int r);
		unsigned int heap_size;
		idxNode<double> * heap;
	public:
		inline void setind(const int i,const int idx);
		inline void setval(const int i,const double val);
		inline int getind(const int i);
		inline double getval(const int i);
		friend ostream & operator<<(ostream & os,const binary_heap & b);
};
inline void binary_heap::setind(const int i,const int idx)
{
	heap[i].index = idx;
}
inline void binary_heap::setval(const int i,const double v)
{
	heap[i].value = v;
}
inline int binary_heap::getind(const int i)
{
	return heap[i].index;
}
inline double binary_heap::getval(const int i)
{
	return heap[i].value;
}
#endif // _BINARY_HEAP_HPP
