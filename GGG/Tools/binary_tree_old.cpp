// ~~~~~~~~~~~~~~~
// binary_heap.cpp
// ~~~~~~~~~~~~~~~

#include"binary_heap.hpp"

using namespace std;

int binary_heap::init(const double *values, const unsigned int size)
{  // Check the vector, if defined the destroy
	heap_size = size;
	heap = new idxNode<double>[heap_size];
	for(unsigned int i = 0; i < heap_size; i++)
	{
		heap[i].init(i,values[i]);
	}
	heap_sort();
	if(!heapness())
	{
		cerr << "Binary tree is not sorted in function binary_heap::init()!" << endl;
		exit(EXIT_FAILURE);
	}
	return(EXIT_SUCCESS);
}
// insert free(void)
void binary_heap::heap_sort(void)
{
	for(int i = heap_size-1; i >= 0; i--) sift_down(i,heap_size-1);
	for(int i = heap_size-1; i >  0; i--)
	{
		idxNode<double> a = heap[0];
		heap[0] = heap[i];
		heap[i] = a;
		sift_down(0,i-1);
	}
	//if(!heapness())
	//{
	//	cerr << "Binary tree is not sorted in function binary_heap::init()!" << endl;
	//	exit(EXIT_FAILURE);
	//}
}

void binary_heap::sift_down(const int l, const int r)
{
	idxNode<double> a = heap[l];
	int jold = l;
	int jnew = l+1;
	while(jnew <= r)
	{
		//if((jnew < r)  &&  (heap[jnew] < heap[jnew+1])) jnew++;
		//if(a >= heap[jnew]) break;
		if((jnew < r)  &&  (heap[jnew].value < heap[jnew+1].value)) jnew++;
		if(a.value >= heap[jnew].value) break;
		heap[jold] = heap[jnew];
		jold = jnew;
		jnew = 2*jnew+1;
	}
	heap[jold] = a;
}

void binary_heap::heap_updt(const double old_time,const double new_time,const unsigned int position)
{
	if(new_time < old_time)
	{
		cerr << "Attempted to update to smaller time in function binary_heap::update_top!" << endl;
		exit(EXIT_FAILURE);
	}
	heap[position].setval(new_time);
	//if(!heapness())
	//{
	//	cout << heap ;
	//}
}

bool binary_heap::heapness(void)
{
	int level;
	int tmp = heap_size+1;
	for(int i = 0;i < heap_size;i++)
	{
		level = i;
		tmp  /= 2;
		if(tmp == 0) break;
	}
	int width = 1;
	bool is_heap = true;
	for(int lvl = 0;lvl < level;lvl++)
	{
		for(int i = 0; i < width; i++)
			for(int j = 0; j < 2*width; j++)
			{
				if(2*width+j-1 >= heap_size) break;
				//if(heap[width+i-1] > heap[2*width-1+j])
				if(heap[width+i-1].value > heap[2*width-1+j].value)
				{
					is_heap = false;
					break;
				}
			}
		width *= 2;
	}
	return is_heap;
}

ostream & operator<<(ostream & os, const binary_heap & b)
{
	int k=0;
	do {
		for(int j=k; j<=2*k; j++)
		{
			if (j >= b.heap_size) break;
			os << b.heap[j].index << "," << b.heap[j].value << " ";
		}
		os << endl ;
		k = 2*k+1;
	} while(k < b.heap_size);
	return os;
}
