// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// r a n d o m _ v a r i b l e . c p p
// this version was corrected by T. Erdmann
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "random_variable.hpp"

using namespace std;

random_var::random_var(void)
{
	max_moments = 3;
	events = 0;
	weight = 0.0;
	moments = new double[max_moments];
	for (unsigned int i = 0; i< max_moments; i++) moments[i] = 0.0;
}

random_var::random_var(int size)
{// need three to calculate skewness
	if (size >= 3)
		max_moments = size;
	else
		max_moments = 3;
	events  = 0;
	weight  = 0.0;
	moments = new double[max_moments];
	for (unsigned int i = 0; i< max_moments; i++) moments[i] = 0.0;
}

// THAT WAS IN THE OLD VERSION:
//void random_var::rnd_add(const double vl,const double wght)
//{
//	events++;
//	weight += wght;
//	for (int i = 0; i < max_moments; i++) moments[i] += wght*pow(vl,i+1);
//}

void random_var::rnd_add(const double vl,const double wght)
{
	events++;
	weight += wght;
	double summand = wght;
	for (int i = 0; i < max_moments; i++)
	{
		summand *= vl;
		moments[i] += summand;
	}
}

double random_var::moment(const int i)
{
	return moments[i-1]/weight;
}

double random_var::mean(void)
{
	return moments[0]/weight;
}

double random_var::variance(void)
{
	double av    = moments[0]/weight;
	double av_sq = moments[1]/weight;
	return sqrt(av_sq - av*av);
}

double random_var::skewness(void)
{
	double av    = moments[0]/weight;
	double av_sq = moments[1]/weight;
	double av_cb = moments[2]/weight;
	double sk = av_cb - (3.0*av_sq - 2.0*av*av)*av;
	if (sk > 0.0)
	{
		sk =   pow(sk,0.33333333333);
	}
	else if (sk < 0.0)
	{
		sk = - pow(-sk,0.33333333333);
	}
	else
	{
		sk = 0.0;
	}
	return sk;
}

// ~~~~~~~~~
// Operators
// ~~~~~~~~~

void random_var::operator=(random_var & rnd)
{
	if (this == &rnd) return;
	delete [] moments ;
	max_moments = rnd.max_moments;
	events = rnd.events;
	weight = rnd.weight;
	moments = new double[max_moments];
	for (unsigned int i = 0; i < max_moments; i++) moments[i] = rnd.moments[i];
	return;
}

ostream & operator<<(ostream & os, const random_var & rv)
{
	for(unsigned int moment = 0; moment < rv.max_moments; moment++)
	{
		os << rv.moments[moment] << " ";
	}
	os << rv.events << " " << rv.weight;	
	return os;
}

istream & operator>>(istream & is, random_var & rv)
{
	for(unsigned int moment = 0; moment < rv.max_moments; moment++)
	{
		is >> rv.moments[moment];
	}
	is >> rv.events >> rv.weight;	
	return is;
}
