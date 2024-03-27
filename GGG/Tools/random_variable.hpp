// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// r a n d o m _ v a r i b l e . h p p : Enthaelt ein array, das
// die Summe der Potenzen der Werte der einzelnen Ereignisse
// speichert, die zahle der ereignisse und das gesamte Gewicht
// der ereignisse. Daraus werden die Momente der Zufallsverteilung
// ermittelt und die Varianz und skewness koennen berechnet werden.                        |
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef _RANDOM_VARIABLE_HPP
#define _RANDOM_VARIABLE_HPP

#include<cstdlib>
#include<cmath>
#include<algorithm>

class random_var
{
	public:
		random_var(void);
		random_var(int size);
		~random_var(void) {delete [] moments;};
		double moment(const int);
		double mean(void);
		double variance(void);
		double skewness(void);
		void rnd_add(const double,const double);
		void operator=(random_var &);
		bool add_everything(random_var&);
		bool operator+=(random_var&);
	public:
		unsigned int max_moments;
		unsigned int events;
		double weight;
		double* moments;
};

#endif // _RANDOM_VARIABLE_HPP
