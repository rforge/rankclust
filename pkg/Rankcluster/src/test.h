#ifndef TEST_H_
#define TEST_H_

#include <vector>
#include <string>
#include <set>
#include <iostream>

#include "functions.h"

struct Rank
{
	std::vector<int> rank;
	bool isPartial;
	std::vector<int> missingIndex;
	std::set<int> missingNumber;
};


void simulMixtureISR(std::vector<std::vector<int> > &simul,int const& n,int const& m,std::vector<std::vector<int> > const& mu,std::vector<double> const& p,std::vector<double> const& prop);

double khi2(std::vector<std::vector<int> > const& data,std::vector<double> const& p,std::vector<double> const& prop,std::vector<std::vector<int> > const& mu,int const& nBoot);

double khi2partial(std::vector<Rank > &data,std::vector<double> const& p,std::vector<double> const& prop,std::vector<std::vector<int> > const& mu,int const& nBoot);

void updateD(double &divKL,std::vector<int> &index, std::vector<std::vector<std::vector<double> > > const& p1,std::vector<std::vector<std::vector<double> > >  const& p2,int const& d,int const& g,
		std::vector<double> const& proportion1,std::vector<double> const& proportion2);

void updateIndex(std::vector<int> &index,int i,std::vector<int> const& m,bool &stop);

void computePQ(std::vector<std::vector<std::vector<double> > > &p, std::vector<std::vector<std::vector<double> > > &q,std::vector<std::vector<std::vector<int> > > const& mu1,
		std::vector<std::vector<std::vector<int> > > const& mu2,std::vector<std::vector<double> > const& p1,std::vector<std::vector<double> > const& p2,std::vector<int> const& m,int d, int g);

double divKL(std::vector<int> const& m,std::vector<std::vector<std::vector<int> > > const& mu1,std::vector<std::vector<std::vector<int> > > const& mu2,
		std::vector<std::vector<double> > const& p1, std::vector<std::vector<double> >  const& p2,std::vector<double> const& proportion1,std::vector<double> const& proportion2);




#endif /* TEST_H_ */

