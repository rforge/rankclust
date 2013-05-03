#include "runEM.h"
#include <iostream>
#include "EM.h"
#include "functions.h"

#include <tuple>
#include <vector>
#include <unordered_set>
#include <unordered_map>

using namespace Rcpp ;
using namespace std ;
using namespace Eigen ;

RcppExport SEXP EMmelR(SEXP X,SEXP freq,SEXP m,SEXP g,SEXP eps,SEXP maxIt,SEXP detail)
{
    //convert données pour algo
    NumericVector freqR(freq);
    vector<int> frequence=as<vector<int> > (freqR);
    NumericMatrix XR(X);
    int const mC=as<int>(m),n(XR.nrow());
    int const maxItC=as<int>(maxIt);
    int const gC=as<int>(g);
	bool detailC=as<bool>(detail);
    double const epsC=as<double>(eps);
    vector<vector<int> > donnees(n,vector<int> (mC));

    for(int i(0);i<n;i++)
        for(int j(0);j<mC;j++)
            donnees[i][j]=XR[i+j*n];

    tuple<ArrayXd,ArrayXd,vector<vector<int> >,ArrayXXd,double,double,double,ArrayXd,ArrayXd> resEMmel;
    resEMmel=EMmel(donnees,frequence,gC,maxItC,epsC,detailC);

	double icl(get<5>(resEMmel));
	ArrayXd entropie(n);
	ArrayXXd t;
	t=get<3>(resEMmel);
	for(int i(0);i<n;i++)
	{
		entropie(i)=0;
		for(int j(0);j<gC;j++)
		{
			if(t(i,j)!=0)
				entropie(i)-=(double) frequence[i]*2*t(i,j)*log(t(i,j));
		}
		icl+=entropie(i);
	}

    return List::create(Named("referenceRank")=wrap(get<2>(resEMmel)),Named("pi")=wrap(get<0>(resEMmel)),Named("proportion")=wrap(get<1>(resEMmel)),
                        Named("tik")=wrap(get<3>(resEMmel)),Named("loglikelihood")=wrap(get<4>(resEMmel)),
                        Named("bic")=wrap(get<5>(resEMmel)),Named("icl")=wrap(icl),Named("entropie")=wrap(entropie),Named("partition")=wrap(get<8>(resEMmel)),Named("maxProb")=wrap(get<7>(resEMmel)));


}

RcppExport SEXP EMmelmultR(SEXP X,SEXP freq,SEXP m,SEXP g,SEXP eps,SEXP maxIt,SEXP detail)
{
    //convert données pour algo
    NumericVector freqR(freq),mR(m);
    vector<int> frequence=as<vector<int> > (freqR);
    vector<int> M=as<vector<int> > (mR);
    NumericMatrix XR(X);
    int const n(XR.nrow()),d(M.size());
    int const gC=as<int>(g);
    int const maxItC=as<int>(maxIt);
    double const epsC=as<double>(eps);
	bool detailC=as<bool>(detail);

    vector<vector<vector<int> > > donnees(d,vector<vector<int> > (n));
    vector<int> indM(d+1,0);
    for(int i(0);i<d;i++)
        indM[i+1]=indM[i]+M[i];

    for(int i(0);i<d;i++)
        for(int j(0);j<n;j++)
            donnees[i][j].resize(M[i]);


    int indice(0);

    for(int j(0);j<n;j++)
        for(int i(0);i<d;i++)
        {
            indice=0;
            for(int k(indM[i]);k<indM[i+1];k++)
            {
                donnees[i][j][indice]=XR[j+k*n];
                indice++;
            }
        }



    tuple<ArrayXXd,ArrayXd,vector<vector<vector<int> > >,ArrayXXd,double,double,double,ArrayXd,ArrayXd> resEMmel;
    resEMmel=EMmelMulti(donnees,frequence,gC,maxItC,epsC,detailC);

	double icl(get<5>(resEMmel));
	ArrayXd entropie(n);
	ArrayXXd t;
	t=get<3>(resEMmel);
	for(int i(0);i<n;i++)
	{
		entropie(i)=0;
		for(int j(0);j<gC;j++)
		{
			if(t(i,j)!=0)
				entropie(i)-=(double) frequence[i]*2*t(i,j)*log(t(i,j));
		}
		icl+=entropie(i);
	}

    return List::create(Named("referenceRank")=wrap(get<2>(resEMmel)),Named("pi")=wrap(get<0>(resEMmel)),Named("proportion")=wrap(get<1>(resEMmel)),
                        Named("tik")=wrap(get<3>(resEMmel)),Named("loglikelihood")=wrap(get<4>(resEMmel)),
                        Named("bic")=wrap(get<5>(resEMmel)),Named("icl")=wrap(icl),Named("entropie")=wrap(entropie),Named("partition")=wrap(get<8>(resEMmel)),Named("maxProb")=wrap(get<7>(resEMmel)));

}
