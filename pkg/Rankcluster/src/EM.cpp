/**
 * \file      EMmel.cpp
 * \author    Grimonprez Quentin
 * \version   1.0
 * \date      23 aout 2012
 * \brief     contient les fonctions spécifiques à l'algorithme EM pour modèle de mélange
 */


#include "functions.h"
#include "EM.h"

#include <ctime>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;


void Emelange(ArrayXXd &Tau,ArrayXXd &t,vector<ArrayXXd> const& exposantPi,vector<ArrayXXd> const& exposantUnMoinsPi,int const& m,
              int const& n,int const& g,ArrayXd const& p,ArrayXd const& prop,vector<int> const& indMu)
{
    ArrayXXd temp(g,factorial(m)/2);
    double deno(0);
    for(int i(0);i<n;i++)
    {
        for(int k(0);k<g;k++)
            temp.row(k)=prop(k)*(exposantPi[indMu[k]].row(i)*log(p(k))+exposantUnMoinsPi[indMu[k]].row(i)*log(1-p(k))).exp();

        deno=temp.sum();
        t.row(i)=temp.rowwise().sum()/deno;
        Tau.row(i)=temp.colwise().sum()/deno;
    }
}

double Mmelange(vector<ArrayXXd> const& exposantPi,vector<ArrayXXd> const& exposantUnMoinsPi,int const& g,ArrayXd &p,ArrayXd &prop,
                vector<int> &indMu,ArrayXXd const& Tau, ArrayXXd const& t,ArrayXd const& freq,VectorXd const& freqb,int const& N,
                double const& rpi)
{
    const int n(freq.size());
    double s1(0),s2(0);
    int ind(0);
    MatrixXd matA(1,n);
    vector<double> P(n),L(n);

    ArrayXd calculInter=ArrayXd::Zero(n);

    //calcul proportion
    prop=freqb.transpose()*t.matrix()/N;

    for(int k(0);k<g;k++)
    {
        //calcul pour log
        calculInter.setZero(n);
        for(int i(0);i<g;i++)
        {
            if(i!=k)
                calculInter+=((exposantPi[indMu[i]]*log(p(i))+exposantUnMoinsPi[indMu[i]]*log(1-p(i))).exp()).rowwise().sum()*prop(i);
        }

        //précalcul pour p
        matA=(freq*t.col(k)).matrix().transpose();

        for(int i(0);i<n;i++)//parcours des mu
        {
            //calcul p
            s1=(matA*(Tau*exposantPi[i]).matrix()).sum();
            s2=(matA*(Tau*(exposantPi[i]+exposantUnMoinsPi[i])).matrix()).sum();

            p(k)=s1/s2;
            P[i]=p(k);

            indMu[k]=i;
            //calcul log
            L[i]=loglikMel(calculInter,exposantPi[i],exposantUnMoinsPi[i],p(k),prop(k),rpi,freqb);

        }

        //choix du meilleur mu
        ind=max_element(L.begin(),L.end())-L.begin();
        indMu[k]=ind;
        p(k)=P[ind];

    }

    return L[ind];
}

double loglikMel(ArrayXd calculInter,ArrayXXd const& exposantPi,ArrayXXd const& exposantUnMoinsPi,double const& p, double const& prop,
                 double const& rpi,VectorXd const& freqDonnees)
{
    double ll(0);
    calculInter+=((exposantPi*log(p)+exposantUnMoinsPi*log(1-p)).exp()).rowwise().sum()*prop;
    ll=freqDonnees.dot((calculInter*rpi).log().matrix());

    return ll;
}



tuple<ArrayXd,ArrayXd,vector<vector<int> >,ArrayXXd,double,double,double,ArrayXd,ArrayXd,ArrayXd,ArrayXd,vector<vector<int> > > EMmel
(vector<vector<int> > const& donnees,vector<int> const& frequence,int const& g,int const& maxIt,double const& eps,bool const& detail)
{
	if(detail)
    	cout<<endl<<"---------------------VARIATIONAL EM ALGORITHM---------------------"<<endl;
    double t1(0),t2(0),t3(0),t4(0),tE(0),tM(0),tEM(0),tExp(0),tTotal1(0),tTotal2(0);

    tTotal1=clock();

    const int n(donnees.size()),m(donnees[0].size());
    int compteur(0);
    vector<int> tabFact(tab_factorial(m));

    //liste des ordre de présentation
    vector<int> Yindex;
    vector<vector<int> > Y(tabFact[m-1]*0.5,vector<int>(m,0));

    Yindex=listeSigma(m,tabFact);
    for (int i(0);i<tabFact[m-1]*0.5;i++)
        Y[i]=index2rank(Yindex[i],m,tabFact);

    const int nbSigma(Y.size());
    const double rpi((double) 1/nbSigma);

    //crétion vecteur frequence
    int N(frequence[0]);//N:fréquence totale
    VectorXd freqDonnees(n);
    ArrayXd freq(n);

    freq(0)=N;
    freqDonnees(0)=N;
    for (int i(1);i<n;i++)
    {
        N+=frequence[i];
        freqDonnees(i)=frequence[i];
        freq(i)=frequence[i];
    }


    t1=clock();
    //calcul des exposants
	if(detail)
    	cout<<"Pre-computing of exponents"<<endl;
    vector<ArrayXXd>  exposantPi(n,ArrayXXd(n,nbSigma)),exposantUnMoinsPi(n,ArrayXXd(n,nbSigma));
    vector<int> coeff(2,0);


    for (int k(0);k<n;k++)
    {
        for (int i(0);i<n;i++)
        {
            for (int j(0);j<nbSigma;j++)
            {
                coeff=comparaison(donnees[i],Y[j],donnees[k]);
                exposantPi[k](i,j)=coeff[1];
                exposantUnMoinsPi[k](i,j)=coeff[0]-coeff[1];
            }
        }
    }

    t2=clock();
    tExp=(double) (t2-t1)/CLOCKS_PER_SEC;
	if(detail)
    	cout<<"Computing time of exponent: "<<tExp<<" s"<<endl;

    //algo EM
	if(detail)
    	cout<<"EM ALGORITHM"<<endl<<endl<<"Iteration:"<<endl;
    double lim(-numeric_limits<double>::max()),i,loglik(0);
    ArrayXd prop(g),p(g);
    ArrayXXd Tau(n,nbSigma),t(n,g);//t: proba d'appartenance aux classes, Tau: proba que l'ordre de présentation yi soit y
    vector<double> L(maxIt+1,lim);

    //Initialisation
    prop.setConstant((double)1/g);

    vector<int> indMu(g,0);
    {
        for(int j(0);j<g;j++)
        {
            p(j)=(double) rand()/RAND_MAX*0.5+0.5;
            indMu[j]=floor((double) rand()/RAND_MAX*n);
        }
    }

	ArrayXd initProp(prop), initP(p);
	vector<vector<int> > initMu(g,vector<int> (m));
	for(int i(0); i<g;i++)
		initMu[i]=donnees[indMu[i]];


	//start EM
    t1=clock();

    for(i=0;i<maxIt;i++)
    {
		if(detail)
        	cout<<"*";
        t3=clock();
        Emelange(Tau,t,exposantPi,exposantUnMoinsPi,m,n,g,p,prop,indMu);
        t4=clock();
        tE+=t4-t3;

        t3=clock();
        L[i+1]=Mmelange(exposantPi,exposantUnMoinsPi,g,p,prop,indMu,Tau,t,freq,freqDonnees,N,rpi);
        t4=clock();
        tM+=t4-t3;


        loglik=L[i+1];
        if(fabs(L[i]-L[i+1])<eps)
            break;
        compteur++;
    }

    t2=clock();
    tEM=(t2-t1)/CLOCKS_PER_SEC;
	if(detail)
   		cout<<endl<<"Computing time EM: "<<tEM <<" s"<<endl;

    vector<vector<int> > resMu(g,vector<int>(m));
    for(int k(0);k<g;k++)
    {
        resMu[k]=donnees[indMu[k]];
        if(p(k)<0.5)
        {
            p(k)=1-p(k);
            inverseRang(resMu[k],m);
        }
    }
    tTotal2=clock();
    tTotal2-=tTotal1;
    tTotal2/=CLOCKS_PER_SEC;
	tE/=CLOCKS_PER_SEC;
	tM/=CLOCKS_PER_SEC;
    double bic(-2*loglik+(3*g-1)*log(N)),icl(bic-2*(t.log()*t).sum());

	if(detail)
	{
		cout<<"Total computing time:"<<tTotal2<<endl;
		cout<<"Total computing time for step E: "<<tE<<" s. For one step E: "<<(double) tE/compteur<<"s"<<endl;
		cout<<"Total computing time for step M: "<<tM<<" s. For one step M: "<<(double) tM/compteur<<"s"<<endl;
	}
/*
	cout<<"Results:"<<endl;
    cout<<"pi:"<<p<<endl;
    cout<<"proportion:"<<prop<<endl;
    cout<<"mu: "<<endl<<resMu<<endl;
    cout<<"LL: "<<loglik<<endl;
    cout<<"BIC: "<<bic<<endl;
    cout<<"ICL: "<<icl<<endl;*/

	ArrayXd z=ArrayXd::Zero(n),prob=ArrayXd::Ones(n);
	double normconst((double) 2/factorial(m));
    if(g>1)
	{//calcul partition
		double max;
		for(int i(0);i<n;i++)
		{
			max=t(i,0);
			for(int k(1);k<g;k++)
			{
				if(t(i,k)>max)
				{
					max=t(i,k);
					z(i)=k;
				}
			}
		}
	}

    for(int i(0);i<n;i++)
		prob(i)=(exposantPi[indMu[z(i)]].row(i)*log(p(z(i)))+exposantUnMoinsPi[indMu[z(i)]].row(i)*log(1-p(z(i)))).exp().sum()*normconst;


    tuple<ArrayXd,ArrayXd,vector<vector<int> >,ArrayXXd,double,double,double,ArrayXd,ArrayXd,ArrayXd,ArrayXd,vector<vector<int> > > resultat(p,prop,resMu,t,loglik,bic,icl,prob,z,initProp,initP,initMu);

    return(resultat);
}


void EmelangeMult(vector<ArrayXXd> &Tau,ArrayXXd &t,vector<vector<ArrayXXd> > const& exposantPi,
                  vector<vector<ArrayXXd> > const& exposantUnMoinsPi,vector<int> const& m,int const& n,int const& g,
                  ArrayXXd const& p,ArrayXd const& prop,vector<int> const& indMu)
{
    int const d(m.size());

    for(int i(0);i<n;i++)//pour chaque xi
    {
        ArrayXd temp2=ArrayXd::Constant(g,1);
        for(int j(0);j<d;j++)//pour chaque dimension
        {
            ArrayXXd temp(g,factorial(m[j])/2);
            for(int k(0);k<g;k++)//parcours composante
            {
                temp.row(k)=(exposantPi[j][indMu[k]].row(i)*log(p(k,j))+exposantUnMoinsPi[j][indMu[k]].row(i)*log(1-p(k,j))).exp();
                temp2(k)*=temp.row(k).sum();
                temp.row(k)*=prop(k);
            }
            Tau[j].row(i)=temp.colwise().sum()/temp.sum();
        }
        t.row(i)=(temp2*prop)/(temp2*prop).sum();
    }

}


double MmelangeMult(vector<vector<ArrayXXd> > const& exposantPi,vector<vector<ArrayXXd> > const& exposantUnMoinsPi,int const& g,
                    ArrayXXd &p,ArrayXd &prop,vector<int> &indMu,vector<ArrayXXd> const& Tau, ArrayXXd const& t,ArrayXd const& freq,
                    VectorXd const& freqb,int const& N,double const& rpi)
{
    double s1(0),s2(0);
    int const n(freq.size()),d(p.cols());
    int ind(0);
    MatrixXd matA(1,n);
    vector<double> L(n);
    ArrayXXd P(n,d);
    ArrayXd calculInter=ArrayXd::Zero(n);
    ArrayXd calculInter2(n);

    //maj proportion
    prop=freqb.transpose()*t.matrix()/N;

    for(int k(0);k<g;k++)//parcours des composantes
    {
        calculInter.setZero(n);
        for(int i(0);i<g;i++)
        {
            calculInter2.setOnes(n);
            if(i!=k)
            {
                for(int j(0);j<d;j++)
                    calculInter2*=((exposantPi[j][indMu[i]]*log(p(i,j))+exposantUnMoinsPi[j][indMu[i]]*log(1-p(i,j))).exp()).rowwise().sum();
                calculInter+=calculInter2*prop(i);
            }
        }

        matA=(freq*t.col(k)).matrix().transpose();

        for(int i(0);i<n;i++)//parcours des mu
        {
            for(int j(0);j<d;j++)//parcours des dim
            {
                s1=(matA*(Tau[j]*exposantPi[j][i]).matrix()).sum();
                s2=(matA*(Tau[j]*(exposantPi[j][i]+exposantUnMoinsPi[j][i])).matrix()).sum();

                p(k,j)=s1/s2;
                P(i,j)=p(k,j);
            }

            indMu[k]=i;

            L[i]=loglikMelMult(calculInter,exposantPi,exposantUnMoinsPi,indMu,p,prop,freqb,rpi,k);
        }

        ind=max_element(L.begin(),L.end())-L.begin();
        indMu[k]=ind;
        p.row(k)=P.row(ind);
    }

    return L[ind];
}


double loglikMelMult(ArrayXd calculInter,vector<vector<ArrayXXd> > const& exposantPi,vector<vector<ArrayXXd> > const& exposantUnMoinsPi,
                     vector<int> const& indiceMu,ArrayXXd const& p,ArrayXd const& prop,VectorXd const& freqDonnees,double const& rpi,
                     int const& numComp)
{
    double ll(0);
    int const g(prop.size()),d(p.cols());

    ArrayXXd li=ArrayXXd::Constant(freqDonnees.size(),g,1);

    ArrayXd calculInter2=ArrayXd::Ones(freqDonnees.size());
    for(int j(0);j<d;j++)
        calculInter2*=((exposantPi[j][indiceMu[numComp]]*log(p(numComp,j))+exposantUnMoinsPi[j][indiceMu[numComp]]*log(1-p(numComp,j))).exp()).rowwise().sum();
    calculInter+=calculInter2*prop(numComp);

    ll=freqDonnees.dot((calculInter/rpi).log().matrix());
    return ll;
}


tuple<ArrayXXd,ArrayXd,vector<vector<vector<int> > >,ArrayXXd,double,double,double,ArrayXd,ArrayXd,ArrayXd,ArrayXXd,vector<vector<vector<int> > > >
EMmelMulti(vector<vector<vector<int> > > const& donnees,vector<int> const& frequence,int const& g,int const& maxIt,double const& eps,bool const& detail)
{
    cout<<endl<<"---------------------VARIATIONAL EM ALGORITHM---------------------"<<endl;
    const int d(donnees.size()),n(donnees[0].size());
    int compteur(0);
    double t1(0),t2(0),t3(0),t4(0),tE(0),tM(0),tEM(0),tExp(0);//pour mesure du temps de calcul
    vector<int> m(d);
    for(int i(0);i<d;i++)
        m[i]=donnees[i][0].size();
    int M(*max_element(m.begin(),m.end()));
    vector<int> tabFact(tab_factorial(M));

    //----------liste des ordre de présentation
    vector<vector<vector<int> > > Y(d,vector<vector<int> > (n));

    for(int j(0);j<d;j++)
    {
        vector<int> listeIndex(tabFact[m[j]-1]*0.5);
        vector<vector<int> > vecTemp(tabFact[m[j]-1]*0.5);
        listeIndex=listeSigma(m[j],tabFact);
        for (int i(0);i<tabFact[m[j]-1]*0.5;i++)
            vecTemp[i]=index2rank(listeIndex[i],m[j],tabFact);
        Y[j]=vecTemp;
    }


    //crétion vecteur frequence
    int N(frequence[0]);//N:fréquence totale
    VectorXd freqDonnees(n);
    ArrayXd freq(n);

    freq(0)=N;
    freqDonnees(0)=N;
    for (int i(1);i<n;i++)
    {
        N+=frequence[i];
        freqDonnees(i)=frequence[i];
        freq(i)=freqDonnees(i);
    }


    //------------------------------calcul des exposants
	if(detail)
    	cout<<"Pre-computing of exponents"<<endl;
    t1=clock();
    vector<vector<ArrayXXd> >  exposantPi(d,vector<ArrayXXd>(n)),exposantUnMoinsPi(d,vector<ArrayXXd>(n));
    vector<int> coeff(2,0);

    for(int l(0);l<d;l++)
    {
        ArrayXXd expPiTemp(n,Y[l].size()),expUnMoinsPiTemp(n,Y[l].size());
        for(int k(0);k<n;k++)//parcours des mu
        {
            for(int i(0);i<n;i++)//parcours des x
            {
                for(unsigned int j(0);j<Y[l].size();j++)//parcours des y
                {
                    coeff=comparaison(donnees[l][i],Y[l][j],donnees[l][k]);
                    expPiTemp(i,j)=coeff[1];
                    expUnMoinsPiTemp(i,j)=coeff[0]-coeff[1];
                }
            }
            exposantPi[l][k]=expPiTemp;
            exposantUnMoinsPi[l][k]=expUnMoinsPiTemp;
        }
    }
    t2=clock();
    tExp=(t2-t1)/CLOCKS_PER_SEC;
	if(detail)
    	cout<<"Computing time of exponent: "<<tExp<<" s"<<endl;


    //-------------------------------algo EM
	if(detail)
    	cout<<"EM Algorithm"<<endl;
    int lim(numeric_limits<double>::min()),i;
    pair<vector<ArrayXXd>,ArrayXXd> resE;
    tuple<ArrayXd,ArrayXXd,double,vector<int> > resM;
    ArrayXd prop=ArrayXd::Constant(g,(double) 1/g);
    ArrayXXd p(g,d);
    vector<ArrayXXd> Tau(d);
    ArrayXXd t(n,g);
    for(int i(0);i<d;i++)
        Tau[i].resize(n,tabFact[m[i]-1]/2);

    double loglik(0),rpi(1);
    for(int i(0);i<d;i++)
        rpi*=tabFact[m[i]-1]/2;


    vector<double> L(maxIt+1,lim);
    t1=clock();
    //Initialisation
    vector<int> indMu(g,0);
    for(int i(0);i<g;i++)
    {
        for(int j(0);j<d;j++)
            p(i,j)=(double) rand()/RAND_MAX*0.5+0.5;
        indMu[i]=floor((double) rand()/RAND_MAX*n);
    }

	ArrayXd initProp(prop);
	ArrayXXd initP(p);
	vector<vector<vector<int> > > initMu(g,vector<vector<int> > (d) ); 

	for(int j(0);j<d;j++)
    {
        for(int k(0);k<g;k++)
        {
            initMu[k][j]=donnees[j][indMu[k]];
        }
    }


	//start EM
    for(i=0;i<maxIt;i++)
    {
        cout<<"*";
        t3=clock();
        EmelangeMult(Tau,t,exposantPi,exposantUnMoinsPi,m,n,g,p,prop,indMu);
        t4=clock();
        tE+=t4-t3;


        t3=clock();
        L[i+1]=MmelangeMult(exposantPi,exposantUnMoinsPi,g,p,prop,indMu,Tau,t,freq,freqDonnees,N,rpi);
        t4=clock();
        tM+=t4-t3;

        loglik=L[i+1];
        if(fabs(L[i]-L[i+1])<eps)
            break;
        compteur++;

    }
    cout<<endl;
    t2=clock();
    tEM=(t2-t1)/CLOCKS_PER_SEC;
	if(detail)
    	cout<<"Computing time of EM: "<<tEM <<" s"<<endl;

    vector<vector<vector<int> > > resMu(g,vector<vector<int> >(d));

    for(int j(0);j<d;j++)
    {
        vector<int> rang(m[j]+1);
        for(int k(0);k<g;k++)
        {
            resMu[k][j]=donnees[j][indMu[k]];
            if(p(k,j)<0.5)
            {
                p(k,j)=1-p(k,j);
                inverseRang(resMu[k][j],m[j]);
            }
        }
    }
    double bic(BIC(loglik,N,2*g*d+g-1)),icl(bic-2*(t.log()*t).sum());

	tE/=CLOCKS_PER_SEC;
    tM/=CLOCKS_PER_SEC;
	if(detail)
	{
    	cout<<"Total computing time for step E: "<<tE<<" s. For one E step: "<<(double) tE/compteur<<"s"<<endl;
    	cout<<"Total computing time for step M: "<<tM<<" s. For one M step: "<<(double) tM/compteur<<"s"<<endl;
	}
/*
	cout<<"Results:"<<endl;

    cout<<"pi:"<<endl<<p<<endl;
    cout<<"proportion:"<<endl<<prop<<endl;
    cout<<"mu: "<<endl<<resMu<<endl;
    cout<<"LL:"<<loglik<<endl;
    cout<<"BIC: "<<bic<<endl;
    cout<<"ICL: "<<icl<<endl;
	cout<<"tik:"<<endl<<t<<endl;*/


	ArrayXd z=ArrayXd::Zero(n);
    if(g>1)
	{//calcul partition
		double max;

		for(int i(0);i<n;i++)
		{
			max=t(i,0);
			z(i)=0;
			for(int k(1);k<g;k++)
			{
				if(t(i,k)>max)
				{
					max=t(i,k);
					z(i)=k;
				}
			}
		}
	}

	ArrayXd prob(n);
	double normconst(1);
	for(int i(0);i<d;i++)
		normconst*=factorial(m[i])/2;

	for(int i(0);i<n;i++)
	{
		prob(i)=1;
		for(int j(0);j<d;j++)
			prob(i)*=(exposantPi[j][indMu[z(i)]].row(i)*log(p(z(i),j))+exposantUnMoinsPi[j][indMu[z(i)]].row(i)*log(1-p(z(i),j))).exp().sum();

		prob(i)/=normconst;
	}



    tuple<ArrayXXd,ArrayXd,vector<vector<vector<int> > >,ArrayXXd,double,double,double,ArrayXd,ArrayXd,ArrayXd,ArrayXXd,vector<vector<vector<int> > > >
 resultat(p,prop,resMu,t,loglik,bic,icl,prob,z,initProp,initP,initMu);

    return(resultat);
}

