/**
 * \file      EMmel.h
 * \author    Grimonprez Quentin
 * \version   1.0
 * \date      23 aout 2012
 * \brief     contient les prototypes de EM.cpp
 */

#ifndef EM_H_INCLUDED
#define EM_H_INCLUDED

#include <vector>
#include <Eigen/Dense>


/**
 *---------------�tape E pour mod�le de m�lange----------------
 * \fn void Emelange(Eigen::ArrayXXd &Tau,Eigen::ArrayXXd &t,std::vector<Eigen::ArrayXXd> const& exposantPi,
              std::vector<Eigen::ArrayXXd> const& exposantUnMoinsPi,int const& m,int const& n,int const& g,Eigen::ArrayXd const& p,
              Eigen::ArrayXd const& prop,std::vector<int> const& indMu);
 * \brief Etape E de l'algorithme EM variationnel pour modele de melange univari�
 * \param Tau tableau contenant les probabilit�s que y soit l'ordre de pr�sentation de xi pour chauqe individu xi et chauqe y
 * \param t tableau contenant les probabilit�s d'appartenance aux classes pour chaque xi
 * \param exposantPi tableau du nombre de bonne comparaison pour tt les rangs des donn�es et tt les y de la liste
 * \param exposantUnMoinsPi tableau de A(x,y)-G(x,y,mu) pour tt les rangs x des donn�es et tt les y de la liste
 * \param m taille du rang
 * \param n nombre derangs diff�rents dans les donn�es
 * \param g nombre de groupe
 * \param p vecteur de taille g contenant les proportions du m�langes
 * \param prop vecteur de taille g contenant le param�tre p de chacune des composantes
 * \param indMu vecteur de taille g contenant les indices des param�tres mu de chaque composante
 *
 * la fonction ne retourne rien
 * les parametres Tau et t sont modifies par la fonction
 *
 */
void Emelange(Eigen::ArrayXXd &Tau,Eigen::ArrayXXd &t,std::vector<Eigen::ArrayXXd> const& exposantPi,
              std::vector<Eigen::ArrayXXd> const& exposantUnMoinsPi,int const& m,int const& n,int const& g,Eigen::ArrayXd const& p,
              Eigen::ArrayXd const& prop,std::vector<int> const& indMu);

/**
 *---------------�tape M pour mod�le de m�lange----------------
 * \fn double Mmelange(std::vector<Eigen::ArrayXXd> const& exposantPi,std::vector<Eigen::ArrayXXd> const& exposantUnMoinsPi,int const g,
                Eigen::ArrayXd &p,Eigen::ArrayXd &prop,std::vector<int> &indMu,Eigen::ArrayXXd const& Tau, Eigen::ArrayXXd const& t,
                Eigen::ArrayXd const& freq,Eigen::VectorXd const& freqb,int const& N,double const& rpi);
 * \brief etape M de l'algorithme EM variationnel pour mod�le de mlange univari�
 * \param exposantPi tableau du nombre de bonne comparaison pour tt les rangs des donn�es et tt les y de la liste
 * \param exposantUnMoinsPi tableau de A(x,y)-G(x,y,mu) pour tt les rangs x des donn�es et tt les y de la liste
 * \param g nombre de groupe
 * \param p vecteur de taille g contenant les proportions du m�langes
 * \param prop vecteur de taille g contenant le param�tre p de chacune des composantes
 * \param indMu vecteur de taille g contenant les indices des param�tres mu de chaque composante
 * \param Tau tableau \see Emelange
 * \param t tableau \see Emelange
 * \param freq tableau des fr�quences des donn�es
 * \param freqb tableau des fr�quences des donn�es sous un autre type de stockage
 * \param N frequence totale
 * \param rpi m1!*...*md! (pour la constante de normalisation)
 * \return la logvraisemblance du mod�le
 *
 * la fonction mets � jour les parametres p,prop,indMu  et retourne la logvraisemblance asoci�e
 *
 */
double Mmelange(std::vector<Eigen::ArrayXXd> const& exposantPi,std::vector<Eigen::ArrayXXd> const& exposantUnMoinsPi,int const& g,
                Eigen::ArrayXd &p,Eigen::ArrayXd &prop,std::vector<int> &indMu,Eigen::ArrayXXd const& Tau, Eigen::ArrayXXd const& t,
                Eigen::ArrayXd const& freq,Eigen::VectorXd const& freqb,int const& N,double const& rpi);

/**-----------Calcul de la logvraisemblance pour mod�le de m�lange------------
 * \fn double loglikMel(Eigen::ArrayXd calculInter,Eigen::ArrayXXd const& exposantPi,Eigen::ArrayXXd const& exposantUnMoinsPi,
                 double const& p,double const& prop,double const& rpi,Eigen::VectorXd const& freqDonnees);
 * \brief calcul de la logvraisemblance dans l'algorithme EM pour modeles de melange univarie
 * \param calculInter calcul inermediaire effectue dans l'etape M ( \see Mmelange)
 * \param exposantPi G(x,y,mu)( \see comparaison) pour chaque rang x des donn�es, chaque y de liste Sigma( \see listeSigma), et tous les mu possibles ( \see listeMu)
 * \param exposantUnMoinsPi A(x,y)-G(x,y,mu)( \see comparaison) pour chaque rang x des donn�es, chaque y de liste Sigma( \see listeSigma), et tous les mu possibles ( \see listeMu)
 * \param m taille du rang
 * \param p vecteur probabilit� d'une bonne comparaison pour chaque classe
 * \param prop vecteur des proportions de chaque classe
 * \param rpi m!/2 (pour la constante de normalisation)
 * \param freqDonnees frequence des rangs des donn�es
 * \return logvraisemblance du mod�le
 *
 * cette fonction est propre a l'algorithme et ne peut etre utilisee en dehors
 *
 */
double loglikMel(Eigen::ArrayXd calculInter,Eigen::ArrayXXd const& exposantPi,Eigen::ArrayXXd const& exposantUnMoinsPi, double const& p,
                 double const& prop,double const& rpi,Eigen::VectorXd const& freqDonnees);


/**
 *---------------ALGO EM mod�le de m�lange----------------
 * \fn std::tuple<Eigen::ArrayXd,Eigen::ArrayXd,std::vector<std::vector<int> >,Eigen::ArrayXXd,double,double,double>
EMmel(std::vector<std::vector<int> > const& donnees,std::vector<int> const& frequence,int const& g,double const& eps=0.0001,
      int const& maxIt=40);

 * \brief Algorithme EM pour modeles de melange univarie
 * \param donnees vecteur o� chaque composante est un tableau de donnees de rang import�es( \see import)
 * \param frequence frequence des donn�es
 * \param g nombre de groupe
 * \param eps seuil pour crit�re d'arr�t sur la vraisemblance
 * \param maxIt nombre maximal d'it�ration pour un EM sur un mu
 * \return un tuple contenant prop,p,mu,les tik
 */
std::tuple<Eigen::ArrayXd,Eigen::ArrayXd,std::vector<std::vector<int> >,Eigen::ArrayXXd,double,double,double,Eigen::ArrayXd,Eigen::ArrayXd,Eigen::ArrayXd,Eigen::ArrayXd,std::vector<std::vector<int> > >
EMmel(std::vector<std::vector<int> > const& donnees,std::vector<int> const& frequence,int const& g,int const& maxIt=40,
      double const& eps=0.0001,bool const& detail=false);


/**
 *---------------etape E mod�le de m�lange----------------
 * \fn void EmelangeMult(std::vector<Eigen::ArrayXXd> &Tau,Eigen::ArrayXXd &t,std::vector<std::vector<Eigen::ArrayXXd> > const& exposantPi,
                  std::vector<std::vector<Eigen::ArrayXXd> > const& exposantUnMoinsPi,std::vector<int> const& m,int const& n,int const& g,
                  Eigen::ArrayXXd const& p,Eigen::ArrayXd const& prop,std::vector<int> const& indMu);
 * \brief Etape E de l'algorithme EM variationnel pour modele de melange multivari�
 * \param Tau tableau contenant les probabilit�s que y soit l'ordre de pr�sentation de xi pour chauqe individu xi et chauqe y pour hauqe dimension
 * \param t tableau contenant les probabilit�s d'appartenance aux classes pour chaque xi
 * \param exposantPi G(x,y,mu)( \see comparaison) pour chaque rang x des donn�es, chaque y de liste Sigma( \see listeSigma), et tous les mu possibles ( \see listeMu)
 * \param exposantUnMoinsPi A(x,y)-G(x,y,mu)( \see comparaison) pour chaque rang x des donn�es, chaque y de liste Sigma( \see listeSigma), et tous les mu possibles ( \see listeMu)
 * \param m taille des rangs pour chaque dimension
 * \param n nombre de rangs diff�rents
 * \param g nombre de groupe
 * \param p probabilite d'une bonne comparaison pour chaque composante et dimension
 * \param prop proportion du melange
 * \param indMu tableau avec les indices des mu actuels
 *
 * La fonction ne retourne rien
 * Les parametre Tau et t sont modifies
 *
 */
void EmelangeMult(std::vector<Eigen::ArrayXXd> &Tau,Eigen::ArrayXXd &t,std::vector<std::vector<Eigen::ArrayXXd> > const& exposantPi,
                  std::vector<std::vector<Eigen::ArrayXXd> > const& exposantUnMoinsPi,std::vector<int> const& m,int const& n,int const& g,
                  Eigen::ArrayXXd const& p,Eigen::ArrayXd const& prop,std::vector<int> const& indMu);

/**
 *---------------etape m mod�le de m�lange----------------
 * \fn double MmelangeMult(std::vector<std::vector<Eigen::ArrayXXd> > const& exposantPi,
                    std::vector<std::vector<Eigen::ArrayXXd> > const& exposantUnMoinsPi,int const& g,Eigen::ArrayXXd &p,
                    Eigen::ArrayXd &prop,std::vector<int> &indMu,std::vector<Eigen::ArrayXXd> const& Tau, Eigen::ArrayXXd const& t,
                    Eigen::ArrayXd const& freq,Eigen::VectorXd const& freqb,int const& N,double const& rpi);
 * \brief etape M de l'algorithme EM variationnel pour mod�le de mlange multivari�
 * \param exposantPi tableau du nombre de bonne comparaison pour tt les rangs des donn�es et tt les y de la liste
 * \param exposantUnMoinsPi tableau de A(x,y)-G(x,y,mu) pour tt les rangs x des donn�es et tt les y de la liste
 * \param g nombre de groupe
 * \param p vecteur de taille g contenant les proportions du m�langes
 * \param prop vecteur de taille g contenant le param�tre p de chacune des composantes
 * \param indMu vecteur de taille g contenant les indices des param�tres mu de chaque composante
 * \param Tau tableau \see EmelangeMult
 * \param t tableau \see EmelangeMult
 * \param freq tableau des fr�quences des donn�es
 * \param freqb tableau des fr�quences des donn�es sous un autre type de stockage
 * \param N frequence totale
 * \param rpi m1!/2*...*md!/2 (pour la constante de normalisation
 * \return la logvraisemblance du mod�le
 *
 * la fonction mets � jour les parametres p,prop,indMu  et retourne la logvraisemblance asoci�e
 *
 */
double MmelangeMult(std::vector<std::vector<Eigen::ArrayXXd> > const& exposantPi,
                    std::vector<std::vector<Eigen::ArrayXXd> > const& exposantUnMoinsPi,int const& g,Eigen::ArrayXXd &p,
                    Eigen::ArrayXd &prop,std::vector<int> &indMu,std::vector<Eigen::ArrayXXd> const& Tau, Eigen::ArrayXXd const& t,
                    Eigen::ArrayXd const& freq,Eigen::VectorXd const& freqb,int const& N,double const& rpi);

/**-----------Calcul de la logvraisemblance pour mod�le de m�lange------------
 * \fn double loglikMelMult(Eigen::ArrayXd calculInter,std::vector<std::vector<Eigen::ArrayXXd> > const& exposantPi,
                     std::vector<std::vector<Eigen::ArrayXXd> > const& exposantUnMoinsPi,
                     std::vector<int> const& indiceMu,Eigen::ArrayXXd const& p,Eigen::ArrayXd const& prop,
                     Eigen::VectorXd const& freqDonnees,double const& rpi,int const& numComp);
 * \brief calcul de la logvraisemblance dans l'algorithme EM pour modeles de melange multivarie
 * \param calculInter calcul inermediaire effectue dans l'etape M ( \see MmelangeMult)
 * \param exposantPi G(x,y,mu)( \see comparaison) pour chaque rang x des donn�es, chaque y de liste Sigma( \see listeSigma), et tous les mu possibles ( \see listeMu)
 * \param exposantUnMoinsPi A(x,y)-G(x,y,mu)( \see comparaison) pour chaque rang x des donn�es, chaque y de liste Sigma( \see listeSigma), et tous les mu possibles ( \see listeMu)
 * \param p vecteur probabilit� d'une bonne comparaison pour chaque classe
 * \param prop vecteur des proportions de chaque classe
 * \param freqDonnees frequence des rangs des donn�es
 * \param rpi m1!/2*...*md!/2 (pour la constante de normalisation)
 * \param numComp indice de la composante actuelle dans l'etape M \see MmelangeMult
 * \return logvraisemblance du mod�le
 *
 * cette fonction est propre a l'algorithme et ne peut etre utilisee en dehors
 *
 */
double loglikMelMult(Eigen::ArrayXd calculInter,std::vector<std::vector<Eigen::ArrayXXd> > const& exposantPi,
                     std::vector<std::vector<Eigen::ArrayXXd> > const& exposantUnMoinsPi,
                     std::vector<int> const& indiceMu,Eigen::ArrayXXd const& p,Eigen::ArrayXd const& prop,
                     Eigen::VectorXd const& freqDonnees,double const& rpi,int const& numComp);

/**
 *---------------ALGO EM mod�le de m�lange multi dim----------------
 * \fn std::tuple<Eigen::ArrayXXd,Eigen::ArrayXd,std::vector<std::vector<std::vector<int> > >,Eigen::ArrayXXd>
melangeMulti(std::vector<std::vector<std::vector<int> > > const& donnees,std::vector<int> const& frequence,int const& g,
             double const& eps=0.01,int const& maxIt=40);
 * \brief Algorithme EM pour modeles de melange multivarie
 * \param donnees vecteur o� chaque composante est un tableau de donnees de rang import�es( \see import)
 * \param frequence frequence des donn�es
 * \param g nombre de groupe
 * \param eps seuil pour crit�re d'arr�t sur la vraisemblance
 * \param maxIt nombre maximal d'it�ration pour un EM sur un mu
 * \return un tuple contenant prop,p,mu,les tik
 *
 */
std::tuple<Eigen::ArrayXXd,Eigen::ArrayXd,std::vector<std::vector<std::vector<int> > >,Eigen::ArrayXXd,double,double,double,Eigen::ArrayXd,Eigen::ArrayXd,
Eigen::ArrayXd,Eigen::ArrayXXd,std::vector<std::vector<std::vector<int> > > >
EMmelMulti(std::vector<std::vector<std::vector<int> > > const& donnees,std::vector<int> const& frequence,int const& g,
             int const& maxIt=40,double const& eps=0.0001,bool const& detail=false);


#endif // EMMEL_H_INCLUDED
