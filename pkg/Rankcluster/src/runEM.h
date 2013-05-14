#ifndef RUNEM_H_INCLUDED
#define RUNEM_H_INCLUDED

#include <RcppEigen.h>


RcppExport SEXP EMmelR(SEXP X,SEXP freq,SEXP m,SEXP g,SEXP eps,SEXP maxIt, SEXP run,SEXP detail);
RcppExport SEXP EMmelmultR(SEXP X,SEXP freq,SEXP m,SEXP g,SEXP eps,SEXP maxIt,SEXP run, SEXP detail);


#endif // RUNEM_H_INCLUDED
