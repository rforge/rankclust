#ifndef RUN_H_
#define RUN_H_

#include "RankCluster.h"
#include <RcppEigen.h>

RcppExport SEXP semR(SEXP X,SEXP m,SEXP K,SEXP Qsem,SEXP Bsem,SEXP Ql,SEXP Bl,SEXP RjSE,SEXP RjM,SEXP maxTry,SEXP run,SEXP detail);

#endif /* RUN_H_ */
