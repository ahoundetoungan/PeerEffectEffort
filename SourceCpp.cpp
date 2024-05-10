// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;
using namespace std;


// Function for foptim
//[[Rcpp::export]]
List fsigma2eps(const double& tau,
                const double& rho,
                const arma::vec& resids,
                List& W,
                List& WW,
                List& I,
                const arma::vec& nst,
                const int& sumn,
                const int& S){
  int n1(0), n2(-1);
  double sig2eps(0);
  List lOm(S);
  for(int s(0); s < S; ++ s){
    n1            = n2 + 1;
    n2            = n1 + nst(s) - 1;
    arma::vec es  = resids.subvec(n1, n2);
    arma::mat Ws  = W[s];
    arma::mat WWs = WW[s];
    arma::mat Is  = I[s];
    arma::mat Om  = Is + tau*tau*WWs + rho*tau*(Ws + arma::trans(Ws));
    lOm[s]        = Om;
    sig2eps      += arma::cdot(es, arma::solve(Om, es));
  }
  
  return List::create(Named("se2") = sig2eps/sumn, Named("Om") = lOm);
}

//[[Rcpp::export]]
double floglike(const arma::vec& x,
                const arma::vec& resids,
                List& W,
                List& WW,
                List& I,
                const arma::vec& nst,
                const int& sumn,
                const int& S){
  double tau      = exp(x(0));
  double rho      = (exp(x(1)) - 1)/(exp(x(1)) + 1);
  //Rprintf("tau: %f\n", tau);
  //Rprintf("rho: %f\n", rho);

  double sig2eps, llh(0), val, sign;
  List lOm;
  
  {List tmp1      = fsigma2eps(tau, rho, resids, W, WW, I, nst, sumn, S);
    double tmp2   = tmp1[0];
    List tmp3     = tmp1[1];
    sig2eps       = tmp2;
    lOm           = tmp3;}
  

  //Rprintf("sigma^2_epsilon: %f\n", sig2eps);
  for(int s(0); s < S; ++ s){
    arma::mat Om  = lOm[s];
    log_det(val, sign, Om); 
    llh          += (-0.5*nst(s)*(log(sig2eps) + 1) - 0.5*(val + log(sign)));
  }
  //Rprintf("log(likelihood): %f\n", llh);
  //Rprintf("*****************\n");
  return -llh;
}

//[[Rcpp::export]]
List fdataFs(List& J){
  int S           = J.length();
  arma::vec eigval;
  arma::uvec ind;
  arma::mat eigvec;
  List F(S);
  for(int s(0); s < S; ++ s){
    //Rprintf("Find Eigenvalues and Eigenvectors: %i/%i\n", s + 1, S);
    arma::mat Js = J[s];
    arma::eig_sym(eigval, eigvec, Js, "std");
    ind          = arma::find(eigval > 0.999);
    eigvec       = eigvec.cols(ind);
    F[s]         = eigvec.t();
  }
  return F;
}


//[[Rcpp::export]]
List ftoolsml(const arma::vec& resids, 
              List& F,
              List& network,
              const double& lambda,
              const int& S){
  arma::vec es, n(S);
  arma::uvec ind;
  arma::mat Is, Ws, WWs, As;
  int n1(0), n2(-1), ns, nss;
  List I(S), W(S), WW(S), lresid(S);
  for(int s(0); s < S; ++ s){
    arma::mat Gs = network[s];
    ns           = Gs.n_rows;
    As           = arma::eye(ns, ns) - lambda*Gs;
    arma::mat Fs = F[s];
    nss          = Fs.n_rows;
    Is           = arma::eye(nss, nss);
    Ws           = Fs*As;
    WWs          = Ws*Ws.t();
    Ws           = Ws*Fs.t();
    n1           = n2 + 1;
    n2           = n1 + Fs.n_rows - 1;
    //cout<<n2<<endl;
    es           = resids.subvec(n1, n2);
    n(s)         = nss;
    I[s]         = Is;
    W[s]         = Ws;
    WW[s]        = WWs;
    lresid[s]    = es;
  }
  return List::create(Named("I") = I, Named("W") = W, Named("WW") = WW, Named("lresid") = lresid, Named("n") = n);
}

// // Function for foptimParallel
// // This function computes the elements of Omega and sigma^2_epsilon
// //[[Rcpp::export]]
// List f1s(const double& tau,
//          const double& rho,
//          const arma::vec& es, 
//          const arma::mat& Ws, 
//          const arma::mat& WWs, 
//          const arma::mat& Is){
//   arma::mat Oms = Is + tau*tau*WWs + rho*tau*(Ws + arma::trans(Ws));
//   double se2s   = arma::cdot(es, arma::solve(Oms, es));
//   return List::create(Named("se2s") = se2s, Named("Oms") = Oms);  
// }
// 
// //[[Rcpp::export]]
// double f1sig2(const double& tau,
//               const double& rho,
//               const arma::vec& es, 
//               const arma::mat& Ws, 
//               const arma::mat& WWs, 
//               const arma::mat& Is){
//   arma::mat Oms = Is + tau*tau*WWs + rho*tau*(Ws + arma::trans(Ws));
//   double se2s   = arma::cdot(es, arma::solve(Oms, es));
//   return se2s;  
// }
// 
// // This function computes the elements of Omega and sigma^2_epsilon
// //[[Rcpp::export]]
// double f2s(const double& sig2eps,
//            const arma::mat& Oms, 
//            const int& ns){
//   double val, sign;
//   log_det(val, sign, Oms); 
//   return -0.5*ns*(log(sig2eps) + 1) - 0.5*(val + log(sign));  
// }

