// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(cpp11)]]
#include <string>
#include <RcppArmadillo.h>
#include <limits>

void UpdateCentroids(
    arma::mat &Data,
    arma::mat &Cluster,
    arma::mat &Centroids,
    int &m
    ) {
    size_t n = Data.n_rows;
    size_t j = Data.n_cols;
    size_t C = Cluster.n_cols;
    arma::Col<double> newcents(n,arma::fill::zeros);
    double tw=0;
    for(size_t k = 0; k < C; k++) {
        tw=0;
        for(size_t g = 0; g < j; g++) {
            newcents += pow(Cluster(g,k),m) * Data.col(g); 
            tw += pow(Cluster(g,k),m); 
        } //   Data %*% Clusterk
        Centroids.col(k) = newcents/tw;
        for(size_t s = 0; s < n; s++) newcents(s) = 0;
    }
}

void UpdateCluster(
    arma::mat &Data,
    arma::mat &Cluster,
    arma::mat &Centroids,
    int &m,
    bool &going,
    double &thresh
    ) {
    size_t n = Data.n_rows;
    size_t j = Data.n_cols;
    size_t C = Centroids.n_cols;
    arma::mat Dist(j,C, arma::fill::zeros);
    arma::Row<double> thisclust(C, arma::fill::zeros);
    bool keep_going = false;
    double sum_denom=0;
    for(size_t g = 0; g < j; g++) {
        for(size_t k = 0; k < C; k++) {
            for(size_t nn = 0; nn < n; nn++) {
                Dist(g,k) += pow(Data(nn,g)-Centroids(nn,k),2);
            }
            Dist(g,k) = pow(Dist(g,k),0.5);
        }
    }
    for(size_t g = 0; g < j; g++) {
        for (size_t k = 0; k < C; k++) thisclust(k) = 0;
        if (m==1){
            size_t which_min_dist = 0;
            double min_dist = Dist(g,0);
            thisclust(0) = 1;
            for (size_t k = 1; k < C; k++) {
                if(Dist(g,k)<min_dist){
                    min_dist<-Dist(g,k);
                    thisclust(which_min_dist) = 0;
                    thisclust(k) = 1;
                    which_min_dist = k;
                }
            }
            for (size_t k = 0; k < C; k++) {
                if(abs(thisclust(k)-Cluster(g,k))>thresh) keep_going = true;
            }
        }else{
            for(size_t k = 0; k < C; k++) {
                sum_denom=0;
                for(size_t k2 = 0; k2 < C; k2++) {
                    if(k2==k){
                      sum_denom += 1;
                    }else{
                      sum_denom += pow(Dist(g,k)/Dist(g,k2),2/(m-1));
                    }
                }
                thisclust(k) = 1/sum_denom;
                if(abs(thisclust(k)-Cluster(g,k))>thresh) keep_going = true;
            }
        }
        Cluster.row(g) = thisclust;
    }
    going = keep_going;
}


void fca(arma::mat &Data,
    arma::mat &Cluster,
    arma::mat &Centroids,
    int &m,
    double &thresh,
    bool &usecentroids
    ){
      size_t n = Data.n_rows;
      size_t j = Data.n_cols;
      size_t C = Centroids.n_cols;
      arma::mat newCentroids(n,C, arma::fill::zeros);
      bool going = true;
      if (usecentroids){
        for(size_t nn = 0; nn < n; nn++) {
          for(size_t k = 0; k < C; k++) {
            newCentroids(nn,k) += Centroids(nn,k);
          }
        }
      }else{
        arma::uvec samples(j,arma::fill::zeros);
        samples = arma::randperm(j);
        for(size_t cl = 0; cl < C; cl++) {
          newCentroids.col(cl) = Data.col(samples(cl)); 
        }
      }
      UpdateCluster(Data, Cluster, newCentroids, m, going, thresh);
      int count = 0;
      while (going){
        UpdateCentroids(Data,Cluster,newCentroids,m);
        UpdateCluster(Data,Cluster,newCentroids,m,going,thresh);
        if (count>10000) going = false;
        count++;
      }
      if (!usecentroids){
        for(size_t nn = 0; nn < n; nn++) {
          for(size_t k = 0; k < C; k++) {
            Centroids(nn,k) = newCentroids(nn,k);
          }
        }
      }
    }         

// [[Rcpp::export]]
Rcpp::List FuzzyClusterCpp(
    arma::cube Exp,
    arma::mat Medx, 
    int C, int m,
    double thresh
    ) {
      // Initialization
      size_t n = Exp.n_rows;
      size_t j = Exp.n_cols;
      size_t r = Exp.n_slices;
      if(C<=1) C = 2;
      if(C>=j) C = j-1;
      arma::mat Centroids(n,C, arma::fill::zeros);
      arma::mat Cluster(j,C, arma::fill::zeros);
      arma::cube Fuzzy(j, C, r, arma::fill::zeros);
      bool usecents = false;
      //
      fca(Medx, Cluster, Centroids, m, thresh, usecents);//initiates clusters and centroids
      usecents = true;
      for(size_t s = 0; s < r; s++) {
        fca(Exp.slice(s), Cluster, Centroids ,m,thresh, usecents);//update centroids internally and Clusters globally
        Fuzzy.slice(s) = Cluster;
      }
      Rcpp::NumericVector FuzzyRcpp(Rcpp::wrap(Fuzzy));
      FuzzyRcpp.attr("dim") = Rcpp::IntegerVector::create(
        Rcpp::_["j"] = Fuzzy.n_rows,
        Rcpp::_["C"] = Fuzzy.n_cols,
        Rcpp::_["r"] = Fuzzy.n_slices
      );
      return Rcpp::List::create(FuzzyRcpp,Centroids,Cluster);
    }



