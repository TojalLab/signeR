
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(cpp11)]]
#include <string>
#include <RcppArmadillo.h>
#include <limits>

// cube sum along rows
arma::mat cube_sum_i(const arma::cube &c) {
    arma::mat ss(c.n_cols, c.n_slices, arma::fill::zeros);
    for(size_t i = 0; i < c.n_rows; i++) {
    for(size_t j = 0; j < c.n_cols; j++) {
    for(size_t k = 0; k < c.n_slices; k++) {
        ss(j,k) += c(i,j,k);
    }}}
    return ss;
}

// cube sum along columns
arma::mat cube_sum_j(const arma::cube &c) {
    arma::mat ss(c.n_rows, c.n_slices, arma::fill::zeros);
    for(size_t i = 0; i < c.n_rows; i++) {
    for(size_t j = 0; j < c.n_cols; j++) {
    for(size_t k = 0; k < c.n_slices; k++) {
        ss(i,k) += c(i,j,k);
    }}}
    return ss;
}

// single vector multinomial draw, uses R function
arma::Col<int> one_multinom(int size, arma::rowvec prob) {
    int k = prob.size();
    arma::Col<int> v_ans(k);
    R::rmultinom(size, prob.begin(), k, v_ans.begin());
    return v_ans;
}

// draw a value from a gamma distribution
inline double one_gamma_dist(double shape, double rate) {
    shape = std::max(pow(10,-160),shape);
    double scale = std::max(pow(10,-160),1/rate);
    return R::rgamma(shape, scale);
}

// turns 3D/4D index into 1D
inline size_t adr(const size_t * dim, size_t i=0, size_t j=0, size_t k=0, size_t l=0) {
    return i + dim[0]*(j + dim[1]*(k + l*dim[2]));
}

// creates a array with given dims, use adr() function for element access
Rcpp::NumericVector create_array(Rcpp::IntegerVector dims) {
    size_t k = 1;
    for(int i = 0; i < dims.size(); i++) k *= dims[i];
    Rcpp::NumericVector array(k);
    array.attr("dim") = dims;
    return array;
}

// copy the cube values to a 4D array cube-slice, like: dest[,,,l] <- src
inline void copy_cube_to_4Darray(const arma::cube &src, const size_t l,
    Rcpp::NumericVector &dest
    ) {
    const size_t dim[3] = {src.n_rows, src.n_cols, src.n_slices};
    for(size_t i = 0; i < src.n_rows; i++) {
    for(size_t j = 0; j < src.n_cols; j++) {
    for(size_t k = 0; k < src.n_slices; k++) {
        dest[adr(dim, i,j,k,l)] = src(i,j,k);
    }}}
}

// converts a arma cube to rcpp and add dim attribs
Rcpp::NumericVector cube_to_rcpp(const arma::cube &c,
    const std::vector<std::string> &names
    ) {
    Rcpp::NumericVector v(Rcpp::wrap(c));
    v.attr("dim") = Rcpp::IntegerVector::create(
        Rcpp::_[names[0]] = c.n_rows,
        Rcpp::_[names[1]] = c.n_cols,
        Rcpp::_[names[2]] = c.n_slices
    );
    return v;
}

// step 1: update Z
void gibbs_step1(
    const arma::mat &M, const arma::mat &P, const arma::mat &E,
    arma::cube &Z, arma::cube &Fi
    ) {
    arma::mat PE = P * E;
    for(size_t m = 0; m < Z.n_slices; m++) {
        Fi.slice(m) = (P.col(m) * E.row(m)) / PE;
    }
    arma::rowvec v(Z.n_slices, arma::fill::zeros);
    for(size_t s = 0; s < Z.n_rows; s++) {
        for(size_t g = 0; g < Z.n_cols; g++) {
            for(size_t m = 0; m < Z.n_slices; m++) { v[m] = Fi(s,g,m); }
            arma::Col<int> u = one_multinom(M(s,g), v);
            for(size_t m = 0; m < Z.n_slices; m++) {
                Z(s,g,m) = u[m];
            }
        }
    }
}

// step 2: update P
void gibbs_step2(
    const arma::mat &W, const arma::cube &Z, const arma::mat &E,
    const arma::mat &Ap, const arma::mat &Bp,
    arma::mat &P
    ) {
    arma::mat corrE = W * E.t();
    arma::mat Zin = cube_sum_j(Z);
    for(size_t s = 0; s < Z.n_rows; s++) {
        for(size_t m = 0; m < Z.n_slices; m++) {
            double shape = Ap(s,m) + 1 + Zin(s,m);
            double rate = Bp(s,m) + corrE(s,m);
            P(s,m) = one_gamma_dist(shape, rate);
        }
    }
}

// step 3: update E
void gibbs_step3(
    const arma::mat &W, const arma::cube &Z, const arma::mat &P,
    const arma::mat &Ae, const arma::mat &Be,
    arma::mat &E
    ) {
    arma::mat corrP = P.t() * W;
    arma::mat Znj = cube_sum_i(Z).t();
    for(size_t m = 0; m < Z.n_slices; m++) {
        for(size_t g = 0; g < Z.n_cols; g++) {
            double shape = Ae(m,g) + 1 + Znj(m,g);
            double rate = Be(m,g) + corrP(m,g);
            E(m,g) = one_gamma_dist(shape, rate);
        }
    }
}

// step 4: update Bp
void gibbs_step4(
    const arma::cube &Z, const arma::mat &P, const arma::mat &Ap,
    const double ap, const double bp,
    const double var_ap,
    arma::mat &Bp
    ) {
    for(size_t s = 0; s < Z.n_rows; s++) {
        for(size_t m = 0; m < Z.n_slices; m++) {
            double shape = Ap(s,m) + ap + 1;
            double rate = P(s,m) + bp;
            Bp(s,m) = std::max(pow(10,-160), one_gamma_dist(shape,rate));
        }
    }
}

// step 5: update Be
void gibbs_step5(
    const arma::cube &Z, const arma::mat &E, const arma::mat &Ae,
    const double ae, const double be,
    const double var_ae,
    arma::mat &Be
    ) {
    for(size_t m = 0; m < Z.n_slices; m++) {
        for(size_t g = 0; g < Z.n_cols; g++) {
            double shape = Ae(m,g) + ae + 1;
            double rate = E(m,g) + be;
            Be(m,g) = std::max(pow(10,-160), one_gamma_dist(shape,rate));
        }
    }
}

// step 6: update Ap
void gibbs_step6(
    const arma::cube &Z, const arma::mat &P, const arma::mat Bp,
    const double lp, const double var_ap,
    arma::mat &Ap
    ) {
    arma::mat Y(Z.n_rows, Z.n_slices, arma::fill::zeros);
    for(size_t s = 0; s < Z.n_rows; s++) {
        for(size_t m = 0; m < Z.n_slices; m++) {
            double shape = pow(Ap(s,m), 2) / var_ap;
            double rate = Ap(s,m) / var_ap;
            Y(s,m) = std::max(pow(10,-160), one_gamma_dist(shape,rate));
        }
    }

    arma::mat LnProb = (Y-Ap) % (log(Bp)+log(P)-lp) +
        arma::lgamma(Ap+1) - arma::lgamma(Y+1) +
        arma::lgamma(pow(Ap,2)/var_ap) - arma::lgamma(pow(Y,2)/var_ap) +
        ((pow(Y,2) - pow(Ap,2)) / var_ap) % log((Ap % Y)/var_ap) +
        log(Y) - log(Ap);

    arma::mat ProbCh = exp(LnProb);
    arma::mat U = arma::randu(Z.n_rows, Z.n_slices);

    for(size_t i = 0; i < ProbCh.n_rows; i++) {
        for(size_t j = 0; j < ProbCh.n_cols; j++) {
            if(P(i,j) == 0) {
                ProbCh(i,j) = Y(i,j) < Ap(i,j);
            }
        }
    }

    for(size_t i = 0; i < Ap.n_rows; i++) {
        for(size_t j = 0; j < Ap.n_cols; j++) {
            if(U(i,j) <= ProbCh(i,j)) {
                Ap(i,j) = Y(i,j);
            }
        }
    }
}

// step 7: update Ae
void gibbs_step7(
    const arma::cube &Z, const arma::mat &E, const arma::mat &Be,
    const double le, const double var_ae,
    arma::mat &Ae
    ) {
    arma::mat Y(Z.n_slices, Z.n_cols, arma::fill::zeros);
    for(size_t m = 0; m < Z.n_slices; m++) {
        for(size_t g = 0; g < Z.n_cols; g++) {
            double shape = pow(Ae(m,g), 2) / var_ae;
            double rate = Ae(m,g) / var_ae;
            Y(m,g) = std::max(pow(10,-160), one_gamma_dist(shape,rate));
        }
    }

    arma::mat LnProb = (Y - Ae) % (log(Be) + log(E) - le) +
        arma::lgamma(Ae + 1) - arma::lgamma(Y + 1) +
        arma::lgamma(pow(Ae,2)/var_ae) - arma::lgamma(pow(Y,2)/var_ae) +
        ((pow(Y,2) - pow(Ae,2))/var_ae) % log((Ae % Y)/var_ae) +
        log(Y) - log(Ae);

    arma::mat ProbCh = exp(LnProb);
    arma::mat U = arma::randu(Z.n_slices, Z.n_cols);

    for(size_t i = 0; i < ProbCh.n_rows; i++) {
        for(size_t j = 0; j < ProbCh.n_cols; j++) {
            if(E(i,j) == 0) {
                ProbCh(i,j) = Y(i,j) < Ae(i,j);
            }
        }
    }

    for(size_t i = 0; i < Ae.n_rows; i++) {
        for(size_t j = 0; j < Ae.n_cols; j++) {
            if(U(i,j) <= ProbCh(i,j)) {
                Ae(i,j) = Y(i,j);
            }
        }
    }
}

// [[Rcpp::export]]
Rcpp::List GibbsSamplerCpp(
    arma::mat M, arma::mat W, arma::cube Z, arma::mat P,
    arma::mat E, arma::mat Ap, arma::mat Bp, arma::mat Ae, arma::mat Be,
    double ap, double bp, double ae, double be, double lp, double le,
    double var_ap,
    double var_ae,
    int burn, int eval,
    bool Pfixed,
    bool Zfixed, bool Thetafixed, bool Afixed, bool keep_par
    ) {

    bool Bfixed = false;
    if(Zfixed) Thetafixed = true;
    if(Thetafixed) Afixed = Bfixed = true;

    //Initialize variables:
    size_t i = Z.n_rows;
    size_t j = Z.n_cols;
    size_t n = Z.n_slices;
    arma::cube Fi(i, j, n, arma::fill::zeros);

    //Initialize max_mlhood_num as -Infinity:
    double max_mlhood_num = std::numeric_limits<double>::lowest();

    Rcpp::NumericVector Zs;
    if(!Zfixed && keep_par) Zs = create_array(Rcpp::IntegerVector::create(
        Rcpp::_["i"] = i, Rcpp::_["j"] = j,
        Rcpp::_["n"] = n, Rcpp::_["r"] = 1
    ));

    Rcpp::NumericVector Fis;
    if(!Zfixed && Thetafixed) Fis = create_array(Rcpp::IntegerVector::create(
        Rcpp::_["i"] = i, Rcpp::_["j"] = j,
        Rcpp::_["n"] = n, Rcpp::_["r"] = eval
    ));
    arma::cube Ps(i, n, eval, arma::fill::zeros);
    arma::cube Es(n, j, eval, arma::fill::zeros);
    arma::cube Aps;
    arma::cube Aes;
    if(!Afixed && keep_par) {
        Aps = arma::cube(i, n, eval, arma::fill::zeros);
        Aes = arma::cube(n, j, eval, arma::fill::zeros);
    }
    arma::cube Bps;
    arma::cube Bes;
    if(!Bfixed) {
        Bps = arma::cube(i, n, eval, arma::fill::zeros);
        Bes = arma::cube(n, j, eval, arma::fill::zeros);
    }
    arma::Row<int> flist = {1, 2, 3, 4, 5, 6, 7};

    for(int k = 0; k < burn+eval; k++) {
        flist = arma::shuffle(flist);
        for(size_t f = 0; f < flist.size(); f++) {
            int step = flist[f];
            switch(step) {
                case 1: if(!Zfixed) gibbs_step1(M,P,E, Z,Fi); break;
                case 2: if(!Pfixed) gibbs_step2(W,Z,E,Ap,Bp, P); break;
                case 3: gibbs_step3(W,Z,P,Ae,Be, E); break;
                case 4: if(!Pfixed && !Bfixed) gibbs_step4(Z,P,Ap,ap,bp,var_ap, Bp); break;
                case 5: if(!Bfixed) gibbs_step5(Z,E,Ae,ae,be,var_ae, Be); break;
                case 6: if(!Pfixed && !Afixed) gibbs_step6(Z,P,Bp,lp,var_ap, Ap); break;
                case 7: if(!Afixed) gibbs_step7(Z,E,Be,le,var_ae, Ae); break;
            }
        }

        if(k >= burn) {
            // save samples
            int k2 = k - burn;
            if(!Zfixed && keep_par) copy_cube_to_4Darray(Z, 0, Zs);
            if(!Zfixed && Thetafixed) copy_cube_to_4Darray(Fi, k2, Fis);
            Ps.slice(k2) = P;
            Es.slice(k2) = E;
            if(!Afixed && keep_par) {
                Aps.slice(k2) = Ap;
                Aes.slice(k2) = Ae;
            }
            if(!Bfixed) {
                Bps.slice(k2) = Bp;
                Bes.slice(k2) = Be;
            }

            //Compute mlhood_numerator
            //LogPosteriorZ=log(P(Z,M|E,P,Theta,Omega))
            arma::cube EPz(Z.n_rows, Z.n_cols, Z.n_slices, arma::fill::zeros);
            for (size_t kk=0; kk<n; kk++){
                EPz.slice(kk) = (P.col(kk) * E.row(kk)) % W;
                //EPz keeps expected values for Z, given P and E
            }
            double logPosteriorZ = arma::accu(Z%log(EPz)-EPz -arma::lgamma(Z+1));
            EPz.clear();
            //LogPriorP=log(P(Phat|Thetahat))
            double logPriorP = arma::accu( (Ap+1)%log(Bp) - arma::lgamma(Ap+1) +
                Ap%log(P) -(Bp)%P );

            //LogPriorE=log(P(Ehat|Thetahat))
            double logPriorE = arma::accu( (Ae+1)%log(Be) - arma::lgamma(Ae+1) +
                Ae%log(E) -(Be)%E );

            //LogPriorAp=log(P(Aphat|Omega))
            double logPriorAp= i*n*log(lp) -lp*arma::accu(Ap);
            //LogPriorBp=log(P(Bphat|Omega))
            double logPriorBp= i*n*(ap*log(bp) -lgamma(ap)) +
                (ap-1)*arma::accu(log(Bp)) -bp*arma::accu(Bp);

            //LogPriorAe=log(P(Aehat|Omega))
            double logPriorAe= n*j*log(le) -le*arma::accu(Ae);
            //LogPriorBe=log(P(Behat|Omega))
            double logPriorBe= n*j*(ae*log(be) -lgamma(ae)) +
                (ae-1)*arma::accu(log(Be)) -be*arma::accu(Be);

            //Likelihood numerator
            double mlhood_num = logPosteriorZ+logPriorP+logPriorE+logPriorAp+
                logPriorBp+logPriorAe+logPriorBe;

            //Save instance if mlhood_num is maximal
            if(mlhood_num > max_mlhood_num){
                max_mlhood_num = mlhood_num;
            }
        }
    }

    if (Zfixed) {
        return Rcpp::List::create(NA_LOGICAL,NA_LOGICAL,
            cube_to_rcpp(Ps, {"i","n","r"}),
            cube_to_rcpp(Es, {"n","j","r"}));
    }
    else if (Thetafixed) {
        arma::cube meanlogFis(i, j, n, arma::fill::zeros);
        double logsum=0;
        const size_t dimsl[3] = {i, j, n};

        //meanlogFis = mean.tensor(log(Fis),along='r'):
        for(size_t s=0; s<i ; s++) {
            for(size_t g=0; g<j ; g++) {
                for(size_t m=0; m<n ; m++) {
                    logsum = 0;
                    for(int r = 0; r < eval; r++) {
                        logsum += log(Fis[adr(dimsl, s,g,m,r)]);
                    }
                    meanlogFis(s,g,m) = logsum/eval;
                }
            }
        }
        return Rcpp::List::create(NA_LOGICAL,
            cube_to_rcpp(meanlogFis, {"i","j","n"}));
    }
    else if (Afixed) {
        return Rcpp::List::create(NA_LOGICAL,NA_LOGICAL,
            cube_to_rcpp(Ps, {"i","n","r"}),
            cube_to_rcpp(Es, {"n","j","r"}));
    }
    else if (keep_par) {
        return Rcpp::List::create(Zs,NA_LOGICAL,
            cube_to_rcpp(Ps, {"i","n","r"}),
            cube_to_rcpp(Es, {"n","j","r"}),
            cube_to_rcpp(Aps, {"i","n","r"}),
            cube_to_rcpp(Bps, {"i","n","r"}),
            cube_to_rcpp(Aes, {"n","j","r"}),
            cube_to_rcpp(Bes, {"n","j","r"}));
    }
    else {
        return Rcpp::List::create(NA_LOGICAL,NA_LOGICAL,
            cube_to_rcpp(Ps, {"i","n","r"}),
            cube_to_rcpp(Es, {"n","j","r"}),
            NA_LOGICAL,
            cube_to_rcpp(Bps, {"i","n","r"}),
            NA_LOGICAL,
            cube_to_rcpp(Bes, {"n","j","r"}));
    }
}


