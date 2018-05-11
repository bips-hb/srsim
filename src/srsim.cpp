#ifdef _OPENMP
  #include <omp.h>
#endif

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

//' Valid Report
//' 
//' Checks whether the given report is 'valid'. A report should contain at least
//' one drug and at least one event (otherwise it would never been sent to the 
//' spontaneous reporitng sytem anyway).
//' 
//' @param report A logical matrix with one row and \code{n_drugs + n_events} columns
//' @param n_drugs The total number of drugs 
//' @param n_events The total number of events 
//' 
//' @return \code{TRUE} when the report is valid, \code{FALSE} otherwise
//' @export
// [[Rcpp::export]] 
bool validReport (Rcpp::LogicalMatrix report, int n_drugs, int n_events) {
  bool valid = false ; 
  for (int i = 0; i < n_drugs; i ++) {
    if (report(0, i)) {
      valid = true ; 
    } 
  }
  
  if (!valid) {
    return false ; 
  }

  for (int i = n_drugs; i < n_drugs + n_events; i ++) {
    if (report(0, i)) {
      return true ; 
    } 
  }

  return false ; 
}

//' Correlations and Marginal Probabilities 
//' 
//' Returns \code{TRUE} when two binary variables with the given marginal probabilities
//' \code{p_x} and \code{p_y} can have a correlation of \code{rho}. Otherwise it returns 
//' \code{FALSE}. 
//' 
//' Suppose two binary variables, \eqn{X} and \eqn{Y}, have the marginal probabilities
//' \eqn{p(X = 1) = p_X} and \eqn{p(Y = 1) = p_Y}. Their common probability is equal to 
//' \eqn{p(X = 1, Y = 1) = p_{XY}}. Their correlation, \eqn{\rho}, can then be written as
//' \deqn{\rho = \frac{p_{XY} - p_{X}p_{Y}}{\sqrt{p_{X} (1 - p_{X}) p_{Y} (1 - p(Y))}}}
//' One can see that not every correlation is admissable given the marginal probabilities. 
//' In fact, the correlation must lie within the interval
//' \deqn{max(-\alpha_X \alpha_Y, -1/(\alpha_X \alpha_Y)) \le \rho \le min(\alpha_i/\alpha_j, \alpha_j / \alpha_i)}
//' where \eqn{\alpha_X = \sqrt{p_X (1 - p_X)}} and \eqn{\alpha_Y = \sqrt{p_Y (1 - p_Y)}}.
//' These constraints are a result of the Fréchet inequalities. 
//' 
//' @param rho The correlation (\eqn{\rho})
//' @param p_x The probability \eqn{p(X = 1)}
//' @param p_y The probability \eqn{p(Y = 1)}
//' 
//' @return \code{TRUE} when the two binary variables can exhibit the given correlation, \code{FALSE} otherwise
//' @export
// [[Rcpp::export]]
bool validCorrelation (double rho, double p_x, double p_y) {
  
  // first check the upper bound 
  double constraint = sqrt((p_x * (1 - p_y)) / ((1 - p_x) * p_y)) ; 
  if (constraint > sqrt(((1 - p_x) * p_y) / (p_x * (1 - p_y)))) { 
    constraint = sqrt(((1 - p_x) * p_y) / (p_x * (1 - p_y))) ; 
  }
  
  if (rho > constraint) { 
    return false ; 
  }
  
  // check the lower bound 
  constraint = -1.0 * sqrt((p_x * p_y) / ((1 - p_x) * (1 - p_y))) ; 
  
  if (constraint < -1.0 * sqrt(((1 - p_x) * (1 - p_y))/(p_x * p_y))) {
    constraint = -1.0 * sqrt(((1 - p_x) * (1 - p_y))/(p_x * p_y)) ; 
  }
  
  if (rho < constraint) {
    return false ;
  }
  
  return true ; 
}

//' Randows Rows 
//' 
//' Returns a random sample of rows from an integer matrix (without replacement). 
//' 
//' @section Warning:
//' The matrix can contain integers only. 
//' 
//' @param mat A integer matrix 
//' @param nrows The number of rows of the matrix \code{mat}
//' @param n_samples The number of rows one wants to sample 
//' 
//' @return An integer matrix with \code{n_samples} rows
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix sampleRandomRows (Rcpp::IntegerMatrix mat, int nrows, int n_samples) {

  // create matrix that will contain the randomly samples row  
  Rcpp::IntegerMatrix random_sample (n_samples, mat.ncol()); 
  
  // create a index vector
  Rcpp::IntegerVector rows(nrows) ;
  for (int i = 0; i < nrows; i ++) {
    rows[i] = i ; 
  }
  
  // randomly sample n_samples
  rows = RcppArmadillo::sample(rows, n_samples, false, NumericVector::create()) ; 
  
  for (int i = 0 ; i < n_samples; i ++) {
    random_sample(i, _ ) = mat(rows[i], _ ) ; 
  }
  
  return random_sample; 
} 

//' Random Drug-Event Pairs 
//' 
//' Returns a number of randomly select drug-event pairs that can exhibit a desired
//' correlation (\code{rho}). 
//' 
//' @param prob_drugs List with the marginal probabilities of the drugs 
//' @param prob_events List with the marginal probabilities of the events
//' @param n_wanted_pairs The number of randomly selected pairs 
//' @param rho The correlation one wants these pairs to have
//' 
//' @return A matrix with \code{n_wanted_pairs} rows and two column. Each row are the 
//'         indices of a drug-event pair (the first index is for the drug, the second
//'         for the event)
//' 
//' @seealso \code{\link{validCorrelation}}, \code{\link{sampleRandomRows}}, 
//'          \code{\link{returnRandomPairsRcpp}}
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix returnRandomDrugEventPairsRcpp (Rcpp::NumericVector prob_drugs, Rcpp::NumericVector prob_events, int n_wanted_pairs, double rho) {
  int n_drugs  = prob_drugs.size() ;
  int n_events = prob_events.size() ;
  int n        = n_drugs + n_events ;
  
  // will keep the total number of pairs 
  Rcpp::IntegerMatrix pairs (n*n, 2) ; 
  int n_pairs = 0 ; // keeps track of the number of pairs found so far
  
  // go over all possible pairs
  for (int i = 0 ; i < n_drugs; i ++) {
    for (int j = 0; j < n_events; j ++) {
      if (validCorrelation(rho, prob_drugs[i], prob_events[j])) {
        pairs(n_pairs, 0) = i ; 
        pairs(n_pairs, 1) = j ; 
        n_pairs ++ ; 
      }
    }
  }
  
  pairs = sampleRandomRows(pairs, n_pairs, n_wanted_pairs) ; 
  
  return(pairs); 
}

//' Random Pairs 
//' 
//' Returns a number of randomly select (drug-drug or event-event) pairs that can 
//' exhibit a desired correlation (\code{rho}). 
//' 
//' @param margprob List with the marginal probabilities 
//' @param n_wanted_pairs The number of randomly selected pairs 
//' @param rho The correlation one wants these pairs to exhibit 
//' 
//' @return A matrix with \code{n_wanted_pairs} rows and two column. Each row are the 
//'         indices of a pair
//' 
//' @seealso  \code{\link{validCorrelation}}, \code{\link{sampleRandomRows}}, 
//'           \code{\link{returnRandomDrugEventPairsRcpp}}
// [[Rcpp::export]]
Rcpp::IntegerMatrix returnRandomPairsRcpp (Rcpp::NumericVector margprob, int n_wanted_pairs, double rho) {
  int n = margprob.size() ; 
  
  // will keep the total number of pairs 
  Rcpp::IntegerMatrix pairs (n*(n - 1)/2, 2) ; 
  int n_pairs = 0 ; // keeps track of the number of pairs found so far
  
  // go over all possible pairs
  for (int i = 0 ; i < n; i ++) {
    for (int j = i+1; j < n; j ++) {
      if (validCorrelation(rho, margprob[i], margprob[j])) {
        pairs(n_pairs, 0) = i ; 
        pairs(n_pairs, 1) = j ; 
        n_pairs ++ ; 
      }
    }
  }
  
  pairs = sampleRandomRows(pairs, n_pairs, n_wanted_pairs) ; 
  
  return(pairs); 
}

//' Correlation Matrix 
//' 
//' Generates a correlation matrix for the SR data set. It randomly selects 
//' a number of drug-event, drug-drug, and event-event pairs to exhibit a 
//' given correlation. It takes the marginal probabilities of the drugs and 
//' events into account, to make sure that the pairs can indeed show the 
//' wanted correlation (see the function \code{\link{validCorrelation}}).
//' 
//' @param prob_drugs List with the marginal probabilities of the drugs 
//' @param prob_events List with the marginal probabilities of the events
//' @param n_correlated_drugs The number of drug-drug pairs that will have correlation \code{rho_drugs}
//' @param rho_drugs The correlation for the drug-drug pairs 
//' @param n_correlated_events The number of event-event pairs that will have correlation \code{rho_events}
//' @param rho_events The correlation for the event-event pairs 
//' @param n_correlated_pairs The number of drug-event pairs that will have correlation \code{rho}
//' @param rho The correlation for the drug-event pairs 
//' 
//' @return A matrix 
//' 
//' @section Warning: 
//' The matrix that is returned is not guaranteed to be a correlation matrix. 
//' One still needs to check whether the matrix is indeed postive definite. 
//' 
//' @seealso \code{\link{returnRandomDrugEventPairsRcpp}}, 
//'          \code{\link{returnRandomDPairsRcpp}},
//'          \code{\link{validCorrelation}}
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix generateCorrelationMatrixRcpp (Rcpp::NumericVector prob_drugs, 
                                                   Rcpp::NumericVector prob_events, 
                                                   int n_correlated_drugs, 
                                                   double rho_drugs, 
                                                   int n_correlated_events,
                                                   double rho_events,
                                                   int n_correlated_pairs, 
                                                   double rho) {

  int i, j ; 
  
  int n_drugs  = prob_drugs.size() ; 
  int n_events = prob_events.size() ;
  int n        = n_drugs + n_events ; 
  
  // initialize correlation matrix
  Rcpp::NumericMatrix corrmat (n) ; 
  for (i = 0; i < n; i ++) {
    for (j = 0; j < n; j ++) {
      corrmat(i, j) = 0 ; 
    }
  }
    
  // set diagonal
  for (i = 0; i < n; i ++) {
    corrmat(i, i) = 1 ; 
  }
    
  // select random pairs 
  Rcpp::IntegerMatrix dd_pairs = returnRandomPairsRcpp(prob_drugs, n_correlated_drugs, rho_drugs) ; 
  Rcpp::IntegerMatrix ee_pairs = returnRandomPairsRcpp(prob_events, n_correlated_events, rho_events) ; 
  Rcpp::IntegerMatrix de_pairs = returnRandomDrugEventPairsRcpp(prob_drugs, prob_events, n_correlated_pairs, rho) ; 
  
  // fill in the matrix with the pairs  
  for (i = 0; i < n_correlated_drugs; i ++) {
    corrmat(dd_pairs(i, 0), dd_pairs(i, 1)) = rho_drugs ; 
    corrmat(dd_pairs(i, 1), dd_pairs(i, 0)) = rho_drugs ; 
  }
  
  for (i = 0; i < n_correlated_events; i ++) {
    corrmat(ee_pairs(i, 0) + n_drugs, ee_pairs(i, 1) + n_drugs) = rho_events ; 
    corrmat(ee_pairs(i, 1) + n_drugs, ee_pairs(i, 0) + n_drugs) = rho_events ; 
  }
  
  for (i = 0; i < n_correlated_pairs; i ++) {
    corrmat(de_pairs(i, 0), de_pairs(i, 1) + n_drugs) = rho ; 
    corrmat(de_pairs(i, 1) + n_drugs, de_pairs(i, 0)) = rho ; 
  }
  
  return corrmat ;   
}

//' Block Correlation Matrix
//' 
//' @param margprob List with the marginal probabilities
//' @param blocksize The size of the blocks
//' 
//' @return A matrix
//' 
//' @section Warning: 
//' The matrix that is returned is not guaranteed to be a correlation matrix. 
//' One still needs to check whether the matrix is indeed postive definite. 
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix generateBlockCorrelationMatrixRcpp (Rcpp::NumericVector margprob, 
                                                           int blocksize,
                                                           double rho) {
  int i, j, index1, index2; 
  int n = margprob.size() ;
  
  // initialize correlation matrix
  Rcpp::NumericMatrix corrmat (n) ; 
  for (i = 0; i < n; i ++) {
    for (j = 0; j < n; j ++) {
      corrmat(i, j) = 0 ; 
    }
  }
  
  double p_x, p_y;
  double lower_thres, upper_thres ; 
  
  // go through all the blocks
  for (int b = 0 ; b < (n / blocksize); b++) {
    for (i = 0; i < blocksize; i ++) {
      index1 = b*blocksize + i ; 
      p_x = margprob[index1] ; 
      for (j = 0; j < blocksize; j ++) {
        index2 = b*blocksize + j ; 
        p_y = margprob[index2] ; 
        upper_thres = sqrt((p_x * (1 - p_y)) / ((1 - p_x) * p_y)) ; 
        lower_thres = -1.0 * sqrt((p_x * p_y) / ((1 - p_x) * (1 - p_y))) ;
        if (rho > upper_thres) { 
          corrmat(index1, index2) = upper_thres ; 
          corrmat(index2, index1) = upper_thres ; 
        } else {
          corrmat(index1, index2) = rho ; 
          corrmat(index2, index1) = rho ;
        }
      }
    }
  }
  
  for (i = 0; i < n; i ++) {
     corrmat(i,i) = 1.0 ;  
  }
  
  return corrmat ; 
}

//' Simple Random Correlation Matrix
//' 
//' @param margprob List with the marginal probabilities
//'
//' @return A matrix 
//' 
//' @section Warning: 
//' The matrix that is returned is not guaranteed to be a correlation matrix. 
//' One still needs to check whether the matrix is indeed postive definite. 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix generateSimpleRandomCorrelationMatrix (Rcpp::NumericVector margprob) {
  int i, j; 
  int n = margprob.size() ;
  
  // initialize correlation matrix
  Rcpp::NumericMatrix corrmat (n) ; 
  for (i = 0; i < n; i ++) {
    for (j = 0; j < n; j ++) {
      corrmat(i, j) = 0 ; 
    }
  }
  
  double p_x, p_y;
  double lower_thres, upper_thres ; 
  
  // go over all drug-drug pairs and select a valid random correlation
  for (i = 0; i < n; i ++) { 
    p_x = margprob[i] ; 
    for (j = i+1; j < n; j ++) {
      p_y = margprob[j] ; 
      upper_thres = sqrt((p_x * (1 - p_y)) / ((1 - p_x) * p_y)) ; 
      lower_thres = -1.0 * sqrt((p_x * p_y) / ((1 - p_x) * (1 - p_y))) ;
      if (upper_thres > 1.0) { 
        upper_thres = 1.0 ; 
      } 
      if (lower_thres < -1.0) {
        lower_thres = -1.0 ;  
      }
      
      corrmat(i, j) = R::runif(lower_thres, upper_thres) ; 
      corrmat(j, i) = corrmat(i, j) ; 
    }
  }
  
  // set diagonal
  for (i = 0; i < n; i ++) {
    corrmat(i, i) = 1 ; 
  }
  
  return corrmat ; 
}

//' Random Correlation Matrix 
//' 
//' Generates a correlation matrix for the SR data set. It randomly selects 
//' a number of drug-event, drug-drug, and event-event pairs to exhibit a 
//' given correlation. It takes the marginal probabilities of the drugs and 
//' events into account, to make sure that the pairs can indeed show the 
//' wanted correlation (see the function \code{\link{validCorrelation}}).
//' 
//' @param prob_drugs List with the marginal probabilities of the drugs 
//' @param prob_events List with the marginal probabilities of the events
//' @param n_correlated_pairs The number of drug-event pairs that will have correlation \code{rho}
//' @param rho The correlation for the drug-event pairs 
//' 
//' @return \item{corrmat}{A (potential) correlation matrix}
//'         \item{de_pairs}{A integer matrix with the drug-event pairs}
//' 
//' @section Warning: 
//' The matrix that is returned is not guaranteed to be a correlation matrix. 
//' One still needs to check whether the matrix is indeed postive definite. 
//' 
//' @seealso \code{\link{returnRandomDrugEventPairsRcpp}}, 
//'          \code{\link{returnRandomDPairsRcpp}},
//'          \code{\link{validCorrelation}}
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix generateRandomCorrelationMatrixRcpp (Rcpp::NumericVector prob_drugs, 
                                                         Rcpp::NumericVector prob_events, 
                                                         int n_correlated_pairs, 
                                                         double rho) {
  
  int i, j; 
  
  int n_drugs  = prob_drugs.size() ; 
  int n_events = prob_events.size() ;
  int n        = n_drugs + n_events ; 
  
  // initialize correlation matrix
  Rcpp::NumericMatrix corrmat (n) ; 
  for (i = 0; i < n; i ++) {
    for (j = 0; j < n; j ++) {
      corrmat(i, j) = 0 ; 
    }
  }
  
  // select random pairs 
  Rcpp::IntegerMatrix de_pairs = returnRandomDrugEventPairsRcpp(prob_drugs, prob_events, n_correlated_pairs, rho) ; 
  
  // fill in the matrix with the pairs  
  for (i = 0; i < n_correlated_pairs; i ++) {
    corrmat(de_pairs(i, 0), de_pairs(i, 1) + n_drugs) = rho ; 
    corrmat(de_pairs(i, 1) + n_drugs, de_pairs(i, 0)) = rho ; 
  }
  
  double p_x, p_y;
  double lower_thres, upper_thres ; 
  
  // go over all drug-drug pairs and select a valid random correlation
  for (i = 0; i < n_drugs; i ++) { 
    p_x = prob_drugs[i] ; 
    for (j = i+1; j < n_drugs; j ++) {
      p_y = prob_drugs[j] ; 
      upper_thres = sqrt((p_x * (1 - p_y)) / ((1 - p_x) * p_y)) ; 
      lower_thres = -1.0 * sqrt((p_x * p_y) / ((1 - p_x) * (1 - p_y))) ;
      
      corrmat(i, j) = R::runif(lower_thres, upper_thres) ; 
      corrmat(j, i) = corrmat(i, j) ; 
    }
  }
  
  // go over all event-event pairs and select a valid random correlation
  for (i = 0; i < n_events; i ++) { 
    p_x = prob_events[i] ; 
    for (j = i+1; j < n_events; j ++) {
      p_y = prob_events[j] ; 
      upper_thres = sqrt((p_x * (1 - p_y)) / ((1 - p_x) * p_y)) ; 
      lower_thres = -1.0 * sqrt((p_x * p_y) / ((1 - p_x) * (1 - p_y))) ;
      
      corrmat(n_drugs + i, n_drugs + j) = R::runif(lower_thres, upper_thres) ; 
      corrmat(n_drugs + j, n_drugs + i) = corrmat(n_drugs + i, n_drugs + j) ; 
    }
  }
  
  // set diagonal
  for (i = 0; i < n; i ++) {
    corrmat(i, i) = 1 ; 
  }
  
  return corrmat ;   
}



//' From Correlation To Covariance Matrix
//' 
//' Converts a correlation matrix for binary variates to a covariance matrix given the
//' variables marginal probabilities. 
//' 
//' @param corrmat The correlation matrix 
//' @param margprob A list with the marginal probabilities 
//' 
//' @return The covariance matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix corr2cov (Rcpp::NumericMatrix corrmat, Rcpp::NumericVector margprob) {
  int n = margprob.size() ; 
  Rcpp::NumericMatrix covmat(n) ; 
  
  // compute the off-diagonal entries
  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j ++) {
      covmat(i,j) = corrmat(i,j) * sqrt(margprob[i] * (1 - margprob[i]) * margprob[j] * (1 - margprob[j])) ;  
      covmat(j,i) = covmat(i,j) ; 
    }
  }
  
  // compute the diagonal entries
  for (int i = 0; i < n; i ++) {
    covmat(i,i) = margprob[i] * (1 - margprob[i]) ;
  }
  
  return covmat ; 
}

//' Sampling Multivariate Normal 
//' 
//' Returns a sample of multivariate normal data. The code is based on the code 
//' from \link{http://gallery.rcpp.org/articles/simulate-multivariate-normal/} 
//' 
//' @param n The number of samples
//' @param mu Vector with the means 
//' @param L The Cholesky decomposition of the covariance matrix 
//' @param ncols The number of variates 
//' 
//' @return A matrix, where each row is a sample 
//' 
//' @examples
//' arma::mat L = arma::chol(covmat) ;   // covmat is the covariance matrix
//' arma::mat sample = mvrnormArma(n, means, L, covmat.n_cols) ; 
//'
// [[Rcpp::export]] 
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat L, int ncols) {
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * L;
}


//' Larger than Zero 
//' 
//' Turns a numerical matrix into a logical matrix, by checking for
//' every item whether it is (strictly) larger than zero. 
//' 
//' @param raw The numerical matrix 
//' 
//' @return Logical matrix 
// [[Rcpp::export]] 
Rcpp::LogicalMatrix threshold (arma::mat raw) {
  
  Rcpp:LogicalMatrix res(raw.n_rows, raw.n_cols) ; 
  
  // go over the matrix
  for (int i = 0; i < raw.n_rows; i ++) {
    for (int j = 0; j < raw.n_cols; j ++) {
      if (raw(i,j) > 0) {
        res(i, j) = true ; 
      } else {
        res(i, j) = false ; 
      }
    }
  }
  
  return res ; 
}

//' Generating Reports
//' 
//' Generates a number of reports. First, we sample from a multivariate normal distribution.
//' We then turn the sampled data into binary data, by checking whether the values are 
//' larger than zero or not.
//' \cr\cr
//' A report normally contains a least one drug and at least one event. In order to make 
//' sure that is the case, we check whether the generated report is indeed valid 
//' (see the function \code{\link{validReport}}). One can turn this off by setting
//' \code{valid_reports} to \code{FALSE}
//' 
//' @param n_reports The number of reports generated
//' @param means The mean of the multivariate normal distribution
//' @param L The Cholesky decomposition of the covariance matrix of the multivariate normal distribution
//' @param valid_reports When \code{TRUE}, only valid reports are added 
//' @param n_drugs The number of drugs 
//' @param n_events The number of events
//' @param verbose Verbosity. If \code{TRUE}, the function prints the 
//'                 number of reports that have been processed. Works only when valid_reports is \code{TRUE}
//' 
//' @return A logical matrix. Each row is a report. The number of columns is equal to \code{n_drugs + n_events}
//' @seealso \code{\link{validReport}}
// [[Rcpp::export]] 
Rcpp::LogicalMatrix generateReports (int n_reports, Rcpp::NumericVector means, arma::mat L, bool valid_reports, int n_drugs, int n_events, bool verbose) {
  
  if (!valid_reports) { // If we don't care whether the reports are 'valid'
    arma::mat raw = mvrnormArma(n_reports, means, L, L.n_cols) ; 
    return(threshold(raw)) ; 
    
  } else {
    
    // generate reports until we have enough valid ones
    int n_reports_generated = 0 ;
    
    Rcpp::LogicalMatrix reports(n_reports, means.size()) ; // will keep the final reports
    
    // variables to hold the current report 
    arma::mat report_raw;
    Rcpp::LogicalMatrix report ; 
    
    Progress p(n_reports, verbose) ;
    
    while (n_reports_generated < n_reports) {
      
      // generate one report
      report_raw = mvrnormArma(1, means, L, L.n_cols) ; 
      report = threshold(report_raw) ; 
      
      // check whether it is valid
      if (validReport(report, n_drugs, n_events)) {
        reports(n_reports_generated, _ ) = report(0, _ ); 
        n_reports_generated ++ ;
        p.increment() ; 
      }
    }
    
    return(reports) ; 
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector simulateReport(int n_drugs, int n_events,
                                         Rcpp::IntegerVector n_parents,
                                         Rcpp::NumericVector beta0, 
                                         Rcpp::NumericVector beta1, 
                                         Rcpp::IntegerVector parent_id, 
                                         bool verbose) {
  
  int i; 
  int n = n_drugs + n_events ; 
  double logit ; 
  
  // vector to hold the current report 
  Rcpp::IntegerVector report(n, -1) ; 
  
  // Keeps track of which indices in the report have already been set
  Rcpp::LogicalVector drawn(n, false) ;
  
  bool all_drawn ; // true only when the report is completely filled
  bool parent_drawn ; 
  
  do {
    // go over all nodes
    for (i = 0; i < n; i++) {
      if (!drawn[i]) { // check whether the node is already drawn
        
        logit = beta0[i] ; 
        
        parent_drawn = true ;
        if (n_parents[i] == 1) { 
          if (drawn[parent_id[i]-1]) {
            logit += beta1[i] * report[parent_id[i]-1] ; 
          } else { 
            parent_drawn = false ;  
          }
        }
        
        if (parent_drawn) { 
          report[i] = rbinom(1, 1, exp(logit) / (1 + exp(logit)))[0] ;  
          drawn[i] = true ; 
        } 
      }
    }
    
    // check whether all variates are drawn
    all_drawn = true ; 
    for (i = 0; i < n; i++) {
      if (!drawn[i]) { 
        all_drawn = false ;
        break ; 
      }
    }
  } while (!all_drawn) ; 
  
  return(report) ; 
  
}

//' Simulate a single report from a DAG
//' 
//' Returns a single report based on a DAG structure. 
//' 
//' @param n_drugs The number of drugs in the SR
//' @param n_events The number of events in the SR
//' @param id A list with the indices of the drugs and the events 
//' 
// [[Rcpp::export]] 
Rcpp::IntegerVector simulateReportDAG(int n_drugs, int n_events, 
                                      Rcpp::NumericVector beta0, 
                                      Rcpp::NumericMatrix betas, 
                                      bool verbose) {
  int i, j; 
  int n = n_drugs + n_events ; 
  double logit ; 

  // vector to hold the current report 
  Rcpp::IntegerVector report(n, -1) ; 
  
  // Keeps track of which indices in the report have already been set
  Rcpp::LogicalVector drawn(n, false) ;

  bool all_drawn ; // true only when the report is completely filled
  bool parents_drawn ; 
  
  do {
    // go over all nodes
    for (i = 0; i < n; i++) {
      if (!drawn[i]) { // check whether the node is already drawn
      
        logit = beta0[i] ; 
      
        parents_drawn = true ; 
        for (j = 0; j < n; j++) {
          if (fabs(betas(j, i)) > 0.00001) { // a parent 
            if (drawn[j]) {
                logit += betas(j, i) * report[j] ; 
            } else {
                parents_drawn = false ; 
                break ; 
            }
          }
        }
      
        if (parents_drawn) { 
          report[i] = rbinom(1, 1, exp(logit) / (1 + exp(logit)))[0] ;  
          drawn[i] = true ; 
        } 
      }
    }
    
    // check whether all variates are drawn
    all_drawn = true ; 
    for (i = 0; i < n; i++) {
      if (!drawn[i]) { 
        all_drawn = false ;
        break ; 
      }
    }
  } while (!all_drawn) ; 
 
  return(report) ; 
}


//' Create 2 x 2 Tables based on DAG generated SR data
//' 
//' Creates a data frame containing all 2 x 2 contingency tables 
//' from the results generated by generateSRDAG. 
//' See the R wrapper function \code{\link{create2x2TablesDAG}}
//' for more information. 
//' 
//' @param reports A binary matrix. Each row is a report
//' @param oddsratios The matrix with the odds ratio increase for each pair of variates
//' @param n_parents List with the number of parents for each node
//' @param prob_drugs A list with the marginal probabilities of the drugs 
//' @param prob_events A list with the marginal probabilities of the events
//' 
//' @return A dataframe. A description of the columns can be found in the commentary
//'         for the function \code{\link{create2x2TablesDAG}}
//' 
//' @seealso \code{\link{create2x2TablesDAG}}
// [[Rcpp::export]]
Rcpp::DataFrame create2x2TablesDAGRcpp (Rcpp::IntegerMatrix reports, 
                                        Rcpp::NumericVector prob_drugs, 
                                        Rcpp::NumericVector prob_events,
                                        Rcpp::IntegerVector n_parents,
                                        Rcpp::IntegerVector parent_id,
                                        Rcpp::NumericVector beta1) {
  int n_drugs = prob_drugs.size() ; 
  int n_events = prob_events.size() ; 
  int n_pairs = n_drugs * n_events ; 
  int k ; 
  
  // create vectors that will make up the data frame 
  Rcpp::IntegerVector drug_id (n_pairs) ;
  Rcpp::IntegerVector event_id (n_pairs) ;
  Rcpp::NumericVector prob_drug (n_pairs) ;
  Rcpp::NumericVector prob_event (n_pairs) ;
  Rcpp::NumericVector OR (n_pairs, 1.0) ;
  Rcpp::LogicalVector associated (n_pairs, false) ; 
  Rcpp::IntegerVector a (n_pairs, 0) ;
  Rcpp::IntegerVector b (n_pairs, 0) ;
  Rcpp::IntegerVector c (n_pairs, 0) ;
  Rcpp::IntegerVector d (n_pairs, 0) ;
  
  // run over all pairs and fill the vectors
  for (int i = 0; i < n_drugs; i ++) {
    for (int j = 0; j < n_events; j ++) {
      k = i*n_events + j ; // current pair index
      drug_id[k]        = i+1 ; 
      event_id[k]       = j+1 ; 
      prob_drug[k]      = prob_drugs[i] ; 
      prob_event[k]     = prob_events[j] ; 
      
      if (n_parents[j + n_drugs] == 1) {
        if (i == (parent_id[j + n_drugs]-1)) {
          OR[k] = exp(beta1[j + n_drugs]) ; 
          associated[k] = true ; 
        }
      }
    }
  }
  
  // compute the actual tables
  // go over all the reports
  bool drug, event ; 
  for (int r = 0; r < reports.nrow(); r ++) {
    // go over all drug-event pairs 
    for (int i = 0; i < n_drugs; i ++) {
      for (int j = 0; j < n_events; j ++) {
        k = i*n_events + j ; // pair index
        drug = (reports(r, i) == 1) ; 
        event = (reports(r, n_drugs + j) == 1) ; 
        if (drug) {
          if (event) {
            a[k] ++ ; 
          } else {
            c[k] ++ ; 
          }
        } else {
          if (event) {
            b[k] ++ ; 
          } else {
            d[k] ++ ; 
          }
        }
      } 
    }
  }
  
  return Rcpp::DataFrame::create( Named("drug_id") = drug_id, 
                                  Named("event_id") = event_id, 
                                  Named("prob_drug") = prob_drug, 
                                  Named("prob_event") = prob_event,
                                  Named("or") = OR,
                                  Named("associated") = associated,
                                  Named("a") = a, 
                                  Named("b") = b,
                                  Named("c") = c,
                                  Named("d") = d
  ) ; 
}


//' Create 2 x 2 Tables 
//' 
//' Creates a data frame containing all 2 x 2 contingency tables 
//' from the results generated by generateSRData. 
//' See the R wrapper function \code{\link{create2x2Tables}}
//' for more information. 
//' 
//' @param reports A binary matrix. Each row is a report
//' @param corrmat The correlation matrix
//' @param prob_drugs A list with the marginal probabilities of the drugs 
//' @param prob_events A list with the marginal probabilities of the events
//' 
//' @return A dataframe. A description of the columns can be found in the commentary
//'         for the function \code{\link{create2x2Tables}}
//' 
//' @seealso \code{\link{create2x2Tables}}
// [[Rcpp::export]] 
Rcpp::DataFrame create2x2TablesRcpp (Rcpp::IntegerMatrix reports, Rcpp::NumericMatrix corrmat, Rcpp::NumericVector prob_drugs, Rcpp::NumericVector prob_events) {
  
  int n_drugs = prob_drugs.size() ; 
  int n_events = prob_events.size() ; 
  int n_pairs = n_drugs * n_events ; 
  int k ; 
  
  // create vectors that will make up the data frame 
  Rcpp::IntegerVector drug_id (n_pairs) ;
  Rcpp::IntegerVector event_id (n_pairs) ;
  Rcpp::NumericVector prob_drug (n_pairs) ;
  Rcpp::NumericVector prob_event (n_pairs) ;
  Rcpp::IntegerVector n_corr_drugs (n_pairs, 0) ;
  Rcpp::IntegerVector n_corr_events (n_pairs, 0) ;
  Rcpp::NumericVector max_corr_drug (n_pairs, 0.0) ;
  Rcpp::NumericVector max_corr_event (n_pairs, 0.0) ;  
  Rcpp::NumericVector corr (n_pairs) ;
  Rcpp::NumericVector cov (n_pairs, 0.0) ;
  Rcpp::NumericVector commonprob (n_pairs) ;
  Rcpp::LogicalVector associated (n_pairs, false) ; 
  Rcpp::IntegerVector a (n_pairs, 0) ;
  Rcpp::IntegerVector b (n_pairs, 0) ;
  Rcpp::IntegerVector c (n_pairs, 0) ;
  Rcpp::IntegerVector d (n_pairs, 0) ;
  
  // get for every drug how often they correlate with other drugs and the maximum 
  // correlation they show. 
  // do the same for the events
  Rcpp::IntegerVector n_corr_drugs_temp (n_drugs, 0) ;
  Rcpp::IntegerVector n_corr_events_temp (n_events, 0) ;
  Rcpp::NumericVector max_corr_drug_temp (n_drugs, 0.0) ;
  Rcpp::NumericVector max_corr_event_temp (n_events, 0.0) ;
  
  // walk over all unique pairs of drugs 
  for (int i = 0 ; i < n_drugs; i ++) {
    for (int j = i + 1; j < n_drugs; j ++) {
      if (corrmat(i,j) != 0.0) {
        n_corr_drugs_temp[i] ++ ; 
        n_corr_drugs_temp[j] ++ ; 
        if (corrmat(i,j) > max_corr_drug_temp[i]) {
          max_corr_drug_temp[i] = corrmat(i,j) ; 
        } 
        if (corrmat(i,j) > max_corr_drug_temp[j]) {
          max_corr_drug_temp[j] = corrmat(i,j) ; 
        } 
      }
    }
  }
  
  // do the same for the events
  double cor  ; 
  for (int i = 0 ; i < n_events; i ++) {
    for (int j = i + 1; j < n_events; j ++) {
      cor = corrmat(i + n_drugs, j + n_drugs) ; 
      if (cor != 0.0) {
        n_corr_events_temp[i] ++ ; 
        n_corr_events_temp[j] ++ ; 
        if (cor > max_corr_event_temp[i]) {
          max_corr_event_temp[i] = cor ; 
        } 
        if (cor > max_corr_event_temp[j]) {
          max_corr_event_temp[j] = cor ; 
        } 
      }
    }
  }
  
  
  // run over all pairs and fill the vectors
  for (int i = 0; i < n_drugs; i ++) {
    for (int j = 0; j < n_events; j ++) {
      k = i*n_events + j ; // current pair index
      drug_id[k]        = i+1 ; 
      event_id[k]       = j+1 ; 
      prob_drug[k]      = prob_drugs[i] ; 
      prob_event[k]     = prob_events[j] ; 
      n_corr_drugs[k]   = n_corr_drugs_temp[i] ; 
      n_corr_events[k]  = n_corr_events_temp[j] ; 
      max_corr_drug[k]  = max_corr_drug_temp[i] ; 
      max_corr_drug[k]  = max_corr_event_temp[j]; 
      corr[k]           = corrmat(i, j + n_drugs) ; 
      cov[k]            = corr[k] * sqrt(prob_drug[k] * (1 - prob_drug[k]) * prob_event[k] * (1 - prob_event[k])) ; 
      commonprob[k]     = cov[k] + prob_drug[k]  * prob_event[k] ; 
      associated[k]     = (corr[k] != 0.0) ; 
    }
  }
  
  // compute the actual tables
  // go over all the reports
  bool drug, event ; 
  for (int r = 0; r < reports.nrow(); r ++) {
    // go over all drug-event pairs 
    for (int i = 0; i < n_drugs; i ++) {
      for (int j = 0; j < n_events; j ++) {
        k = i*n_events + j ; // pair index
        drug = (reports(r, i) == 1) ; 
        event = (reports(r, n_drugs + j) == 1) ; 
        if (drug) {
          if (event) {
            a[k] ++ ; 
          } else {
            c[k] ++ ; 
          }
        } else {
          if (event) {
            b[k] ++ ; 
          } else {
            d[k] ++ ; 
          }
        }
      } 
    }
  }
  
  return Rcpp::DataFrame::create( Named("drug_id") = drug_id, 
                                  Named("event_id") = event_id, 
                                  Named("prob_drug") = prob_drug, 
                                  Named("prob_event") = prob_event, 
                                  Named("n_corr_drugs") = n_corr_drugs, 
                                  Named("n_corr_events") = n_corr_events, 
                                  Named("max_corr_drug") = max_corr_drug,
                                  Named("max_corr_event") = max_corr_event, 
                                  Named("corr") = corr, 
                                  Named("cov") = cov,
                                  Named("commonprob") = commonprob, 
                                  Named("associated") = associated,
                                  Named("a") = a, 
                                  Named("b") = b,
                                  Named("c") = c,
                                  Named("d") = d
  ) ; 
}

