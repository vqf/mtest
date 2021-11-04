#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List traverse(List l) {
  NumericVector desc = as<NumericVector>(l["desc"]);
  NumericVector sdesc = l["sdesc"];
  uint32_t h = desc(0);
  NumericVector t = tail(desc, -1);
  uint32_t sh = sdesc(0);
  NumericVector st = tail(sdesc, -1);
  List result = List::create();
  result.push_front(l);
  if (t.length() == 0){
    return result;
  }
  List nxt = List::create();
  nxt["desc"] = t;
  nxt["sdesc"] = st;
  while(h > 0){
    NumericVector ndesc = as<NumericVector>(nxt["desc"]);
    NumericVector nsdesc = nxt["sdesc"];
    double v = (double) l["val"];
    double si = (double) h;
    double fi = (double) ndesc(0);
    double ssi = (double) sh;
    double sfi = (double) nsdesc(0);
    nxt["val"] = (double) (v * si * (1 + sfi) / ((1 + fi) * ssi));
    NumericVector d0 = nxt["desc"];
    d0(0) = d0(0) + 1;
    NumericVector sd0 = nxt["sdesc"];
    sd0(0) = sd0(0) + 1;
    sh--;
    h--;
    List n = traverse(nxt);
    uint32_t i = 0;
    for (i = 0; i < n.length(); i++){
      List o = n(i);
      NumericVector dsc = o["desc"];
      NumericVector sdsc = o["sdesc"];
      dsc.push_front(h);
      sdsc.push_front(sh);
      List r = List::create(_["val"] = o["val"],
                            _["desc"] = dsc,
                            _["sdesc"] = sdsc);
      result.push_back(r);
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix sl(NumericMatrix from) {
  NumericVector h = head(from, 1);
  int32_t n = from.ncol();
  NumericMatrix t = from(Range(1, n), _);
  if (t.length() == 0){

  }
  return t;
}

double tmtest2(NumericMatrix l, double v, double cutoff, NumericVector nss, NumericVector sus, double vo){
  double result = 0;
  NumericMatrix st = Rcpp::clone(l);
  NumericVector ns = Rcpp::clone(nss);
  NumericVector su = Rcpp::clone(sus);
  uint32_t nd2 = st(1, 0);
  vo = vo * ((ns(1) + 1) * (st(0, 0) + 1) * (su(1) + 1)) / ((ns(0) + 1) * (st(1, 1) + 1) * (su(0) + 1));
  while (st(1, 1) < nd2){;
    v = v + vo / (ns(1) + 1);
    if (v <= cutoff){
      result += v;
    }
    else{
      return result;
    }
    st(1, 0) = st(1, 0) - 1;
    st(1, 1) = st(1, 1) + 1;
    su(0) = su(0) - 1;
    su(1) = su(1) + 1;
    vo = vo * (st(1, 0)+1) * (1 + su(1)) / ((1 + st(1, 1)) * (su(0)+1));
  }
  return result;
}


// [[Rcpp::export]]
double tmtest(NumericMatrix l, double v, double cutoff, NumericVector nss, NumericVector sus, double vo){
  double result = 0;
  if (v <= cutoff){
    result += v;
  }
  NumericMatrix st = Rcpp::clone(l);
  NumericVector ns = Rcpp::clone(nss);
  NumericVector su = Rcpp::clone(sus);
  result += tmtest2(st, v, cutoff, ns, su, vo);
  uint32_t nd1 = (uint32_t) ns(0);
  while (st(0, 0) < nd1){
    v = v + vo / ((double)(ns(0)) + 1);
    if (v <= cutoff){
      result += v;
    }
    else{
      return result;
    }
    st(0, 0) = st(0, 0) + 1;
    st(0, 1) = st(0, 1) - 1;
    su(0) = su(0) + 1;
    su(1) = su(1) - 1;
    vo = vo * ((double)(st(0, 1))+1) * (1 + (double)(su(0))) / ((1 + (double)(st(0, 0))) * ((double)(su(1))+1));
    result += tmtest2(st, v, cutoff, ns, su, vo);
  }
  return result;
}


// [[Rcpp::export]]
List mtest(List l, double cutoff, uint8_t rw = 1, uint8_t cl = 1){
  double v = l["val"];
  NumericMatrix dm = l["desc"];
  NumericVector ds = l["sdesc"];
  NumericMatrix m = Rcpp::clone(dm);
  NumericVector s = Rcpp::clone(ds);
  uint8_t nc = m.ncol();
  uint8_t nr = m.nrow();
  if (m(rw, cl) < 0 || s(cl) < 0){
    stop("Value lower than zero");
  }
  m(rw, cl + 1) = m(rw, cl + 1) + 1;
  s(cl + 1) = s(cl + 1) + 1;
  double newv = v * m(rw, cl) * s(cl + 1) / (m(rw, cl + 1) * s(cl));
  m(rw, cl) = m(rw, cl) - 1;
  s(cl) = s(cl) - 1;
  //Rprintf("%f\n", newv);
  //Rf_PrintValue(m);
  List result = List::create( _["val"] = newv,
                              _["desc"] = m,
                              _["sdesc"] = s,
                              _["r"] = 0
  );
  double r = 0;
  if (newv <= cutoff){
    r = newv;
  }
  if (rw < (nr-1)){
    uint8_t i = 0;
    for (i = (rw + 1); i < nr; i++){
      //Rprintf("a: %u\n", i);
      //Rf_PrintValue(m);
      List tr = mtest(List::create(_["val"] = newv,
                                    _["desc"] = m,
                                    _["sdesc"] = s,
                                    _["r"] = r),
                                    cutoff,
                                    i, 0);
      double newr = tr["r"];
      r += newr;
    }
  }
  if (m(rw, cl) > 0){
    //Rprintf("b: %u\n", m(rw, cl));
    //Rf_PrintValue(m);
    List tr = mtest(List::create(_["val"] = newv,
                                  _["desc"] = m,
                                  _["sdesc"] = s,
                                  _["r"] = r),
                                  cutoff,
                                  rw, cl);
    double newr = tr["r"];
    r += newr;
  }
  if (cl < (uint8_t)(nc-2)){
    //Rprintf("c: %u\n", cl);
    //Rf_PrintValue(m);
    List tr = mtest(List::create(_["val"] = newv,
                                  _["desc"] = m,
                                  _["sdesc"] = s,
                                  _["r"] = r),
                                  cutoff,
                                  rw, cl + 1);
    double newr = tr["r"];
    r += newr;
  }
  result["r"] = r;
  return(result);
}

// [[Rcpp::export]]
double sumlog(double n1, double n2){
  if (n1 == n2){
    return n1 + log(2);
  }
  double l = (n1 > n2) ? n1 : n2;
  double s = (n1 > n2) ? n2 : n1;
  double result = l;
  double Q = s - l;
  double q = exp(Q);
  double qq = q * q;
  double c = 1 - q;
  uint32_t i = 1;
  double k = q;
  double next = q * (1-q/2);
  while (abs(next) > 1e-20){
    i = i + 2;
    result += next;
    k *= qq;
    next = k*((1 + i * c)/(i*(i + 1)));
  }
  return result;
}

// [[Rcpp::export]]
List lmtest(List l, double cutoff, uint8_t rw = 1, uint8_t cl = 1){
  double v = l["val"];
  NumericMatrix dm = l["desc"];
  NumericVector ds = l["sdesc"];
  NumericMatrix m = Rcpp::clone(dm);
  NumericVector s = Rcpp::clone(ds);
  uint8_t nc = m.ncol();
  uint8_t nr = m.nrow();
  if (m(rw, cl) < 0 || s(cl) < 0){
    stop("Value lower than zero");
  }
  m(rw, cl + 1) = m(rw, cl + 1) + 1;
  s(cl + 1) = s(cl + 1) + 1;
  double newv = v + log(m(rw, cl)) + log(s(cl + 1)) - log(m(rw, cl + 1)) - log(s(cl));
  m(rw, cl) = m(rw, cl) - 1;
  s(cl) = s(cl) - 1;
  //Rprintf("%f\t%f\n", v, newv);
  //stop("");
  //Rf_PrintValue(m);
  List result = List::create( _["val"] = newv,
                              _["desc"] = m,
                              _["sdesc"] = s,
                              _["r"] = 0
  );
  double r = 0;
  if (newv <= cutoff){
    r = newv;
  }
  if (rw < (nr-1)){
    uint8_t i = 0;
    for (i = (rw + 1); i < nr; i++){
      //Rprintf("a: %u\n", i);
      //Rf_PrintValue(m);
      List tr = lmtest(List::create(_["val"] = newv,
                                   _["desc"] = m,
                                   _["sdesc"] = s,
                                   _["r"] = r),
                                   cutoff,
                                   i, 0);
      double newr = tr["r"];
      if (r == 0){
        r = newr;
      }
      else{
        if (newr != 0){
          r = sumlog(r, newr);
        }
      }
    }
  }
  if (m(rw, cl) > 0){
    //Rprintf("b: %u\n", m(rw, cl));
    //Rf_PrintValue(m);
    List tr = lmtest(List::create(_["val"] = newv,
                                 _["desc"] = m,
                                 _["sdesc"] = s,
                                 _["r"] = r),
                                 cutoff,
                                 rw, cl);
    double newr = tr["r"];
    if (r == 0){
      r = newr;
    }
    else{
      if (newr != 0){
        r = sumlog(r, newr);
      }
    }
  }
  if (cl < (uint8_t)(nc-2)){
    //Rprintf("c: %u\n", cl);
    //Rf_PrintValue(m);
    List tr = lmtest(List::create(_["val"] = newv,
                                 _["desc"] = m,
                                 _["sdesc"] = s,
                                 _["r"] = r),
                                 cutoff,
                                 rw, cl + 1);
    double newr = tr["r"];
    if (r == 0){
      r = newr;
    }
    else{
      if (newr != 0){
        r = sumlog(r, newr);
      }
    }
  }
  result["r"] = r;
  return(result);
}

