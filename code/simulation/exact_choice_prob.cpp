// Exact choice probabilities for the loss-averse Bayesian attribute model, by
// convolution on a utility grid. Called from piece_choice_logprob() in
// code/simulation/simulate_bayesian_attributes.R, which explains where the
// inputs come from and does the extrapolation over grid levels.
//
// One row of the input matrices is one trial at one candidate sigma; the four
// columns are the trial's boxes:
//   eta     mean of the box's estimate relative to the reference point, w (a - m0)
//   kappa   its SD, w sigma / sqrt(t); 0 for a box that was never seen, which
//           then sits exactly at the reference point and is left out
//   sgn     +1 for a box of the left option, -1 for the right option
//   lambda  loss aversion that applies to the box
// A box's utility relative to the reference point is Y = X^+ + lambda X^- with
// X ~ Normal(eta, kappa^2): a bell curve whose lower half is stretched by
// lambda. An option's utility is the sum of its boxes' Y; sgn only says which
// option a box belongs to.
//
// For each row and each grid level l = 0, ..., levels - 1, with spacing
// h = factor * (smallest kappa of the row) / 2^l:
//   1. each box's probability of landing in every grid cell is its density at
//      the cell midpoint (k + 1/2) h times h, for the cells covering +/- nsd
//      SDs of the box. The reference point (0) is a cell edge, so no cell
//      straddles the kink.
//   2. an option's utility is the sum of its boxes, so its cell probabilities
//      are the convolution of its boxes' (a direct double loop; the vectors
//      are short).
//   3. Pr(left > right) adds, over the right option's cells, that cell's
//      probability times the probability that the left option is in a higher
//      cell, plus half the probability that it is in the same cell.
//      Pr(right > left) is computed the same way (not as 1 minus the other),
//      so that small probabilities keep their relative accuracy.
// Returns an n x (2 levels) matrix: columns 1..levels hold Pr(left) at each
// level, columns levels+1..2 levels hold Pr(right).
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

static const double INV_SQRT_2PI = 0.3989422804014327;

// The cells of one box's Y: writes the cell probabilities to m and returns
// the index of the first cell (cell k covers [k h, (k + 1) h)). Below 0 the
// density is that of X stretched by lambda; above 0 it is X's own.
static long box_cells(double eta, double kappa, double lambda, double h, double nsd, std::vector<double>& m) {
  // Support of Y: the image of eta -/+ nsd kappa under x -> (x if x >= 0, else lambda x), which is increasing.
  double x1 = eta - nsd * kappa, x2 = eta + nsd * kappa;
  double ymin = (x1 >= 0 ? x1 : lambda * x1);
  double ymax = (x2 >= 0 ? x2 : lambda * x2);
  long klo = (long) std::floor(ymin / h);
  long L   = (long) std::ceil(ymax / h) - klo;
  if (L < 1) L = 1;
  m.resize(L);
  for (long j = 0; j < L; ++j) {
    double y   = (klo + j + 0.5) * h;
    double lam = (y < 0) ? lambda : 1.0;
    double u   = (y / lam - eta) / kappa;
    m[j] = h * INV_SQRT_2PI * std::exp(-0.5 * u * u) / (lam * kappa);
  }
  return klo;
}

// An option's cell probabilities and the position of its first cell, in
// half-units of h (a box's cell j sits at (klo + j + 1/2) h, i.e. 2 (klo + j) + 1
// half-units; cells step by 2 half-units; a sum of boxes adds their positions).
struct Side { std::vector<double> m; long pos0; };

static void side_sum(const std::vector<std::vector<double>*>& boxes, const std::vector<long>& klo,
                     Side& s, std::vector<double>& tmp) {
  if (boxes.empty()) { s.m.assign(1, 1.0); s.pos0 = 0; return; }   // no box seen: exactly at the reference point
  s.m = *boxes[0];
  s.pos0 = 2 * klo[0] + 1;
  for (size_t b = 1; b < boxes.size(); ++b) {
    const std::vector<double>& a = s.m;
    const std::vector<double>& c = *boxes[b];
    tmp.assign(a.size() + c.size() - 1, 0.0);
    for (size_t i = 0; i < a.size(); ++i) {
      double ai = a[i];
      if (ai == 0) continue;
      for (size_t j = 0; j < c.size(); ++j) tmp[i + j] += ai * c[j];
    }
    s.m.swap(tmp);
    s.pos0 += 2 * klo[b] + 1;
  }
}

// Pr(L > R) and Pr(R > L), a tie counting half for each. The running sums of
// the left cells are taken upward and downward separately, never as "total
// minus the rest", so a small tail is not lost to rounding.
static void compare(const Side& L, const Side& R, std::vector<double>& pre, std::vector<double>& suf,
                    double& pL, double& pR) {
  long nL = (long) L.m.size();
  pre.assign(nL + 1, 0.0);
  suf.assign(nL + 1, 0.0);
  for (long l = 0; l < nL; ++l) pre[l + 1] = pre[l] + L.m[l];     // pre[l]: left cells before l
  for (long l = nL - 1; l >= 0; --l) suf[l] = suf[l + 1] + L.m[l]; // suf[l]: left cells from l on
  pL = 0.0; pR = 0.0;
  for (long r = 0; r < (long) R.m.size(); ++r) {
    double mr = R.m[r];
    if (mr == 0) continue;
    long d = (R.pos0 + 2 * r) - L.pos0;      // right cell's position relative to the first left cell
    long above, below;
    double tie = 0.0;
    if (d < 0) { above = 0; below = 0; }     // every left cell is higher
    else {
      above = d / 2 + 1;                     // first left cell higher than the right cell
      below = (d + 1) / 2;                   // number of left cells lower than it
      if (d % 2 == 0 && d / 2 < nL) tie = 0.5 * L.m[d / 2];
    }
    if (above > nL) above = nL;
    if (below > nL) below = nL;
    pL += mr * (suf[above] + tie);
    pR += mr * (pre[below] + tie);
  }
}

// [[Rcpp::export]]
NumericMatrix exact_choice_prob(NumericMatrix eta, NumericMatrix kappa, NumericMatrix sgn, NumericMatrix lambda,
                                double factor, int levels, double nsd) {
  int n = eta.nrow();
  NumericMatrix out(n, 2 * levels);
  std::vector<double> m[4], tmp, pre, suf;
  Side L, R;
  for (int i = 0; i < n; ++i) {
    double kmin = R_PosInf;
    int nact = 0;
    for (int p = 0; p < 4; ++p) if (kappa(i, p) > 0) { ++nact; kmin = std::min(kmin, kappa(i, p)); }
    if (nact == 0) {                                   // nothing was seen: the options tie
      for (int lev = 0; lev < levels; ++lev) { out(i, lev) = 0.5; out(i, levels + lev) = 0.5; }
      continue;
    }
    for (int lev = 0; lev < levels; ++lev) {
      double h = factor * kmin / std::pow(2.0, lev);
      std::vector<std::vector<double>*> lb, rb;
      std::vector<long> lk, rk;
      for (int p = 0; p < 4; ++p) {
        if (kappa(i, p) <= 0) continue;
        long klo = box_cells(eta(i, p), kappa(i, p), lambda(i, p), h, nsd, m[p]);
        if (sgn(i, p) > 0) { lb.push_back(&m[p]); lk.push_back(klo); }
        else               { rb.push_back(&m[p]); rk.push_back(klo); }
      }
      side_sum(lb, lk, L, tmp);
      side_sum(rb, rk, R, tmp);
      double pL, pR;
      compare(L, R, pre, suf, pL, pR);
      out(i, lev) = pL;
      out(i, levels + lev) = pR;
    }
  }
  return out;
}
