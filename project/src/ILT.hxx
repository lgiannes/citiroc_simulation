/*----------------------------------------------------------------
  Laplace Inversion Source Code
  Copyright © 2010 James R. Craig, University of Waterloo
----------------------------------------------------------------*/
#include <complex>
#include <functional>
#include <iomanip>
#include <vector>

#include <TGraph.h>
#include <TCanvas.h>
#include <TLine.h>

using namespace std;

namespace ILT {

    typedef complex<double> cmplex;

    const double PI = 3.141592653589793238462643;      // pi
    const int MAX_LAPORDER = 60;

/***********************************************************************
        LaplaceInversion
************************************************************************
Returns Inverse Laplace transform f(t) of function F(s),
where s is complex, evaluated at time t

  f(t) = 1/2*PI*i \intfrmto{gamma-i\infty}{gamma+i\infty} exp(st)*F(s) ds

Based upon De Hoog et al., 1982, An improved method for numerical inversion of
Laplace transforms, SIAM J. Sci. Stat. Comput.

Speed of algorithm is primarily a function of M, the order of the
Taylor series expansion.

Inputs: F(s), a function pointer to some complex function of a single
        complex variable
        t is the desired time of evaluation of f(t)
        tolerance is the required accuracy of f(t)
-----------------------------------------------------------------------*/

//    double LaplaceInversion(cmplex       (*F)(const cmplex &s),
//                            const double &t,
//                            const double tolerance) {

    template<typename Func>
    inline double LaplaceInversion(Func F, const double &t, const double tolerance) {
      ILT::cmplex s;
      ILT::cmplex result = F(s);

//--Variable declaration---------------------------------
      int i, n, m, r;          //counters & intermediate array indices
      int M(40);            //order of Taylor Expansion (must be less than MAX_LAPORDER)
      double DeHoogFactor(4.0);//DeHoog time factor
      double T;                //Period of DeHoog Inversion formula
      double gamma;            //Integration limit parameter
      cmplex h2M, R2M, z, dz;   //Temporary variables

      static cmplex Fctrl[2 * MAX_LAPORDER + 1];
      static cmplex e[2 * MAX_LAPORDER][MAX_LAPORDER];
      static cmplex q[2 * MAX_LAPORDER][MAX_LAPORDER];
      static cmplex d[2 * MAX_LAPORDER + 1];
      static cmplex A[2 * MAX_LAPORDER + 2];
      static cmplex B[2 * MAX_LAPORDER + 2];

//Calculate period and integration limits------------------------------------
      T = DeHoogFactor * t;
      gamma = -0.5 * log(tolerance) / T;

//Calculate F(s) at evalution points gamma+IM*i*PI/T for i=0 to 2*M-1--------
//This is likely the most time consuming portion of the DeHoog algorithm
      Fctrl[0] = 0.5 * F(gamma);
      for (i = 1; i <= 2 * M; i++) {
        s = cmplex(gamma, i * PI / T);
        Fctrl[i] = F(s);
      }

//Evaluate e and q ----------------------------------------------------------
//eqn 20 of De Hoog et al 1982
      for (i = 0; i < 2 * M; i++) {
        e[i][0] = 0.0;
        q[i][1] = Fctrl[i + 1] / Fctrl[i];
      }
      e[2 * M][0] = 0.0;

      for (r = 1; r <= M - 1; r++) //one minor correction - does not work for r<=M, as suggested in paper
      {
        for (i = 2 * (M - r); i >= 0; i--) {
          if ((i < 2 * (M - r)) && (r > 1)) {
            q[i][r] = q[i + 1][r - 1] * e[i + 1][r - 1] / e[i][r - 1];
          }
          e[i][r] = q[i + 1][r] - q[i][r] + e[i + 1][r - 1];
        }
      }

//Populate d vector-----------------------------------------------------------
      d[0] = Fctrl[0];
      for (m = 1; m <= M; m++) {
        d[2 * m - 1] = -q[0][m];
        d[2 * m] = -e[0][m];
      }

//Evaluate A, B---------------------------------------------------------------
//Eqn. 21 in De Hoog et al.
      z = cmplex(cos(PI * t / T), sin(PI * t / T));

      A[0] = 0.0;
      B[0] = 1.0; //A_{-1},B_{-1} in De Hoog
      A[1] = d[0];
      B[1] = 1.0;
      for (n = 2; n <= 2 * M + 1; n++) {
        dz = d[n - 1] * z;
        A[n] = A[n - 1] + dz * A[n - 2];
        B[n] = B[n - 1] + dz * B[n - 2];
      }

//Eqn. 23 in De Hoog et al.
      h2M = 0.5 * (1.0 + z * (d[2 * M - 1] - d[2 * M]));
      R2M = -h2M * (1.0 - sqrt(1.0 + (z * d[2 * M] / h2M / h2M)));

//Eqn. 24 in De Hoog et al.
      A[2 * M + 1] = A[2 * M] + R2M * A[2 * M - 1];
      B[2 * M + 1] = B[2 * M] + R2M * B[2 * M - 1];

//Final result: A[2*M]/B[2*M]=sum [F(gamma+itheta)*exp(itheta)]-------------
      return 1.0 / T * exp(gamma * t) * (A[2 * M + 1] / B[2 * M + 1]).real();
    }
//-----------------------------------------------------------------------
    inline double inverseLaplace(function<cmplex(cmplex)> F, double t, int N=32) {
      if (t <= 0) {
        throw std::invalid_argument("Time t must be positive.");
      }
//      if (t == 0) return 0.0;  // Prevent division by zero issues

      cmplex I(0, 1);  // Imaginary unit
      double pi = M_PI;
      cmplex sum = 0.0;

      // Constants
      double shift = 0.1;
      // Shift contour to the right in case there is a pole on the positive real axis : Note the contour will
      // not be optimal since it was originally devoloped for function with
      // singularities on the negative real axis
      // For example take F(s) = 1/(s-1), it has a pole at s = 1, the contour needs to be shifted with one
      // unit, i.e shift  = 1.
      double sigma = 0.6122;
      double mu = 0.5017;
      double nu = 0.2645;
      double alpha = 0.6407;
      double h = 2 * pi / (double)N;

      std::vector<double> Res, Ims;
      TGraph* contour = new TGraph();

      for (int k = 1; k <= N; k++) {
        double theta_k = - pi + (k - 0.5) * h;
        cmplex s = (double)N/t*( -sigma + theta_k * (mu / tan(alpha*theta_k) + I * nu) ) + shift;
        Res.push_back(s.real());
        Ims.push_back(s.imag());
        contour->SetPoint(k-1, s.real(), s.imag());
        cmplex s_prime = (double)N/t*( mu / tan(alpha*theta_k) + mu * alpha * theta_k/cos(alpha*theta_k)/cos(alpha*theta_k) + I * nu );
        cmplex term = exp(s*t) * F(s) * s_prime;
        std::cout<<k<<"  cos(Im(s*t)): "<<cos(s.imag()*t)<<"  s': "<<s_prime<<" e^(st): "<<exp(s*t)<<" F(s): "<<F(s)<<"  term: "<<term<<std::endl;
        sum += (1.0/(I*(double)N)) * term;
      }

      std::cout<<"sum:     "<<sum<<std::endl;
//      std::cout<<"Re(sum): "<<sum.real()<<std::endl;

      TCanvas* c = new TCanvas("c", "c", 800, 800);
      contour->Draw("APL*");
      TLine* im_axis = new TLine(0, -1000, 0, 1000);
      im_axis->Draw();
      TLine* re_axis = new TLine(-1000, 0, 1000, 0);
      re_axis->Draw();
      c->SaveAs(Form("contour.png"));
      return (sum).real();
    }

    inline double deHoogInverseLaplace(const std::function<cmplex(cmplex)>& F, double t, double tol = 1e-5, int max_iter = 100) {
      const double gamma = 0.0; // Shift parameter
      const double T = 2 * t;   // Scaling parameter
      const cmplex I(0.0, 1.0); // Imaginary unit

      std::vector<cmplex> p(2 * max_iter + 1);
      std::vector<cmplex> f(2 * max_iter + 1);

      // Compute p_k and f_k
      for (int k = 0; k <= 2 * max_iter; ++k) {
        double theta_k = M_PI * k / max_iter;
        p[k] = gamma + (I * theta_k / T);
        f[k] = F(p[k]) * std::exp(p[k] * t);
      }

      // Apply the ε-algorithm for convergence acceleration
      std::vector<std::vector<cmplex>> epsilon(2 * max_iter + 1, std::vector<cmplex>(2 * max_iter + 1, 0.0));

      for (int i = 0; i <= 2 * max_iter; ++i) {
        epsilon[i][0] = f[i];
      }

      for (int j = 1; j <= 2 * max_iter; ++j) {
        for (int i = 0; i <= 2 * max_iter - j; ++i) {
          epsilon[i][j] = epsilon[i + 1][j - 1] + 1.0 / (epsilon[i + 1][j - 1] - epsilon[i][j - 1]);
        }
      }

      cmplex result = epsilon[0][2 * max_iter];

      return result.real();
    }



};

//double testf(const double &t)
//{
//double a=0.333,b=4.0;
//
//switch(test_case){
//case(0):{return a/pow(4*PI*pow(t,3.0),0.5)*exp(-a*a/4.0/t); break;}
//case(1):{return sin(a*t+b);                                 break;}
//case(2):{return exp(-a*t);                                  break;}
//case(3):{return t;                                          break;}
//default:{return 1;                                          break;}
//}
//}
////-----------------------------------------------------------------------
//cmplex testF(const cmplex &s)
//{
//double a=0.333,b=4.0;
//switch(test_case){
//case(0):{return exp(-a*pow(s,0.5));                         break;}
//case(1):{return (sin(b)*s+a*cos(b))/(s*s+a*a);              break;}
//case(2):{return 1.0/(s+a);                                  break;}
//case(3):{return 1.0/s/s;                                    break;}
//default:{return 1.0/s;                                      break;}
//}
//}
////-----------------------------------------------------------------------
//void test()
//{
//  double F,f,max_err(0),t;
//
//  cout <<"Testing inverse Laplace transform..."<<endl;
//  ofstream LAPTEST;
//  LAPTEST.open("LaplaceInversionTest.csv");
//  LAPTEST<<"t,f(t) analytic,f(t) numerical, error"<<endl;
//  LAPTEST.precision(12);
//  for (t=0.0; t<12; t+=0.1)
//  {
//    F=LaplaceInversion(testF,t,1e-8);
//    f=testf(t);
//    upperswap(max_err,fabs(F-f));
//    LAPTEST<<t<<","<<f<<","<<F<<","<<F-f<<endl;
//  }
//  LAPTEST.close();
//  cout <<"...Inverse Laplace transform evaluated. Max Error= "<<max_err<<endl;
//}
////-----------------------------------------------------------------------