// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY:  Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97

#ifndef TYPEMETRIX
#define TYPEMETRIX MetricAnIso
#endif

namespace bamg {

typedef P2<double,double> D2;
typedef P2xP2<double,double> D2xD2;

class MetricAnIso;
class MatVVP2x2;
class MetricIso;

typedef TYPEMETRIX Metric;


class MetricIso {

   friend class MatVVP2x2;
   double h;

 public:
   MetricIso(double a): h(a) { }

   MetricIso(const MetricAnIso M);
   MetricIso(double a11,
             double a21,
             double a22);
   MetricIso(): h(1) { }

   MetricIso(const double    a[3],
             const MetricIso m0,
             const MetricIso m1,
             const MetricIso m2)
           : h(hinterpole ? (a[0]*m0.h+a[1]*m1.h+a[2]*m2.h)
	          : 1/sqrt(a[0]/(m0.h*m0.h)+a[1]/(m1.h*m1.h)+a[2]/(m2.h*m2.h))) { }

   MetricIso(const double    a,
             const MetricIso ma,
             const double    b,
             const MetricIso mb)
           : h(hinterpole ? a*ma.h+b*mb.h : 1/sqrt(a/(ma.h*ma.h)+b/(mb.h*mb.h))) { }

   R2 Orthogonal(const R2 A) const { return R2(-h*A.y,h*A.x); }
   R2 Orthogonal(const I2 A) const { return R2(-h*A.y,h*A.x); }
   double operator()(R2 x) const { return sqrt((x,x))/h; }
   double operator()(R2 x,
                     R2 y) const
   { return ((x,y))/(h*h); }
   int IntersectWith(MetricIso M)
   {
       int r = 0;
       if (M.h<h)
          r = 1, h = M.h;
       return r;
   }
   MetricIso operator*(double c) const { return MetricIso(double(h)/c); } 
   MetricIso operator/(double c) const { return MetricIso(double(h)*c); }

   void Box(double& hx,
            double& hy)
        { hx = hy = h; }

   friend ostream& operator<<(      ostream&   f,
	                          const MetricIso& M)
   {
      f << " h = " << M.h << ";";
      return f;
   }
};


class MetricAnIso {

 public:
  friend class MatVVP2x2;
  double a11, a21, a22;
  MetricAnIso(double a) { a11 = 1/(a*a); a22 = 1/(a*a); a21 = 0; }

  MetricAnIso(double a,
              double b,
              double c)
  { a11 = a; a21 = b; a22 = c; }

  MetricAnIso() { }

  MetricAnIso(const double       a[3],
              const MetricAnIso m0,
              const MetricAnIso m1,
              const MetricAnIso m2);

  R2 mul(const R2 x) const { return R2(a11*x.x+a21*x.y,a21*x.x+a22*x.y); }
  R2 Orthogonal(const R2 x) { return R2(-(a21*x.x+a22*x.y),a11*x.x+a21*x.y); }
  R2 Orthogonal(const I2 x) { return R2(-(a21*x.x+a22*x.y),a11*x.x+a21*x.y); }
  MetricAnIso(      double      a,
              const MetricAnIso ma,
                    double      b,
              const MetricAnIso mb);

  int IntersectWith(const MetricAnIso M);
  MetricAnIso operator*(double c) const { double c2=c*c; return MetricAnIso(a11*c2,a21*c2,a22*c2); }
  MetricAnIso operator/(double c) const { double c2=1/(c*c); return MetricAnIso(a11*c2,a21*c2,a22*c2); } 
  double operator()(R2 x) const { return sqrt(x.x*x.x*a11+2*x.x*x.y*a21+x.y*x.y*a22); }
  double operator()(R2 x,R2 y) const { return x.x*y.x*a11+(x.x*x.y+x.y*y.x)*a21+x.y*y.y*a22; }
  inline void  Box(double &hx, double &hy) const;

  friend ostream& operator <<(      ostream&     f,
	                          const MetricAnIso& M)
  {
     f << " mtr a11 = " << M.a11 << " a21 = a12 = " << M.a21 << " a22 = " << M.a22 << ";";
     return f;
  }

  MetricAnIso(const MatVVP2x2);
};


class MatVVP2x2 
{
  friend class MetricAnIso;
  friend class MetricIso;

 public:
  double lambda1, lambda2;
  D2 v;

  MatVVP2x2(double r1, double r2, const D2 vp1)
     : lambda1(r1), lambda2(r2), v(vp1) {}

  void Abs() { lambda1 = bamg::Abs(lambda1), lambda2 = bamg::Abs(lambda2); }
  void pow(double p) { lambda1 = ::pow(lambda1,p); lambda2 = ::pow(lambda2,p); }
  void Min(double a) { lambda1 = bamg::Min(a,lambda1); lambda2 = bamg::Min(a,lambda2); }
  void Max(double a) { lambda1 = bamg::Max(a,lambda1); lambda2 = bamg::Max(a,lambda2); }
  void Minh(double h) { Max(1.0/(h*h)); }
  void Maxh(double h) { Min(1.0/(h*h)); }
  void Isotrope() { lambda1 = lambda2 = bamg::Max(lambda1,lambda2); }
  friend ostream& operator <<(ostream& f, const MatVVP2x2& c)
   { f << '{' << 1/sqrt(c.lambda1) << ',' << 1/sqrt(c.lambda2) << ','
       << c.v << '}' << std::flush; return f; }
  friend istream& operator >>(istream& f,  MatVVP2x2 & c)
  { f >> c.lambda1 >>c.lambda2 >> c.v.x >> c.v.y ;c.v /= Norme2(c.v); return f; }
  MatVVP2x2(const MetricAnIso );
  MatVVP2x2(const MetricIso M) : lambda1(1/(M.h*M.h)),lambda2(1/(M.h*M.h)),v(1,0) {}
  double hmin() const { return sqrt(1/bamg::Max3(lambda1,lambda2,1e-30)); }
  double hmax() const { return sqrt(1/bamg::Max(bamg::Min(lambda1,lambda2),1e-30)); }
  double lmax() const { return bamg::Max3(lambda1,lambda2,1e-30); }
  double lmin() const { return bamg::Max(bamg::Min(lambda1,lambda2),1e-30); }
  double Aniso2() const { return lmax()/lmin(); }
  inline void BoundAniso2(const double coef);
  double Aniso() const { return sqrt(Aniso2()); }
  void BoundAniso(const double c){ BoundAniso2(1/(c*c)); }
  void operator*=(double coef) { lambda1*=coef; lambda2*=coef; }
};


inline void MatVVP2x2::BoundAniso2(const double coef) 
{
   if (coef<=1.00000000001) 
      if (lambda1<lambda2)
         lambda1 = bamg::Max(lambda1,lambda2*coef);
      else
         lambda2 = bamg::Max(lambda2,lambda1*coef);
   else  // a verifier 
      if (lambda1>lambda2)
         lambda1 = bamg::Min(lambda1,lambda2*coef);
      else
         lambda2 = bamg::Min(lambda2,lambda1*coef);
}


void ReductionSimultanee(MetricAnIso M1,
                         MetricAnIso M2,
                         double&     l1,
                         double&     l2,
                         D2xD2&      V);

inline MetricAnIso::MetricAnIso(const MatVVP2x2 M)  
{
// recompose M in: M = V^t lambda V 
// V = (v,v^\perp) where v^\perp = (-v_1,v_0)
   double v00=M.v.x*M.v.x;
   double v11=M.v.y*M.v.y;
   double v01=M.v.x*M.v.y;
   a11 = v00*M.lambda1 + v11*M.lambda2;
   a21 = v01*(M.lambda1 - M.lambda2);
   a22 = v00*M.lambda2 + v11*M.lambda1;
}


inline void MetricAnIso::Box(double& hx,
                             double& hy) const 
{
   double d=a11*a22-a21*a21;
   hx = sqrt(a22/d);
   hy = sqrt(a11/d);
}


inline MetricIso::MetricIso(const MetricAnIso M)
{
   MatVVP2x2 vp(M);
   h = 1/sqrt(Max(vp.lambda1,vp.lambda2));
}


inline  MetricIso::MetricIso(double a11,
                             double a21,
                             double a22)
{
   MatVVP2x2 vp(MetricAnIso(a11,a21,a22));
   h = 1/sqrt(Max(vp.lambda1,vp.lambda2));
}


class SaveMetricInterpole {

  friend double LengthInterpole(const MetricAnIso, const MetricAnIso, R2);
  friend double abscisseInterpole(const MetricAnIso, const MetricAnIso, R2, double, int);
  int opt;
  double lab, L[1024], S[1024];
};

extern SaveMetricInterpole LastMetricInterpole;

double LengthInterpole(const MetricAnIso Ma,
                       const MetricAnIso Mb,
                             R2          AB);

double abscisseInterpole(const MetricAnIso Ma,
                         const MetricAnIso Mb,
                               R2          AB,
                               double      s,
                               int         optim=0);

inline double LengthInterpole(double la,
                              double lb) 
{
   return (Abs(la-lb)<1.0e-6*Max3(la,lb,1.0e-20)) ? (la+lb)/2 : la*lb*log(la/lb)/(la-lb);
}


inline double abscisseInterpole(double la,
                                double lb,
                                double lab,
                                double s)
{
   return (Abs(la-lb)<1.0e-6*Max3(la,lb,1.0e-20)) ? s : (exp(s*lab*(la-lb)/(la*lb))-1)*lb/(la-lb);
}

}
