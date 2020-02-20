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

#ifndef __R2_H
#define __R2_H

namespace bamg {

template <class R,class RR> class P2xP2;

template <class R,class RR>
class P2 {
  
 public:  
   R x,y;

   P2() : x(0),y(0) { }

   P2(R a, R b) : x(a),y(b) { }

   P2(P2 A, P2 B) : x(B.x-A.x),y(B.y-A.y) { }

   P2<R,RR> operator+(const P2<R,RR> & cc) const
   { return P2<R,RR>(x+cc.x,y+cc.y); }

   P2<R,RR> operator-(const P2<R,RR>& cc) const
   { return P2<R,RR>(x-cc.x,y-cc.y); }

   P2<R,RR> operator-() const
   { return P2<R,RR>(-x,-y); }

   RR operator,(const P2<R,RR>& cc) const
   { return  (RR)x*(RR)cc.x + (RR)y*(RR)cc.y; }

   P2<R,RR> operator*(R cc) const
   { return P2<R,RR>(x*cc,y*cc); }

   P2<R,RR> operator/(R cc) const
   { return P2<R,RR>(x/cc,y/cc); }

   P2<R,RR> operator+=(const P2<R,RR>& cc)
   { x += cc.x; y += cc.y; return *this; }

   P2<R,RR> operator/=(const R r)
   { x /= r;y /= r; return *this; }

   P2<R,RR> operator*=(const R r)
   { x *= r;y *= r; return *this; }

   P2<R,RR>  operator-=(const P2<R,RR>& cc)
   { x -= cc.x; y -= cc.y; return *this; }

};


template <class R,class RR>
class P2xP2 {

  friend ostream& operator <<(ostream& f, const P2xP2<R,RR>& c)
     { f << '[' << c.x << ',' << c.y << ']' << std::flush; return f; }

  friend P2<R,RR> operator*(P2<R,RR> c,P2xP2<R,RR> cc)
     { return P2<R,RR>(c.x*cc.x.x + c.y*cc.y.x, c.x*cc.x.y + c.y*cc.y.y); }

 public:

  P2<R,RR> x, y;

  P2xP2(): x(),y()  {}

  P2xP2(P2<R,RR> a, P2<R,RR> b) : x(a), y(b) {}

  P2xP2 (P2<R,RR> a, P2<R,RR> b, P2<R,RR> c) : x(b-a), y(c-a) {}

  P2xP2 (R xx,R xy,R yx,R yy) : x(xx,xy),y(yx,yy) {}

  P2<R,RR> operator*(const P2<R,RR> c) const 
  { return P2<R,RR>(x.x*c.x + x.y*c.y, y.x*c.x + y.y*c.y); }

  P2xP2<R,RR> operator*(P2xP2<R,RR> c) const 
    { return  P2xP2<R,RR>(x.x*c.x.x + x.y*c.y.x,
			  x.x*c.x.y + x.y*c.y.y,
			  y.x*c.x.x + y.y*c.y.x,
			  y.x*c.x.y + y.y*c.y.y);
    }

  RR det() const
  { return (RR)x.x*(RR)y.y - (RR)x.y*(RR)y.x; }

  P2xP2<R,RR> inv() const
  { RR d = (*this).det(); 
    return P2xP2<R,RR>((R)(y.y/d) ,(R)(-x.y/d),(R)(-y.x/d) ,(R)(x.x/d));
  };
  
  P2xP2<R,RR>t()
  { return P2xP2<R,RR>(x.x,y.x,x.y,y.y); }

  P2<R,RR>tx()
  { return P2<R,RR>(x.x,y.x); }

  P2<R,RR>ty()
  {return P2<R,RR>(x.y,y.y); }
};  


template<class R,class RR>  
inline RR Det(const P2<R,RR> x,
              const P2<R,RR> y)
{
   return (RR)x.x*(RR)y.y - (RR) x.y*(RR)y.x;
}


template<class R,class RR>  
inline RR Area2(const P2<R,RR>& a,
                const P2<R,RR>& b,
                const P2<R,RR>& c)
{
   return Det(b-a,c-a);
}


template<class R,class RR>  
inline R Norme1(const P2<R,RR>& x)
{
   return (Abs(x.x)+Abs(x.y));
}


template <class R,class RR>  
inline R NormeInfini(const P2<R,RR>& x)
{
   return Max(Abs(x.x),Abs(x.y));
} 


template<class R,class RR>
inline RR Norme2_2(const P2<R,RR>& x)
{
   return (RR)x.x*(RR)x.x + (RR)x.y*(RR)x.y;
}


template<class R,class RR> 
inline RR Norme2(const P2<R,RR>& x)
{
   return sqrt((RR)x.x*(RR)x.x + (RR)x.y*(RR)x.y);
}


template<class R,class RR>  
inline P2<R,RR> Orthogonal(const P2<R,RR>& x)
{
   return P2<R,RR>(-x.y,x.x);
}


template<class R,class RR>
inline ostream& operator<<(      ostream&  f,
                           const P2<R,RR>& c)
{
   f << '[' << c.x << ',' << c.y <<']' << std::flush;
   return f;
}

}

#endif
