/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                         Implementation of class 'Grid'

  ==============================================================================*/

#include "mesh/Grid.h"
#include "mesh/MeshUtil.h"
#include "linear_algebra/GraphOfMatrix.h"
#include "linear_algebra/Vect_impl.h"
#include "OFELIException.h"

using std::to_string;

namespace OFELI {

Grid::Grid()
     : _uniform(true), _dim(2), _nb_nodes(0), _nb_dof(1), _n(10), _cm(0), _cM(0)
{
}


Grid::Grid(real_t xm,
           real_t xM,
           size_t npx)
     : _uniform(true), _dim(1), _cm(0), _cM(0)
{
   _n.x = npx; _n.y = _n.z = 0;
   _h.y = _h.z = 0;
   if (xM <= xm)
      throw OFELIException("Grid::Grid(real_t,real_t,size_t): xmin = "+to_string(xm)+" is larger than xmax = "+to_string(xM));
   _xmin.x = xm; 
   _xmax.x = xM;
   _xmin.y = _xmax.y = _xmin.z = _xmax.z = 0;
   _h.x = (_xmax.x-_xmin.x)/_n.x;
   _active.resize(_n.x);
   _nb_nodes = _n.x + 1;
   for (size_t i=0; i<_n.x; i++)
      _active[i] = 1;
   _x.setSize(_n.x+1);
   for (size_t i=0; i<=_n.x; ++i)
      _x[i] = _xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x;
}


Grid::Grid(real_t xm,
           real_t xM,
           real_t ym,
           real_t yM,
           size_t npx,
           size_t npy)
     : _uniform(true), _dim(2), _nb_dof(1), _cm(0), _cM(0)
{
   _n.x = npx; _n.y = npy; _n.z = 0;
   _h.z = 0;
   if (xM <= xm)
      throw OFELIException("Grid::Grid(real_t,real_t,real_t,size_t,size_t): xmin = " + to_string(xm) + 
                           " is larger than xmax = " + to_string(xM));
   if (yM <= ym)
      throw OFELIException("Grid::Grid(real_t,real_t,real_t,size_t,size_t): ymin = " + to_string(ym) + 
                           " is larger than ymax = " + to_string(yM));
   _xmin.x = xm; _xmin.y = ym; _xmin.z = 0; 
   _xmax.x = xM; _xmax.y = yM; _xmax.z = 0;
   _h.x = (_xmax.x-_xmin.x)/_n.x;
   _h.y = (_xmax.y-_xmin.y)/_n.y;
   _active.resize(_n.x*_n.y);
   _nb_nodes = (_n.x+1)*(_n.y+1);
   for (size_t i=0; i<_n.x; i++)
      for (size_t j=0; j<_n.y; j++)
          _active[_n.y*i+j] = 1;
   _x.setSize(_n.x+1);
   _y.setSize(_n.y+1);
   for (size_t i=0; i<=_n.x; ++i)
      _x[i] = _xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x;
   for (size_t j=0; j<=_n.y; ++j)
      _y[j] = _xmin.y+(_xmax.y-_xmin.y)*(j-1)/_n.y;
}


Grid::Grid(Point<real_t> m,
           Point<real_t> M,
           size_t        npx,
           size_t        npy)
      : _uniform(true), _dim(2), _nb_dof(1), _cm(0), _cM(0)
{
   _n.x = npx; _n.y = npy; _n.z = 0;
   _h.z = 0;
   if (M.x <= m.x)
      throw OFELIException("Grid::Grid(Point<real_t>,Point<real_t>,size_t,size_t): xmin = " + to_string(M.x) + 
                           " is larger than xmax = " + to_string(m.x));
   if (M.y <= m.y)
      throw OFELIException("Grid::Grid(Point<real_t>,Point<real_t>,size_t,size_t): ymin = " + to_string(M.y) + 
                           " is larger than ymax = " + to_string(m.y));
   _xmin = m; _xmax = M;
   _xmin.z = _xmax.z = 0;
   _h.x = (_xmax.x-_xmin.x)/_n.x;
   _h.y = (_xmax.y-_xmin.y)/_n.y;
   _nb_nodes = (_n.x+1)*(_n.y+1);
   _active.resize(_n.x*_n.y);
   for (size_t i=0; i<_n.x; i++)
      for (size_t j=0; j<_n.y; j++)
         _active[_n.y*i+j] = 1;
   _x.setSize(_n.x+1);
   _y.setSize(_n.y+1);
   for (size_t i=0; i<=_n.x; ++i)
      _x[i] = _xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x;
   for (size_t j=0; j<=_n.y; ++j)
      _y[j] = _xmin.y+(_xmax.y-_xmin.y)*(j-1)/_n.y;
}


Grid::Grid(real_t xm,
           real_t xM,
           real_t ym,
           real_t yM,
           real_t zm,
           real_t zM, 
           size_t npx,
           size_t npy,
           size_t npz)
     : _uniform(true), _dim(3), _nb_dof(1), _cm(0), _cM(0)
{
   _n.x = npx; _n.y = npy; _n.z = npz;
   if (xM <= xm)
      throw OFELIException("Grid::Grid(real_t,real_t,real_t,real_t,real_t,real_t,size_t,size_t,size_t): xmin = " + 
                            to_string(xM) + " is larger than xmax = " + to_string(xm));
   if (yM <= ym)
      throw OFELIException("Grid::Grid(real_t,real_t,real_t,real_t,real_t,real_t,size_t,size_t,size_t): ymin = " + 
                           to_string(yM) + " is larger than ymax = " + to_string(ym));
   if (zM <= zm)
      throw OFELIException("Grid::Grid(real_t,real_t,real_t,real_t,real_t,real_t,size_t,size_t,size_t): zmin = " + 
                            to_string(zM) + " is larger than zmax = " + to_string(zm));
   _xmin.x = xm; _xmin.y = ym; _xmin.z = zm; 
   _xmax.x = xM; _xmax.y = yM; _xmax.z = zM;
   _h.x = (_xmax.x-_xmin.x)/_n.x;
   _h.y = (_xmax.y-_xmin.y)/_n.y;
   _h.z = (_xmax.z-_xmin.z)/_n.z;
   _nb_nodes = (_n.x+1)*(_n.y+1)*(_n.z+1);
   _active.resize(_n.x*_n.y*_n.z);
   for (size_t i=0; i<_n.x; i++)
      for (size_t j=0; j<_n.y; j++)
         for (size_t k=0; k<_n.z; k++)
            _active[_n.y*_n.z*i+_n.z*j+k] = 1;
   _x.setSize(_n.x+1);
   _y.setSize(_n.y+1);
   _z.setSize(_n.z+1);
   for (size_t i=0; i<=_n.x; ++i)
      _x[i] = _xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x;
   for (size_t j=0; j<=_n.y; ++j)
      _y[j] = _xmin.y+(_xmax.y-_xmin.y)*(j-1)/_n.y;
   for (size_t k=0; k<=_n.z; ++k)
      _z[k] = _xmin.z+(_xmax.z-_xmin.z)*(k-1)/_n.z;
}


Grid::Grid(Point<real_t> m,
           Point<real_t> M,
           size_t        npx,
           size_t        npy,
           size_t        npz)
     : _uniform(true), _dim(3), _nb_dof(1), _cm(0), _cM(0)
{
   _n.x = npx; _n.y = npy; _n.z = npz;
   if (npz==0)
      _dim = 2;
   if (npy==0)
      _dim = 1;
   if (M.x <= m.x)
      throw OFELIException("Grid::Grid(Point<real_t>,Point<real_t>,size_t,size_t,size_t): xmin = " + to_string(M.x) + 
                           " is larger than xmax = " + to_string(m.x));
   if (M.y <= m.y)
      throw OFELIException("Grid::Grid(Point<real_t>,Point<real_t>,size_t,size_t,size_t): ymin = " + to_string(M.y) + 
                           " is larger than ymax = " + to_string(m.y));
   if (M.z <= m.z)
      throw OFELIException("Grid::Grid(Point<real_t>,Point<real_t>,size_t,size_t,size_t): zmin = " + to_string(M.z) + 
                           " is larger than zmax = " + to_string(m.z));
   _xmin = m; _xmax = M;
   _h.x = (_xmax.x-_xmin.x)/_n.x;
   _h.y = (_xmax.y-_xmin.y)/_n.y;
   _h.z = (_xmax.z-_xmin.z)/_n.z;
   _nb_nodes = (_n.x+1)*(_n.y+1)*(_n.z+1);
   _active.resize(_n.x*_n.y*_n.z);
   for (size_t i=0; i<_n.x; i++)
      for (size_t j=0; j<_n.y; j++)
         for (size_t k=0; k<_n.z; k++)
            _active[_n.y*_n.z*i+_n.z*j+k] = 1;
   _x.setSize(_n.x+1);
   _y.setSize(_n.y+1);
   _z.setSize(_n.z+1);
   for (size_t i=0; i<=_n.x; ++i)
      _x[i] = _xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x;
   for (size_t j=0; j<=_n.y; ++j)
      _y[j] = _xmin.y+(_xmax.y-_xmin.y)*(j-1)/_n.y;
   for (size_t k=0; k<=_n.z; ++k)
      _z[k] = _xmin.z+(_xmax.z-_xmin.z)*(k-1)/_n.z;
}


void Grid::setX(const Vect<real_t>& x)
{
   _x = x;
   _uniform = false;
}


void Grid::setXY(const Vect<real_t>& x,
                 const Vect<real_t>& y)
{
   _x = x;
   _y = y;
   _uniform = false;
}


void Grid::setXYZ(const Vect<real_t>& x,
                  const Vect<real_t>& y,
                  const Vect<real_t>& z)
{
   _x = x;
   _y = y;
   _z = z;
   _uniform = false;
}


void Grid::setNbDOF(size_t n)
{
   _nb_dof = n;
}


void Grid::setN(size_t nx,
                size_t ny,
                size_t nz)
{
   _dim = 3;
   _n.x = nx; _n.y = ny; _n.z = nz;
   _h.x = (_xmax.x-_xmin.x)/_n.x;
   _h.y = (_xmax.y-_xmin.y)/_n.y;
   _h.z = (_xmax.z-_xmin.z)/_n.z;
   if (nz==0) {
      _dim = 2;
      _active.resize(_n.x*_n.y);
      for (size_t i=0; i<_n.x; i++)
         for (size_t j=0; j<_n.y; j++)
            _active[_n.y*i+j] = 1;
      _nb_nodes = (_n.x+1)*(_n.y+1);
   }
   if (ny==0) {
      _dim = 1;
      _active.resize(_n.x);
      for (size_t i=0; i<_n.x; i++)
          _active[i] = 1;
      _nb_nodes = (_n.x+1);
   }
   if (_dim==3) {
      _active.resize(_n.x*_n.y*_n.z);
      for (size_t i=0; i<_n.x; i++)
         for (size_t j=0; j<_n.y; j++)
            for (size_t k=0; k<_n.z; k++)
               _active[_n.y*_n.z*i+_n.z*j+k] = 1;
      _nb_nodes = (_n.x+1)*(_n.y+1)*(_n.z+1);
   }
   if (_h==0)
      throw OFELIException("Grid::setN(size_t,size_t,size_t): Grid size must be given first.");
}


void Grid::setDomain(real_t xmin,
                     real_t xmax)
{
   _xmin.x = xmin; _xmin.y = _xmin.z = 0;
   _xmax.x = xmax; _xmax.y = _xmax.z = 0;
   _h = 0;
   if (_n.x>0)
      _h.x = (_xmax.x-_xmin.x)/_n.x;
}


void Grid::setDomain(real_t xmin,
                     real_t xmax,
                     real_t ymin,
                     real_t ymax)
{
   _dim = 2;
   _xmin.x = xmin; _xmax.x = xmax;
   _xmin.y = ymin; _xmax.y = ymax;
   _xmin.z = _xmax.z = 0;
   _h = 0;
   if (_n.x>0)
      _h.x = (_xmax.x-_xmin.x)/_n.x;
   if (_n.y>0)
      _h.y = (_xmax.y-_xmin.y)/_n.y;
}


void Grid::setDomain(real_t xmin,
                     real_t xmax,
                     real_t ymin,
                     real_t ymax,
                     real_t zmin,
                     real_t zmax)
{
   _dim = 3;
   _xmin.x = xmin; _xmax.x = xmax;
   _xmin.y = ymin; _xmax.y = ymax;
   _xmin.z = zmin; _xmax.z = zmax;
   _h = 0;
   if (_n.x>0)
      _h.x = (_xmax.x-_xmin.x)/_n.x;
   if (_n.y>0)
      _h.y = (_xmax.y-_xmin.y)/_n.y;
   if (_n.z>0)
      _h.z = (_xmax.z-_xmin.z)/_n.z;
}


void Grid::setDomain(Point<real_t> xmin,
                     Point<real_t> xmax)
{
   _dim = 3;
   if (_n.z==0)
      _dim = 2;
   if (_n.y==0)
      _dim = 1;
   _xmin = xmin; _xmax = xmax;
   _h = 0;
   if (_n.x>0)
      _h.x = (_xmax.x-_xmin.x)/_n.x;
   if (_n.y>0)
      _h.y = (_xmax.y-_xmin.y)/_n.y;
   if (_n.z>0)
      _h.z = (_xmax.z-_xmin.z)/_n.z;
}


real_t Grid::getX(size_t i) const
{
   return (_xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x);
}


real_t Grid::getY(size_t j) const
{
   return (_xmin.y+(_xmax.y-_xmin.y)*(j-1)/_n.y);
}


real_t Grid::getZ(size_t k) const
{
   return (_xmin.z+(_xmax.z-_xmin.z)*(k-1)/_n.z);
}


Point2D<real_t> Grid::getXY(size_t i,
                            size_t j) const
{
   return Point2D<real_t>(getX(i),getY(j));
}


Point<real_t> Grid::getXYZ(size_t i,
                           size_t j,
                           size_t k) const
{
   return Point<real_t>(getX(i),getY(j),getZ(k));
}


Point<real_t> Grid::getCoord(size_t i) const
{
   return Point<real_t>(_xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x,0,0);
}


Point<real_t> Grid::getCoord(size_t i,
                             size_t j) const
{
   return Point<real_t>(_xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x,
                        _xmin.y+(_xmax.y-_xmin.y)*(j-1)/_n.y,0);
}
    

Point<real_t> Grid::getCoord(size_t i,
                             size_t j,
                             size_t k) const
{
   return Point<real_t>(_xmin.x+(_xmax.x-_xmin.x)*(i-1)/_n.x,
                        _xmin.y+(_xmax.y-_xmin.y)*(j-1)/_n.y,
                        _xmin.z+(_xmax.z-_xmin.z)*(k-1)/_n.z);
}
 

real_t Grid::getCenter(size_t i) const
{
   return 0.5*(getX(i)+getX(i+1));
}
 

Point<real_t> Grid::getCenter(size_t i,
                              size_t j) const
{
  Point<real_t> x1(getXY(i  ,j  ).x,getXY(i  ,j  ).y),
                x2(getXY(i+1,j  ).x,getXY(i+1,j  ).y),
                x3(getXY(i+1,j+1).x,getXY(i+1,j+1).y),
                x4(getXY(i  ,j+1).x,getXY(i  ,j+1).y);
   return 0.25*(x1+x2+x3+x4);
}
 

Point<real_t> Grid::getCenter(size_t i,
                              size_t j,
                              size_t k) const
{
   return 0.125*(getXYZ(i,j,k  )+getXYZ(i+1,j,k  )+getXYZ(i+1,j+1,k  )+getXYZ(i,j+1,k  ) +
                 getXYZ(i,j,k+1)+getXYZ(i+1,j,k+1)+getXYZ(i+1,j+1,k+1)+getXYZ(i,j+1,k+1));
}


void Grid::setCode(int side,
                   int code)
{
   if (side==MIN_X)
      _cm.x = code;
   else if (side==MAX_X)
      _cM.x = code;
   else if (side==MIN_Y)
      _cm.y = code;
   else if (side==MAX_Y)
      _cM.y = code;
   else if (side==MIN_Z)
      _cm.z = code;
   else if (side==MAX_Z)
      _cM.z = code;
   else
      ;
}


int Grid::getCode(int side) const
{
   if (side==MIN_X)
      return _cm.x;
   else if (side==MAX_X)
      return _cM.x;
   else if (side==MIN_Y)
      return _cm.y;
   else if (side==MAX_Y)
      return _cM.y;
   else if (side==MIN_Z)
      return _cm.z;
   else if (side==MAX_Z)
      return _cM.z;
   else
      return 0;
}


int Grid::getCode(size_t i,
                  size_t j) const
{
   if (i==0)
      return _cm.x;
   else if (i==_n.x)
      return _cM.x;
   else if (j==0)
      return _cm.y;
   else if (j==_n.y)
      return _cM.y;
   else
      return 0;
}


int Grid::getCode(size_t i,
                  size_t j,
                  size_t k) const
{
   if (i==0)
      return _cm.x;
   else if (i==_n.x)
      return _cM.x;
   else if (j==0)
      return _cm.y;
   else if (j==_n.y)
      return _cM.y;
   else if (k==0)
      return _cm.z;
   else if (k==_n.z)
      return _cM.z;
   else
      return 0;
}


void Grid::Deactivate(size_t i)
{
   if (i)
      _active[i-1] = 0;
   else {
      for (size_t ii=0; ii<_n.x; ii++)
         _active[ii] = 0;
   }
}


void Grid::Deactivate(size_t i,
                      size_t j)
{
   if (i && j)
      _active[_n.y*(i-1)+j-1] = 0;
   else if (i && !j) {
      for (size_t jj=1; jj<=_n.y; jj++)
         _active[_n.y*(i-1)+jj-1] = 0;
   }
   else if (!i && j) {
      for (size_t ii=1; ii<=_n.x; ii++)
         _active[_n.y*(ii-1)+j-1] = 0;
   }
   else if (!i && !j) {
      for (size_t ii=1; ii<=_n.x; ii++)
         for (size_t jj=1; jj<=_n.y; jj++)
            _active[_n.y*(ii-1)+jj-1] = 0;
   }
}


void Grid::Deactivate(size_t i,
                      size_t j,
                      size_t k)
{
   if (i && j && k) 
     _active[_n.y*_n.z*(i-1)+_n.z*(j-1)+k-1] = 0;
   else if (!i) {
      if (!j && k) {
         for (size_t ii=1; ii<=_n.x; ii++)
            for (size_t jj=1; jj<=_n.y; jj++)
               _active[_n.y*_n.z*(ii-1)+_n.z*(jj-1)+k-1] = 0;
      }
      else if (j && !k) {
         for (size_t ii=1; ii<=_n.x; ii++)
            for (size_t kk=1; kk<=_n.z; kk++)
               _active[_n.y*_n.z*(ii-1)+_n.z*(j-1)+kk-1] = 0;
      }
      else if (j && k) {
         for (size_t ii=1; ii<=_n.x; ii++)
            _active[_n.y*_n.z*(ii-1)+_n.z*(j-1)+k-1] = 0;
      }
   }
   else if (!j) {
      if (i && !k) {
         for (size_t jj=1; jj<=_n.y; jj++)
            for (size_t kk=1; kk<=_n.z; kk++)
               _active[_n.y*_n.z*(i-1)+_n.z*(jj-1)+kk-1] = 0;
      }
      if (i && k) {
         for (size_t jj=1; jj<=_n.y; jj++)
            _active[_n.y*_n.z*(i-1)+_n.z*(jj-1)+k-1] = 0;
      }
   }
   else if (!k) {
      if (i && j) {
         for (size_t kk=1; kk<=_n.z; kk++)
            _active[_n.y*_n.z*(i-1)+_n.z*(j-1)+kk-1] = 0;
      }
   }
}


ostream& operator<<(ostream&    s,
                    const Grid& g)
{
   size_t dim = 3;
   if (g.getNz()==0)
      dim = 2;
   if (g.getNy()==0)
      dim = 1;
   s << "\n\nS T R U C T U R E D    G R I D    D A T A\n";
   s << "=========================================\n\n";
   s << "Space Dimension:       " << setw(6) << dim << endl;
   s << "x-domain:              [" << g.getXMin().x << "," << g.getXMax().x << "]" << endl;
   if (dim > 1)
      s << "y-domain:              [" << g.getXMin().y << "," << g.getXMax().y << "]" << endl;
   if (dim > 2)
      s << "z-domain:              [" << g.getXMin().z << "," << g.getXMax().z << "]" << endl;
   s << "Number of x-intervals: " << setw(6) << g.getNx() << endl;
   if (dim > 1)
      s << "Number of y-intervals: " << setw(6) << g.getNy() << endl;
   if (dim > 2)
      s << "Number of z-intervals: " << setw(6) << g.getNz() << endl;
   s << endl;
   return s;
}

} /* namespace OFELI */
