/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2019 Rachid Touzani

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

                      Implementation of class 'FastMarching2D'

 ==============================================================================*/

#include "equations/interface/FastMarching2D.h"
#include <algorithm>

namespace OFELI {

FastMarching2D::FastMarching2D()
{ 
   _bw = 1;
   _ext = false;
}


FastMarching2D::FastMarching2D(const Grid&   g,
                               Vect<real_t>& ls)
{
   _v = new Vect<real_t>(g.getNx()+1,g.getNy()+1);
   *_v = 1;
   _A = &ls;
   _cg = &g;
   _ext = false;
   Init();
}


FastMarching2D::FastMarching2D(const Grid&   g,
                               Vect<real_t>& ls,
                               Vect<real_t>& F)
{
   _v = &F;
   _A = &ls;
   _cg = &g;
   _ext = true;
   Init();
}


FastMarching2D::~FastMarching2D()
{
   if (!_ext)
      delete _v;
}


void FastMarching2D::Init()
{
   _bw = 1;
   _hx = (_cg->getXMax().x-_cg->getXMin().x)/_cg->getNx();
   _hy = (_cg->getXMax().y-_cg->getXMin().y)/_cg->getNy();
   _g.setXMin(_cg->getXMin()-2*Point<real_t>(_hx,_hy));
   _g.setXMax(_cg->getXMax()+2*Point<real_t>(_hx,_hy));
   _nx = _cg->getNx() + 4;
   _ny = _cg->getNy() + 4;
   _g.setN(_nx,_ny);
   _heap.setSize((_nx+1)*(_ny+1));
   _pos.setSize(_nx+1,_ny+1);
   _sol.setSize(_nx+1,_ny+1);
   _ls.setSize(_nx+1,_ny+1);
   _f.setSize(_nx+1,_ny+1);
   size_t i, j;
   for (i=3; i<=_nx-1; i++) {
      for (j=3; j<=_ny-1; j++) {
         _ls(i,j) = (*_A)(i-2,j-2);
         _f(i,j) = (*_v)(i-2,j-2);
      }
   }
   for (i=3; i<=_nx-1; i++) {
      _ls(i,1) = _ls(i,2) = (*_A)(i-2,1);
      _ls(i,_ny) = _ls(i,_ny+1) = (*_A)(i-2,_ny-3);
   }
   for (j=3; j<=_ny-1; j++) {
      _ls(1,j) = _ls(2,j) = (*_A)(1,j-1);
      _ls(_nx+1,j) = _ls(_nx,j) = (*_A)(_nx-3,j-2);
   }
   _ls(1,1) = _ls(2,2) = _ls(2,1) = _ls(1,2) = (*_A)(1,1);
   _ls(_nx+1,1) = _ls(_nx,2) = _ls(_nx,1) = _ls(_nx+1,2) = (*_A)(_nx-3,1);
   _ls(1,_ny+1) = _ls(2,_ny) = _ls(1,_ny) = _ls(2,_ny+1) = (*_A)(1,_ny-3);
   _ls(_nx+1,_ny+1) = _ls(_nx,_ny) = _ls(_nx+1,_ny) = _ls(_nx,_ny+1) = (*_A)(_nx-3,_ny-3);
   _u.setSize(_nx+1,_ny+1);
   _u = _f;
}


void FastMarching2D::execute()
{
   size_t i, j;
   Vect<real_t> bd(_ls);
   Pt mini;
   _current = 0;
   _sol = -1;

// For each point corresponding to the initial condition, we fix definitely
// its value to black and its neighbors to gray
   _f = _ls;
   ExtendLocalVelocity(bd);

   for (i=1; i<=_nx+1; i++) {
      for (j=1; j<=_ny+1; j++) {
         if (bd(i,j)>-1 && bd(i,j)<1)
            _sol(i,j) = fabs(bd(i,j));
         else
            _f(i,j) =  0;
      }
   }

   for (i=1; i<=_nx+1; i++) {
      for (j=1; j<=_ny+1; j++) {
         if (_sol(i,j)>=0)
            setGray(i,j);
      }
   }

   while (_current>0) {

//    Choose the smallest value of u among the gray ones
      mini = Delete();
      i = mini.ix;
      j = mini.iy;
//    Set it (definitely) to black
      _sol(i,j) = mini.val;
//    Calculate its neighbor values and set them to gray 
      int d1=0, d2=0;
      if (i==1) {
         if (_sol(i+1,j)>=0)
            d1 = 1;
      }
      else if (i==_nx+1) {
         if (_sol(i-1,j)>=0)
            d1 = -1;
      }
      else {
         if (_sol(i-1,j)>=0)
            d1 = -1;
         if ((_sol(i+1,j)>=0 && _sol(i+1,j)<_sol(i-1,j)) || _sol(i-1,j)<0)
            d1 = 1;
      }
      if (j==1) {
         if (_sol(i,j+1)==0)
            d2 = 1;
      }
      else if (j==_ny+1) {
         if (_sol(i,j-1)>=0)
            d2 = -1;
      }
      else {
         if (_sol(i,j-1)>=0)
            d2 = -1;
         if ((_sol(i,j+1)>=0 && _sol(i,j+1)<_sol(i,j-1)) || _sol(i,j-1)<0)
            d2 = 1;
      }
      if (d1==0)
         _f(i,j) = _f(i,j+d2);
      else if (d2==0)
         _f(i,j) = _f(i+d1,j);
      else
         _f(i,j) = (_f(i+d1,j)*(_sol(i,j)-_sol(i+d1,j))+_f(i,j+d2)*(_sol(i,j)-_sol(i,j+d2)))/
                   ((_sol(i,j)-_sol(i+d1,j))+(_sol(i,j)-_sol(i,j+d2)));
      setGray(i,j);
   }
   for (i=1; i<=_nx+1; i++) {
      for (j=1; j<=_ny+1; j++) {
         if (_ls(i,j)>0)
            _ls(i,j) =  _sol(i,j);
         else
            _ls(i,j) = -_sol(i,j);
      }
   }

   for (i=3; i<=_nx-1; i++) {
      for (j=3; j<=_ny-1; j++) {
         (*_A)(i-2,j-2) = _ls(i,j);
         (*_v)(i-2,j-2) = _f(i,j);
      }
   }
}


void FastMarching2D::Check()
{
   Mesh ms(*_cg,TRIANGLE);
   real_t dudx, dudy, errL2=0, errMax=0;
   size_t n = 1;
   for (size_t i=1; i<=_nx-4; i++) {
      for (size_t j=1; j<=_ny-4; j++) {
         dudx = 0.5*((*_A)(i+1,j)-(*_A)(i,j)+(*_A)(i+1,j+1)-(*_A)(i,j+1))/_hx;
         dudy = 0.5*((*_A)(i,j+1)-(*_A)(i,j)+(*_A)(i+1,j+1)-(*_A)(i+1,j))/_hy;
         errL2 += fabs(dudx*dudx + dudy*dudy - 1);
         errMax = std::max(errMax,fabs(dudx*dudx + dudy*dudy - 1));
         n++;
      }
   }
   Vect<real_t> phi(ms,*_A);
   ms.put("phi.m");
   errL2 /= (_nx-4)*(_ny-4);
   cout << endl << "Errors on distance: " << errL2 << "  " << errMax << endl;
}


void FastMarching2D::Add(real_t Uij,
                         size_t a,
                         size_t b)
{
   _current++;
   _heap[0].val = Uij;      // Sentinel
   size_t i = _current;	    // Position of current element not yet filled 
   size_t j = i/2;          // Father of current element
   while (_heap[j].val>Uij) {
      _heap[i] = _heap[j];
      _pos(_heap[j].ix,_heap[j].iy) = i;
      i = j;
      j /= 2;
   }

// Here _heap[j].val <= Uij and the heap[i] is free => inserted at i
   _heap[i].val = Uij;
   _heap[i].ix = a;
   _heap[i].iy = b;
   _pos(a,b) = i;
}


void FastMarching2D::Modify(real_t Uij,
                            size_t a,
                            size_t b)
{
   _heap[0].val = Uij;
   size_t i=_pos(a,b), j=i/2;
   while (_heap[j].val>Uij) {
      _heap[i] = _heap[j];
      _pos(_heap[j].ix,_heap[j].iy) = i;
      i = j;
      j /= 2;
   }
// Here A[j].val <= Uij and A[i] is free => Insert at i
   _heap[i].val = Uij;
   _heap[i].ix = a;
   _heap[i].iy = b;
   _pos(a,b) = i;
}


Pt FastMarching2D::Delete()
{
   Pt aux, ret;
   ret = _heap[1];
   aux = _heap[_current--];
   int i=1, j=2;
   bool finished = (j>_current);
   while (!finished) {
//    If i has exactly 2 sons at j and j+1, look for the one having the smallest value
      if (j<_current) {                          
         if (_heap[j].val>_heap[j+1].val)
            j++;
      }
      if (aux.val<=_heap[j].val)
         finished = true;
      else {
//       The position j is free
         _heap[i] = _heap[j];
         _pos(_heap[j].ix,_heap[j].iy) = i;
         i = j;
         j = 2*i;
         finished = (j>_current);
      }
   }
   _heap[i] = aux;
   _pos(_heap[i].ix,_heap[i].iy) = i;
   return ret;
}


void FastMarching2D::setGray(size_t i,
                             size_t j)
{
   size_t max_x, max_y;
   real_t a, b, c, ug=0;
   bool ps;

   if (i>=3 && i<=_nx-1 && j>=3 && j<=_ny-1) {
      if (_sol(i-1,j)<0 && i>=2 && j>=2 && j<=_ny) {
//       New origin: (i-1,j)
         max_x = 20;                           // max(D-x,-D+x,0)=-D+x
//       Unless the left neighbor is black and its value is smaller than the right one's
         if (_sol(i-2,j)>=0 && _sol(i,j)>_sol(i-2,j))
            max_x = 10;                        // max(D-x,-D+x,0)=D-x

//       If north and south neighbors are not black
         if (_sol(i-1,j+1)<0 && _sol(i-1,j-1)<0) 
            max_y = 3;                         // max(D-y,-D+y,0)=0
//       If north and south neighbors are black
         else if (_sol(i-1,j+1)>=0 && _sol(i-1,j-1)>=0) {
//          Choose the min.
            if (_sol(i-1,j-1)>_sol(i-1,j+1))
               max_y = 2;                      // max(D-y,-D+y,0)=-D+y
            else
               max_y = 1;                      // max(D-y,-D+y,0)=D-y
         }
//       If the north neighbor is black and not the south one
         else if (_sol(i-1,j+1)>=0)
            max_y = 2;                         // max(D-y,-D+y,0)=-D+y
//          If the South neighbor is black and not the north one
         else
            max_y = 1;                         // max(D-y,-D+y,0)=D-y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i-2,j) - 2*_sol(i-1,j-1);
               c = _sol(i-2,j)*_sol(i-2,j) + _sol(i-1,j-1)*_sol(i-1,j-1) - _hx*_hy;
               break;

            case 12:
               a = 2;
               b = -2*_sol(i-2,j) - 2*_sol(i-1,j+1);
               c = _sol(i-2,j)*_sol(i-2,j) + _sol(i-1,j+1)*_sol(i-1,j+1) - _hx*_hy;
               break;

            case 13:
               a = 1;
               b = -2*_sol(i-2,j);
               c = _sol(i-2,j)*_sol(i-2,j) - _hx*_hy;
               break;

            case 21:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i-1,j-1);
//               cout<<"_sol: "<<_sol(i,j)<<" "<<_sol(i-1,j-1)<<endl;
               c = _sol(i,j)*_sol(i,j) + _sol(i-1,j-1)*_sol(i-1,j-1) - _hx*_hy;
               break;

            case 22:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i-1,j+1);
               c = _sol(i,j)*_sol(i,j) + _sol(i-1,j+1)*_sol(i-1,j+1) - _hx*_hy;
               break;

            case 23:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            default:
               ps = false;
               break;
         }

         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i-1,j) was not gray
         if (_sol(i-1,j)!=-2) {
            Add(ug,i-1,j);
            _sol(i-1,j) = -2;
         }
         else
//          If the new value is smaller than the previous, modify it
            if (ug<_heap[_pos(i-1,j)].val)
               Modify(ug,i-1,j);
      }

//    If its right neighbor is not already black and lies in the domain
      if (_sol(i+1,j)<0 && i<=_nx-1 && j>=2 && j<=_ny) {
         max_x = 10;                                   // max(D-x,-D+x,0)=D-x
//       Unless the right neighbor is black and its value is smaller than the left one's
         if (_sol(i+2,j)>=0 && _sol(i,j)>_sol(i+2,j))
            max_x = 20;                                // max(D-x,-D+x,0)=-D+x
//       If the north and south neighbors are not black
         if (_sol(i+1,j+1)<0 && _sol(i+1,j-1)<0)
            max_y = 3;                                 // max(D-y,-D+y,0)=0
//       If the north and south neighbors are black
         else if (_sol(i+1,j+1)>=0 && _sol(i+1,j-1)>=0)
//          We choose the minimum
            if (_sol(i+1,j-1)>_sol(i+1,j+1))
               max_y = 2;                             // max(D-y,-D+y,0)=-D+y
            else
               max_y = 1;                             // max(D-y,-D+y,0)=D-y
//       If the north neighbor is black and not the south one
         else if (_sol(i+1,j+1)>=0)
            max_y = 2;                                // max(D-y,-D+y,0)=-D+y
//          If the south neighbor is black and not the north one
         else
            max_y = 1;                                // max(D-y,-D+y,0)=D-y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i+1,j-1);
               c = _sol(i,j)*_sol(i,j) + _sol(i+1,j-1)*_sol(i+1,j-1) - _hx*_hy;
               break;

            case 12:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i+1,j+1);
               c = _sol(i,j)*_sol(i,j) + _sol(i+1,j+1)*_sol(i+1,j+1) - _hx*_hy;
               break;

            case 13:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            case 21:
               a = 2;
               b = -2*_sol(i+2,j) - 2*_sol(i+1,j-1);
               c = _sol(i+2,j)*_sol(i+2,j) + _sol(i+1,j-1)*_sol(i+1,j-1) - _hx*_hy;
               break;

            case 22:
               a = 2;
               b = -2*_sol(i+2,j) - 2*_sol(i+1,j+1);
               c = _sol(i+2,j)*_sol(i+2,j) + _sol(i+1,j+1)*_sol(i+1,j+1) - _hx*_hy;
               break;

            case 23:
               a = 1;
               b = -2*_sol(i+2,j);
               c = _sol(i+2,j)*_sol(i+2,j) - _hx*_hy;
               break;

            default:
               ps = false;
               break;
         }
         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i+1,j) was not gray
         if (_sol(i+1,j)!=-2) {
            Add(ug,i+1,j);
            _sol(i+1,j) = -2;
         }
         else {
//          If the new value is smaller than the previus one, modify it
            if (ug<_heap[_pos(i+1,j)].val)
               Modify(ug,i+1,j);
         }
      }

//    If north neighbor is not already black and lies in the domain
      if (_sol(i,j+1)<0 && i>=2 && i<=_nx && j<=_ny-1) {
//       New origin: (i,j+1) 
//       If right and left neighbors are not black
         if (_sol(i+1,j+1)<0 && _sol(i-1,j+1)<0)
            max_x = 30;                               // max(D-x,-D+x,0)=0
//       If right and left neighbors are black
         else if (_sol(i+1,j+1)>=0 && _sol(i-1,j+1)>=0) {
//          Choose the minimum
            if (_sol(i+1,j+1)>_sol(i-1,j+1))
               max_x = 10;                            // max(D-x,-D+x,0)=D-x
            else
               max_x = 20;                            // max(D-x,-D+x,0)=-D+x
         }

//       If right neighbor is black and not the left one
         else if (_sol(i+1,j+1)>=0)
            max_x = 20;                               // max(D-x,-D+x,0)=-D+x
//       If left neighbor is black and not the right one
         else
            max_x = 10;                               // max(D-x,-D+x,0)=D-x
         max_y = 1;                                   // max(D-y,-D+y,0)=D-y
//       Unless the north neighbor is black and its value is smaller than the south one's
         if (_sol(i,j+2)>=0 && _sol(i,j)>_sol(i,j+2))
            max_y = 2;                                // max(D-y,-D+y,0)=-D+y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i-1,j+1) - 2*_sol(i,j);
               c = _sol(i-1,j+1)*_sol(i-1,j+1) + _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            case 12:
               a = 2;
               b = -2*_sol(i-1,j+1) - 2*_sol(i,j+2);
               c = _sol(i-1,j+1)*_sol(i-1,j+1) + _sol(i,j+2)*_sol(i,j+2) - _hx*_hy;
               break;

            case 21:
               a = 2;
               b = -2*_sol(i+1,j+1) - 2*_sol(i,j);
               c = _sol(i+1,j+1)*_sol(i+1,j+1) + _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            case 22:
               a = 2;
               b = -2*_sol(i+1,j+1) - 2*_sol(i,j+2);
               c = _sol(i+1,j+1)*_sol(i+1,j+1) + _sol(i,j+2)*_sol(i,j+2) - _hx*_hy;
               break;

            case 31:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            case 32:
               a = 1;
               b = -2*_sol(i,j+2);
               c = _sol(i,j+2)*_sol(i,j+2) - _hx*_hy;
               break;

            default:
               ps = false;
               break;
         }
         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i,j+1) was not gray
         if (_sol(i,j+1)!=-2) {
//          Set to gray
            Add(ug,i,j+1);
            _sol(i,j+1) = -2;
         }
         else
//          If the new value is smaller than the previus one, modify it
            if (ug<_heap[_pos(i,j+1)].val)
               Modify(ug,i,j+1);
      }

//    If its south neighbor is not already black and lies in the domain
      if (_sol(i,j-1)<0 && i>=2 && i<=_nx && j>=3) {
//       The new origin: (i,j-1)
//       If the right and left neighbors are not black
         if (_sol(i+1,j-1)<0 && _sol(i-1,j-1)<0)
            max_x = 30;                              // max(D-x,-D+x,0)=0
//       If the right and left neighbors are black
         else if (_sol(i+1,j-1)>=0 && _sol(i-1,j-1)>=0) {
//          Choose the minimum
            if (_sol(i+1,j-1)>_sol(i-1,j-1))
               max_x = 10;                              // max(D-x,-D+x,0) = D-x
            else
               max_x = 20;                              // max(D-x,-D+x,0)=-D+x
         }

//       If the right neighbor is black and not the left one
         else if (_sol(i+1,j-1)>=0)
            max_x = 20;                                    // max(D-x,-D+x,0)=-D+x
//          If the left neighbor is black and not the right one
         else
            max_x = 10;                                    // max(D-x,-D+x,0)=D-x
         max_y = 2;                                        // max(D-y,-D+y,0)=-D+y

//       Unless the south neighbor is black and its value is smaller than the south one's
         if (_sol(i,j-2)>=0 && _sol(i,j)>_sol(i,j-2))
            max_y = 1;                                     // max(D-y,-D+y,0)=D-y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i-1,j-1) - 2*_sol(i,j-2);
               c = _sol(i-1,j-1)*_sol(i-1,j-1) + _sol(i,j-2)*_sol(i,j-2) - _hx*_hy;
               break;

            case 12:
               a = 2;
               b = -2*_sol(i-1,j-1) - 2*_sol(i,j);
               c = _sol(i-1,j-1)*_sol(i-1,j-1) + _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            case 21:
               a = 2;
               b = -2*_sol(i+1,j-1) - 2*_sol(i,j-2);
               c = _sol(i+1,j-1)*_sol(i+1,j-1) + _sol(i,j-2)*_sol(i,j-2) - _hx*_hy;
               break;

            case 22:
               a = 2;
               b = -2*_sol(i+1,j-1) - 2*_sol(i,j);
               c = _sol(i+1,j-1)*_sol(i+1,j-1) + _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            case 31:
               a = 1;
               b = -2*_sol(i,j-2);
               c = _sol(i,j-2)*_sol(i,j-2) - _hx*_hy;
               break;

            case 32:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) - _hx*_hy;
               break;

            default:
               ps = false;
               break;
         }
         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i,j-1) was not gray
         if (_sol(i,j-1)!=-2) {
//          Set to gray
            Add(ug,i,j-1);
            _sol(i,j-1) = -2;
         }

         else
//          If the new value is smaller than the previus one, modify it
            if (ug<_heap[_pos(i,j-1)].val)
               Modify(ug,i,j-1);
      }
   }
}


void FastMarching2D::setGrayWithObstacle(size_t i,
                                         size_t j)
{
   size_t max_x, max_y;
   real_t a, b, c, ug=0;
   bool ps;

   if (i>3 && i<_nx && j>3 && j<_ny) {
      if (_sol(i-1,j)<0 && i>=4 && i<=_nx && j>=3 && j<=_ny-1) {
//       New origin: (i-1,j)
         max_x = 20;                           // max(D-x,-D+x,0)=-D+x
//       Unless the left neighbor is black and its value is smaller than the right one's
         if (_sol(i-2,j)>=0 && _sol(i,j)>_sol(i-2,j))
            max_x = 10;                        // max(D-x,-D+x,0)=D-x

//       If north and south neighbors are not black
         if (_sol(i-1,j+1)<0 && _sol(i-1,j-1)<0) 
            max_y = 3;                         // max(D-y,-D+y,0)=0
//       If north and south neighbors are black
         else if (_sol(i-1,j+1)>=0 && _sol(i-1,j-1)>=0) {
//          Choose the min.
            if (_sol(i-1,j-1)>_sol(i-1,j+1))
               max_y = 2;                      // max(D-y,-D+y,0)=-D+y
            else
               max_y = 1;                      // max(D-y,-D+y,0)=D-y
         }
//       If the north neighbor is black and not the south one
         else if (_sol(i-1,j+1)>=0)
            max_y = 2;                       // max(D-y,-D+y,0)=-D+y
//          If the South neighbor is black and not the north one
         else
            max_y = 1;                       // max(D-y,-D+y,0)=D-y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i-2,j) - 2*_sol(i-1,j-1);
               c = _sol(i-2,j)*_sol(i-2,j) + _sol(i-1,j-1)*_sol(i-1,j-1) -
                   _hx*_hy*Obstacle(_g.getCoord(i-1,j))*Obstacle(_g.getCoord(i-1,j));
               break;

            case 12:
               a = 2;
               b = -2*_sol(i-2,j) - 2*_sol(i-1,j+1);
               c = _sol(i-2,j)*_sol(i-2,j) + _sol(i-1,j+1)*_sol(i-1,j+1) -
                   _hx*_hy*Obstacle(_g.getCoord(i-1,j))*Obstacle(_g.getCoord(i-1,j));
               break;

            case 13:
               a = 1;
               b = -2*_sol(i-2,j);
               c = _sol(i-2,j)*_sol(i-2,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i-1,j))*Obstacle(_g.getCoord(i-1,j));
               break;

            case 21:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i-1,j-1);
               c = _sol(i,j)*_sol(i,j) + _sol(i-1,j-1)*_sol(i-1,j-1) -
                   _hx*_hy*Obstacle(_g.getCoord(i-1,j))*Obstacle(_g.getCoord(i-1,j));
               break;

            case 22:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i-1,j+1);
               c = _sol(i,j)*_sol(i,j) + _sol(i-1,j+1)*_sol(i-1,j+1) -
                   _hx*_hy*Obstacle(_g.getCoord(i-1,j))*Obstacle(_g.getCoord(i-1,j));
               break;

            case 23:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i-1,j))*Obstacle(_g.getCoord(i-1,j));
               break;

            default:
               ps = false;
               break;
         }

         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i-1,j) was not gray
         if (_sol(i-1,j)!=-2) {
            Add(ug,i-1,j);
            _sol(i-1,j) = -2;
         }
         else
//          If the new value is smaller than the previus then modify it
            if (ug<_heap[_pos(i-1,j)].val)
               Modify(ug,i-1,j);
      }

//    If its right neighbor is not already black and lies in the domain
      if (_sol(i+1,j)<0 && i>=2 && i<=_nx-1 && j>=3 && j<=_ny-1) {
         max_x = 10;                                   // max(D-x,-D+x,0)=D-x
//       Unless the right neighbor is black and its value is smaller than the left one's
         if (_sol(i+2,j)>=0 && _sol(i,j)>_sol(i+2,j))
            max_x = 20;                                // max(D-x,-D+x,0)=-D+x
//       If the north and south neighbors are not black
         if (_sol(i+1,j+1)<0 && _sol(i+1,j-1)<0)
            max_y = 3;                                 // max(D-y,-D+y,0)=0
//       If the north and south neighbors are black
         else if (_sol(i+1,j+1)>=0 && _sol(i+1,j-1)>=0)
//          We choose the minimum
            if (_sol(i+1,j-1)>_sol(i+1,j+1))
               max_y = 2;                             // max(D-y,-D+y,0)=-D+y
            else
               max_y = 1;                             // max(D-y,-D+y,0)=D-y
//       If the north neighbor is black and not the south one
         else if (_sol(i+1,j+1)>=0)
            max_y = 2;                                // max(D-y,-D+y,0)=-D+y
//          If the south neighbor is black and not the north one
         else
            max_y = 1;                                // max(D-y,-D+y,0)=D-y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i+1,j-1);
               c = _sol(i,j)*_sol(i,j) + _sol(i+1,j-1)*_sol(i+1,j-1) -
                   _hx*_hy*Obstacle(_g.getCoord(i+1,j))*Obstacle(_g.getCoord(i+1,j));
               break;

            case 12:
               a = 2;
               b = -2*_sol(i,j) - 2*_sol(i+1,j+1);
               c = _sol(i,j)*_sol(i,j) + _sol(i+1,j+1)*_sol(i+1,j+1) -
                   _hx*_hy*Obstacle(_g.getCoord(i+1,j))*Obstacle(_g.getCoord(i+1,j));
               break;

            case 13:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i+1,j))*Obstacle(_g.getCoord(i+1,j));
               break;

            case 21:
               a = 2;
               b = -2*_sol(i+2,j) - 2*_sol(i+1,j-1);
               c = _sol(i+2,j)*_sol(i+2,j) + _sol(i+1,j-1)*_sol(i+1,j-1) -
                   _hx*_hy*Obstacle(_g.getCoord(i+1,j))*Obstacle(_g.getCoord(i+1,j));
               break;

            case 22:
               a = 2;
               b = -2*_sol(i+2,j) - 2*_sol(i+1,j+1);
               c = _sol(i+2,j)*_sol(i+2,j) + _sol(i+1,j+1)*_sol(i+1,j+1) -
                   _hx*_hy*Obstacle(_g.getCoord(i+1,j))*Obstacle(_g.getCoord(i+1,j));
               break;

            case 23:
               a = 1;
               b = -2*_sol(i+2,j);
               c = _sol(i+2,j)*_sol(i+2,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i+1,j))*Obstacle(_g.getCoord(i+1,j));
               break;

            default:
               ps = false;
               break;
         }
         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i+1,j) was not gray
         if (_sol(i+1,j)!=-2) {
            Add(ug,i+1,j);
            _sol(i+1,j) = -2;
         }
         else {
//          If the new value is smaller than the previus one, modify it
            if (ug<_heap[_pos(i+1,j)].val)
               Modify(ug,i+1,j);
         }
      }

//    If north neighbor is not already black and lies in the domain
      if (_sol(i,j+1)<0 && i>=3 && i<=_nx-1 && j>=2 && j<=_ny-1) {
//       New origin: (i,j+1) 
//       If right and left neighbors are not black
         if (_sol(i+1,j+1)<0 && _sol(i-1,j+1)<0)
            max_x = 30;                               // max(D-x,-D+x,0)=0
//       If right and left neighbors are black
         else if (_sol(i+1,j+1)>=0 && _sol(i-1,j+1)>=0) {
//          Choose the minimum
            if (_sol(i+1,j+1)>_sol(i-1,j+1))
               max_x = 10;                            // max(D-x,-D+x,0)=D-x
            else
               max_x = 20;                            // max(D-x,-D+x,0)=-D+x
         }

//       If right neighbor is black and not the left one
         else if (_sol(i+1,j+1)>=0)
            max_x = 20;                               // max(D-x,-D+x,0)=-D+x
//       If left neighbor is black and not the right one
         else
            max_x = 10;                               // max(D-x,-D+x,0)=D-x
         max_y = 1;                                   // max(D-y,-D+y,0)=D-y
//       Unless the north neighbor is black and its value is smaller than the south one's
         if (_sol(i,j+2)>=0 && _sol(i,j)>_sol(i,j+2))
            max_y = 2;                                // max(D-y,-D+y,0)=-D+y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i-1,j+1) - 2*_sol(i,j);
               c = _sol(i-1,j+1)*_sol(i-1,j+1) + _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j+1))*Obstacle(_g.getCoord(i,j+1));
               break;

            case 12:
               a = 2;
               b = -2*_sol(i-1,j+1) - 2*_sol(i,j+2);
               c = _sol(i-1,j+1)*_sol(i-1,j+1) + _sol(i,j+2)*_sol(i,j+2) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j+1))*Obstacle(_g.getCoord(i,j+1));
               break;

            case 21:
               a = 2;
               b = -2*_sol(i+1,j+1) - 2*_sol(i,j);
               c = _sol(i+1,j+1)*_sol(i+1,j+1) + _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j+1))*Obstacle(_g.getCoord(i,j+1));
               break;

            case 22:
               a = 2;
               b = -2*_sol(i+1,j+1) - 2*_sol(i,j+2);
               c = _sol(i+1,j+1)*_sol(i+1,j+1) + _sol(i,j+2)*_sol(i,j+2) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j+1))*Obstacle(_g.getCoord(i,j+1));
               break;

            case 31:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j+1))*Obstacle(_g.getCoord(i,j+1));
               break;

            case 32:
               a = 1;
               b = -2*_sol(i,j+2);
               c = _sol(i,j+2)*_sol(i,j+2) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j+1))*Obstacle(_g.getCoord(i,j+1));
               break;

            default:
               ps = false;
               break;
         }
         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i,j+1) was not gray
         if (_sol(i,j+1)!=-2) {
//          Set to gray
            Add(ug,i,j+1);
            _sol(i,j+1) = -2;
         }
         else
//          If the new value is smaller than the previus one, modify it
            if (ug<_heap[_pos(i,j+1)].val)
               Modify(ug,i,j+1);
      }

//    If its south neighbor is not already black and lies in the domain
      if (_sol(i,j-1)<0 && i>=3 && i<=_nx-1 && j>=3 && j<=_ny+1) {
//       The new origin: (i,j-1)
//       If the right and left neighbors are not black
         if (_sol(i+1,j-1)<0 && _sol(i-1,j-1)<0)
            max_x = 30;                              // max(D-x,-D+x,0)=0
//       If the right and left neighbors are black
         else if (_sol(i+1,j-1)>=0 && _sol(i-1,j-1)>=0) {
//          Choose the minimum
            if (_sol(i+1,j-1)>_sol(i-1,j-1))
               max_x = 10;                              // max(D-x,-D+x,0) = D-x
            else
               max_x = 20;                              // max(D-x,-D+x,0)=-D+x
         }

//       If the right neighbor is black and not the left one
         else if (_sol(i+1,j-1)>=0)
            max_x = 20;                                    // max(D-x,-D+x,0)=-D+x
//          If the left neighbor is black and not the right one
         else
            max_x = 10;                                    // max(D-x,-D+x,0)=D-x
         max_y = 2;                                        // max(D-y,-D+y,0)=-D+y

//       Unless the south neighbor is black and its value is smaller than the south one's
         if (_sol(i,j-2)>=0 && _sol(i,j)>_sol(i,j-2))
            max_y = 1;                                     // max(D-y,-D+y,0)=D-y

         ps = true;
         switch (max_x+max_y) {

            case 11:
               a = 2;
               b = -2*_sol(i-1,j-1) - 2*_sol(i,j-2);
               c = _sol(i-1,j-1)*_sol(i-1,j-1) + _sol(i,j-2)*_sol(i,j-2) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j-1))*Obstacle(_g.getCoord(i,j-1));
               break;

            case 12:
               a = 2;
               b = -2*_sol(i-1,j-1) - 2*_sol(i,j);
               c = _sol(i-1,j-1)*_sol(i-1,j-1) + _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j-1))*Obstacle(_g.getCoord(i,j-1));
               break;

            case 21:
               a = 2;
               b = -2*_sol(i+1,j-1) - 2*_sol(i,j-2);
               c = _sol(i+1,j-1)*_sol(i+1,j-1) + _sol(i,j-2)*_sol(i,j-2) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j-1))*Obstacle(_g.getCoord(i,j-1));
               break;

            case 22:
               a = 2;
               b = -2*_sol(i+1,j-1) - 2*_sol(i,j);
               c = _sol(i+1,j-1)*_sol(i+1,j-1) + _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j-1))*Obstacle(_g.getCoord(i,j-1));
               break;

            case 31:
               a = 1;
               b = -2*_sol(i,j-2);
               c = _sol(i,j-2)*_sol(i,j-2) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j-1))*Obstacle(_g.getCoord(i,j-1));
               break;

            case 32:
               a = 1;
               b = -2*_sol(i,j);
               c = _sol(i,j)*_sol(i,j) -
                   _hx*_hy*Obstacle(_g.getCoord(i,j-1))*Obstacle(_g.getCoord(i,j-1));
               break;

            default:
               ps = false;
               break;
         }
         if (ps)
            ug = PositiveSol(a,b,c);

//       If (i,j-1) was not gray
         if (_sol(i,j-1)!=-2) {
//          Set to gray
            Add(ug,i,j-1);
            _sol(i,j-1) = -2;
         }

         else
//          If the new value is smaller than the previus one, modify it
            if (ug<_heap[_pos(i,j-1)].val)
               Modify(ug,i,j-1);
      }
   }
}


void FastMarching2D::ExtendLocalVelocity(Vect<real_t>& dis)
{
   if (_nx<5)
      throw OFELIException("FastMarching2D::ExtendLocalVelocity(Vect<real_t>): "
                           "Number of grid intervals in x must be larger than 4");
   if (_ny<5)
      throw OFELIException("FastMarching2D::ExtendLocalVelocity(Vect<real_t>): "
                           "Number of grid intervals in y must be larger than 4");

   Point2D<real_t> x1, x2;
   Vect<real_t> bd(_f), a(_f);
   size_t i, j;

   for (i=1; i<=_nx+1; i++) {
      for (j=1; j<=_ny+1; j++)
         if (a(i,j)!=0)
            a(i,j) = 1;
         else
            _f(i,j) = _u(i,j);
   }

   for (i=2; i<=_nx; i++) {
      for (j=2; j<=_ny; j++) {
         if ((bd(i,j)*bd(i+1,j)<=0) || (bd(i,j)*bd(i-1,j)<=0) ||
             (bd(i,j)*bd(i,j+1)<=0) || (bd(i,j)*bd(i,j-1)<=0)) {

            if (bd(i,j)*bd(i+1,j)<=0) {
               x1 = _g.getCoord(i,j);
               if (bd(i,j)!=0)
                  x1.x += (_g.getX(i+1)-_g.getX(i))*bd(i,j)/(bd(i,j)-bd(i+1,j));
               if (bd(i,j)*bd(i,j-1)<=0) {         // X--X
                  x2 = _g.getXY(i,j);              // |
                  if (bd(i,j)!=0)                  // x
                     x2.y += (_g.getY(j-1)-_g.getY(j))*bd(i,j)/(bd(i,j)-bd(i,j-1));
                  Int(i,j,a,x1,x2);
               }

               if (bd(i+1,j-1)*bd(i,j-1)<=0) {
                  x2.y = _g.getY(j-1);             // X--X
                  if (bd(i+1,j-1)==0)              //
                     x2.x = _g.getX(i+1);          // x--x
                  else
                     x2.x = _g.getX(i) + (_g.getX(i+1)-_g.getX(i))*bd(i,j-1)/(bd(i,j-1)-bd(i+1,j-1));
                  Int(i,j,a,x1,x2);
               }

               if (bd(i+1,j-1)*bd(i+1,j)<=0) {     // X--X
                  x2 = _g.getXY(i+1,j);            //    |
                  if (bd(i+1,j)!=0)                //    x
                     x2.y += (_g.getY(j-1)-_g.getY(j))*bd(i+1,j)/(bd(i+1,j)-bd(i+1,j-1));
                  Int(i,j,a,x1,x2);
               }
            }

            if (bd(i,j)*bd(i,j+1)<=0) {
               x1.x = _g.getX(i);
               if (bd(i,j+1)==0)
                  x1.y = _g.getY(j+1);
               else
                  x1.y = _g.getY(j) + (_g.getY(j+1)-_g.getY(j))*bd(i,j)/(bd(i,j)-bd(i,j+1));

               if (bd(i-1,j)*bd(i,j)<=0) {         //    X
                  x2 = _g.getXY(i,j);              //    |
                  if (bd(i,j)!=0)                  // x--X
                     x2.x += (_g.getX(i-1)-_g.getX(i))*bd(i,j)/(bd(i,j)-bd(i-1,j));
                  Int(i,j,a,x1,x2);
               }

               if (bd(i-1,j)*bd(i-1,j+1)<=0) {     // x  X
                  x2 = _g.getXY(i-1,j);            // |  |
                  if (bd(i-1,j)!=0)                // x--x
                     x2.y += (_g.getY(j+1)-_g.getY(j))*bd(i-1,j)/(bd(i-1,j)-bd(i-1,j+1));
                  Int(i,j,a,x1,x2);
               }

               if (bd(i,j)*bd(i+1,j)<=0) {         // X
                  x2 = _g.getXY(i,j);              // |
                  if (bd(i,j)!=0)                  // X--x
                     x2.x += (_g.getX(j+1)-_g.getX(j))*bd(i,j)/(bd(i,j)-bd(i+1,j));
                  Int(i,j,a,x1,x2);
               }
            }
         }
      }
   }
   for (i=1; i<=_nx+1; i++) {
      for (j=1; j<=_ny+1; j++) {
         if (a(i,j)>_bw*_hx)
            a(i,j) = 1;
      }
   }
   dis = a;
}


void FastMarching2D::Int(size_t           k,
                         size_t           l,
                         Vect<real_t>&    a,
                         Point2D<real_t>& x1,
                         Point2D<real_t>& x2)
{
   for (int bbi=-_bw; bbi<=_bw; bbi++) {
      for (int bbj=-_bw; bbj<=_bw; bbj++) {
         int bi = bbi + k;
         int bj = bbj + l;
         if (bi>0 && bj>0 && bi<=int(_nx) && bj<=int(_ny)) {
            Point2D<real_t> xx = _g.getCoord(bi,bj);
            if (a(bi,bj)!=0) {
               if (a(bi,bj)>Dist(xx,x1,x2)) {
                  a(bi,bj) = std::min(a(bi,bj),Dist(xx,x1,x2));
                  real_t alpha = 0;
                  if (SqrDistance(x1,x2)>0)
                     alpha = ((xx-x1)*(x2-x1))/SqrDistance(x1,x2);
                  if (alpha<0)
                     alpha = 0;
                  if (alpha>1)
                     alpha = 1;
                  _f(bi,bj) = alpha*Interp(x2) + (1-alpha)*Interp(x1);
               }
            }
         }
      }
   }
}


real_t FastMarching2D::Dist(const Point2D<real_t>& x,
                            const Point2D<real_t>& y,
                            const Point2D<real_t>& z)
{
   real_t dist;
   if (y==z)
      dist = Distance(x,y);
   else {
      if ((x-y)*(z-y)>0 && (x-z)*(y-z)>0)
         dist = std::abs((y.x-x.x)*(z.y-x.y)-(y.y-x.y)*(z.x-x.x))/Distance(z,y);
      else
         dist = sqrt(std::min(SqrDistance(x,y),SqrDistance(x,z)));
   }
   return dist;
}


real_t FastMarching2D::Interp(Point2D<real_t>& xx)
{
   real_t ui=0, alpha=0;
   int ifi=0, jfi=0;
   if (xx.x>=_g.getX(_nx))
      throw OFELIException("FastMarching2D::Interp(Point2D<real_t>): "
                           "Interpolation points goes beyond extremal point.");
   if (xx.y>=_g.getY(_ny))
      throw OFELIException("FastMarching2D::Interp(Point2D<real_t>): "
                           "Interpolation points goes beyond extremal point.");

   size_t i=1, j=1;
   while (_g.getX(i)<xx.x) 
      i++;
   i--;
   while (_g.getY(j)<xx.y)
      j++;
   j--;
   if (i==0)
      throw OFELIException("FastMarching2D::Interp(Point2D<real_t>): Illegal index i=0");
   if (j==0)
      throw OFELIException("FastMarching2D::Interp(Point2D<real_t>): Illegal index j=0");

   if (xx.x-_g.getX(i)<4*DBL_EPSILON) {
      ifi = i;
      jfi = 0;
      xx.x = _g.getX(i);
   }
   if (_g.getX(i+1)-xx.x<4*DBL_EPSILON) {
      ifi = i + 1;
      jfi = 0;
      xx.x = _g.getX(i+1);
   }
   if (xx.y-_g.getY(j)<4*DBL_EPSILON) {
      jfi = j;
      ifi = 0;
      xx.y = _g.getY(j);
   }
   if (_g.getY(j+1)-xx.y<4*DBL_EPSILON) {
      jfi = j + 1; 
      ifi = 0;
      xx.y = _g.getY(j+1);
   }
   if (ifi==0) {
      alpha = (xx.x-_g.getX(i))/(_g.getX(i+1)-_g.getX(i));
      ui = alpha*_u(i+1,jfi) + (1-alpha)*_u(i,jfi);
   }
   else if (jfi==0) {
      alpha = (xx.y-_g.getY(j))/(_g.getY(j+1)-_g.getY(j));
      ui = alpha*_u(ifi,j+1) + (1-alpha)*_u(ifi,j);
   }
   else
      throw OFELIException("FastMarching2D::Interp(Point2D<real_t>): Illegal index.");
   return ui;
}


real_t FastMarching2D::Obstacle(Point2D<real_t> a)
{
   real_t ret = 1.0;
   a = 0;

// Cas d'un obstacle rond au centre
/*    !if((a-0.5)*(a-0.5)+(b-0.5)*(b-0.5).lt.0.0625) f=1000 
   
    !       Cas de plusieurs obstacles de types quart de disques 
    !if(((a-0.20)*(a-0.20)+(b-0.5)*(b-0.5).le.0.0625).and.((a-0.20)*(a-0.20)+(b-0.5)*(b-0.5).ge.0.04).and.(a.ge.0.20).and.(b.ge.0.5)) f=1000
    ! if(((a-0.5)*(a-0.5)+(b-0.5)*(b-0.5).le.0.062).and.((a-0.5)*(a-0.5)+(b-0.5)*(b-0.5).ge.0.04).and.(a.le.0.5).and.(b.le.0.5)) f=1000
    ! if(((a-0.45)*(a-0.45)+(b-0.6)*(b-0.6).le.0.062).and.((a-0.45)*(a-0.45)+(b-0.6)*(b-0.6).ge.0.04).and.(a.ge.0.45).and.(b.ge.0.6)) f=1000
    ! if(((a-0.8)*(a-0.8)+(b-0.5)*(b-0.5).le.0.062).and.((a-0.8)*(a-0.8)+(b-0.5)*(b-0.5).ge.0.04).and.(a.le.0.8).and.(b.le.0.5)) f=1000    
*/
   return ret;
}

} /* namespace OFELI */
