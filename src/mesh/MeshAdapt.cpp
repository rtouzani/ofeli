/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

                       Implementation of class 'MeshAdapt'

  ==============================================================================*/

#include "mesh/MeshAdapt.h"
#include "mesh/getMesh.h"


namespace OFELI {


MeshAdapt::MeshAdapt()
{
   setDefault();
}

   
MeshAdapt::MeshAdapt(Domain& dom)
{
   set(dom);
}


MeshAdapt::MeshAdapt(Mesh& ms)
{
   setDefault();
   set(ms);
}


void MeshAdapt::setDefault()
{
   _nb_Jacobi = 3;
   _nb_smooth = 3;
   _abs_error = true;
   _err = 0.004;
   _geo_err = 0.1;
   _verb = 0;
   _omega = 1.8;
   _hmin = 1.e-100;
   _hmax = 1.e17;
   _hmin_aniso = 1e-100;
   _aniso_max = 1.e06;
   _anisotropic = false;
   _max_nb_vertices = 500000;
   _ratio = 0;
   _scaling = true;
   _max_subdiv = 10;
   _power = 1;
   _keep_vertices = true;
   _cut_off = 1.e-5;
   _coef = 1;
   _hessian = false;
   _set_outm = _set_geo = _set_bgm = false;
   _splitbedge = false;
   _set_metric = false;
   _cutoff_radian = -1;
   _set_mbb = _set_mBB = false;
   _set_rbb = _set_rBB = _set_wbb = _set_wBB = false;
   _set_meshr = _set_ometric = false;
   _allquad = 0;
   _iter = 0;
   _domain = NULL;
   _nb_subdiv = 300;
}


void MeshAdapt::set(Domain& dom)
{
   setDefault();
   _domain = &dom;
   _scale_fact = max(_domain->getMaxCoord().x,_domain->getMaxCoord().y);
   *_domain *= 1./_scale_fact;
   setGeoFile("mesh.geo");
   _domain->genGeo(_geo_file);
   setOutputMesh("mesh.0.bamg");
   run();
}


void MeshAdapt::set(Mesh& ms)
{
   _ms.push_back(&ms);
   _scale_fact = max(_domain->getMaxCoord().x,_domain->getMaxCoord().y);
   *_ms[0] *= 1./_scale_fact;
   _nb_nodes = _ms[0]->getNbNodes();
   _nb_elements = _ms[0]->getNbElements();
   setOutputMesh("mesh.0.bamg");
   saveBamg(_output_mesh_file,ms);
   _iter = 1;
}


void MeshAdapt::saveMbb(string              file,
                        const Vect<real_t>& u)
{
   const size_t nb_sol=1;
   ofstream fbb(file.c_str());
   fbb << "  2  " << nb_sol << "  " << u.size() << "  2" << endl;
   for (size_t i=0; i<u.size(); i++) {
      for (size_t n=1; n<=nb_sol; n++)
         fbb << u[i] << endl;
   }
}


void MeshAdapt::setSolution(const Vect<real_t>& u)
{
   _u = &u;
   _mbb_file = _rbb_file = "adapt.rbb";
   _set_mbb = _set_rbb = true;
   const size_t nb_sol=1;
   ofstream fbb("adapt.rbb");
   fbb << "  2  " << nb_sol << "  " << _u->size() << "  2" << endl;
   for (size_t i=0; i<_u->size(); i++) {
      for (size_t n=1; n<=nb_sol; n++)
         fbb << (*_u)[i] << endl;
   }
}


int MeshAdapt::run(const Vect<real_t>& u)
{
   saveMbb("MeshAdapt.bb",u);
   getSolutionMbb("MeshAdapt.bb");
   int ret = run();
   remove("MeshAdapt.bb");
   return ret;
}


int MeshAdapt::run(const Vect<real_t>& u,
                   Vect<real_t>&       v)
{
   saveMbb("MeshAdapt.bb",u);
   getSolutionMbb("MeshAdapt.bb");
   int ret = run();
   v.setMesh(_theMesh);
   MeshToMesh(u,v,_nb_subdiv,_nb_subdiv,0);
   remove("MeshAdapt.bb");
   return ret;
}


int MeshAdapt::run()
{
   verbosity = _verb;
   MeshIstreamErrorHandler = MeshErrorIO;
   hinterpole = 1;
   int fileout=0, NoMeshReconstruction=0;
   real_t costheta=2;
   real_t *solMbb=0, *solMBB=0;
   int *typesolsBB=NULL;
   long nbsolbb=0, lsolbb=0, nbsolBB=0, lsolBB=0;
   int rbbeqMbb=0, rBBeqMBB=0;
   Triangles *Thr=NULL, *Thb=NULL;

   if (_iter>0) {
      setBackgroundMesh(_output_mesh_file);
      setOutputMesh("mesh."+itos(_iter)+".bamg");
      _set_geo = false;
   }

   NoMeshReconstruction = _set_meshr;
   if (!_set_bgm)
      _background_mesh_file = _meshr_file;
   fileout = _set_outm || _set_wbb || _set_wBB;
   if (!fileout && !_set_ometric) {
      cerr << "Error: No output file given " << endl;
      exit(1);
   }
   if (_max_subdiv>10 || _max_subdiv<=1.0) {
      cerr << "maxsubdiv " << _max_subdiv << " is not in (1,100]" << endl;
      exit(3);
   }
   if (!_anisotropic)
      _aniso_max = 1;
   if (!_set_metric && !_set_mbb)
      _nb_smooth = 0;
   if (!_set_rbb)
      _rbb_file = _bb_file;
   if (!_set_rBB)
      _rBB_file = _BB_file;
   if (_set_mbb && _set_rbb)
      rbbeqMbb = _rbb_file==_mbb_file;
   if (_set_mBB && _set_rBB)
      rBBeqMBB = _rBB_file==_mBB_file;
   if (_verb) {
      if (_set_geo && fileout && _verb)
         cout << "Construction of mesh from the geometry file " << _geo_file << endl;
      else if (_set_bgm && fileout && _verb) {
         if (NoMeshReconstruction)
            cout << "Modification of adapted mesh file " << _background_mesh_file << endl;
         else
            cout << "Construction of adapted mesh from background mesh file "
                 << _background_mesh_file << endl;
      }
      else if (_set_bgm && _set_ometric && _verb)
         cout << "Construction of metric file on background mesh "
              << _background_mesh_file << endl;
   }

   if (_set_geo && fileout) {
      if (_verb)
         cout << "Construction of mesh from the geometry file " << _geo_file << endl;
      Geometry Gh(_geo_file.c_str());
      _hmin = max(_hmin,Gh.MinimalHmin());
      _hmax = min(_hmax,Gh.MaximalHmax());
      if (_set_metric)
         Gh.ReadMetric(_metric_file.c_str(),_hmin,_hmax,_coef);
      else {
         for (long iv=0; iv<Gh.nbv; iv++) {
            MetricAnIso M = Gh[iv];
            MatVVP2x2 Vp(M/_coef);
            Vp.Maxh(_hmin);
            Vp.Minh(_hmax);
            Gh._vertices[iv].m = Vp;
         }
      }
      Triangles Th(_max_nb_vertices,Gh);
      if (costheta<=1)
         Th.MakeQuadrangles(costheta);
      if (_allquad)
         Th.SplitElement(_allquad==2);
      if (_splitbedge)
         Th.SplitInternalEdgeWithBorderVertices();
      Th.ReNumberingTheTriangleBySubDomain();
      if (_verb>3)
         Th.ShowHistogram();
      if (_nb_smooth>0)
         Th.SmoothingVertex(_nb_smooth,_omega);
      if (_verb>3 && _nb_smooth>0)
         Th.ShowHistogram();
      if (_set_outm)
         Th.Write(_output_mesh_file.c_str(),Triangles::BDMesh);
      if (_verb>0) {
         cout << "Number of vertices: " << Th.nbv;
         if (Th.nbt-Th.NbOutT-Th.NbOfQuad*2) 
            cout << "Number of Triangles: " << (Th.nbt-Th.NbOutT-Th.NbOfQuad*2);
         cout << endl;
      }
   }

   else if (_set_bgm &&
            (_set_metric || _set_mbb || _set_mBB || NoMeshReconstruction) &&
            (fileout || _set_ometric || (_set_rbb&&_set_wbb) || (_set_rBB&&_set_wBB))) {
      Triangles BTh(_background_mesh_file.c_str(),_cutoff_radian);
      _hmin = max(_hmin,BTh.MinimalHmin());
      _hmax = min(_hmax,BTh.MaximalHmax());
      BTh.MakeQuadTree();
      if (_set_metric)
         BTh.ReadMetric(_metric_file.c_str(),_hmin,_hmax,_coef);
      else {
         Metric Mhmax(_hmax);
         for (long iv=0; iv<BTh.nbv; iv++)
            BTh[iv].m = Mhmax;
      }
      if (_set_mbb) {
         solMbb = ReadbbFile(_mbb_file.c_str(),nbsolbb,lsolbb,2,2);
         if (lsolbb != BTh.nbv) {
            cerr << "fatal error nbsol " << nbsolbb << " " << lsolbb << " =! "
		 << BTh.nbv << endl;
            cerr << "size of sol incorrect" << endl;
            MeshError(99);
         }
         assert(lsolbb==BTh.nbv);
         BTh.IntersectConsMetric(solMbb,nbsolbb,0,_hmin,_hmax,sqrt(_err)*_coef,1e30,
                                 _abs_error?0.0:_cut_off,_nb_Jacobi,_scaling,_power,_hessian);
         if (!rbbeqMbb)
            delete [] solMbb, solMbb = NULL;
      }
      if (_set_mBB) {
         solMBB = ReadBBFile(_mbb_file.c_str(),nbsolBB,lsolBB,typesolsBB,2,2);
         if (lsolBB != BTh.nbv) {
            cerr << "fatal error nbsol " << nbsolBB << " " << lsolBB << " =! " << BTh.nbv << endl;
            cerr << "size of sol incorrect" << endl;
            MeshError(99);
         }
         assert(lsolBB==BTh.nbv);
         BTh.IntersectConsMetric(solMBB,nbsolBB,0,_hmin,_hmax,sqrt(_err)*_coef,1e30,
                                 _abs_error?0.0:_cut_off,_nb_Jacobi,_scaling,_hessian);
         if (!rBBeqMBB)
            delete [] solMBB, solMBB = NULL;
      }
      BTh.IntersectGeomMetric(_geo_err,!_anisotropic);
      if (_ratio)
         BTh.SmoothMetric(_ratio);
      BTh.MaxSubDivision(_max_subdiv);

      if (!_anisotropic)
         _aniso_max = 1;
      BTh.BoundAnisotropy(_aniso_max,_hmin_aniso);
      if (_set_ometric) {
         ofstream f(_ometric_file.c_str());
         if (f)
            BTh.WriteMetric(f,!_anisotropic);
      }
      
      if (fileout) {
         Triangles &Th( *(NoMeshReconstruction
                        ? new Triangles(*Thr,&Thr->Gh,Thb,_max_nb_vertices)
                        : new Triangles(_max_nb_vertices,BTh,_keep_vertices)));
         if (Thr != &BTh)
            delete Thr;
         if (costheta<=1.0)
            Th.MakeQuadrangles(costheta);
         if (_allquad)
            Th.SplitElement(_allquad==2);
         if (_splitbedge)
            Th.SplitInternalEdgeWithBorderVertices();
         Th.ReNumberingTheTriangleBySubDomain();
         if (_verb>3)
            Th.ShowHistogram();
         if (_nb_smooth>0)
            Th.SmoothingVertex(_nb_smooth,_omega);
         if (_verb>3 && _nb_smooth>0)
            Th.ShowHistogram();
         Th.BTh.ReMakeTriangleContainingTheVertex();
         if (_set_outm)
            Th.Write(_output_mesh_file.c_str(),Triangles::BDMesh);
         if ((_set_rbb&&_set_wbb) || (_set_rBB&&_set_wBB)) {
            if (_verb>1) {
               if (_set_rbb)
                  cout << "Interpolation P1 files: " << _rbb_file << ", " << _wbb_file << endl;
               if (_set_rBB)
                  cout << "Interpolation P1 files: " << _rBB_file << ", " << _wBB_file << endl;
            }
            const int dim = 2;
            real_t *solbb=NULL, *solBB=NULL;
            if (_set_rbb)
               solbb = rbbeqMbb ? solMbb : ReadbbFile(_rbb_file.c_str(),nbsolbb,lsolbb,2,2);
            if (_set_rBB)
               solBB = rBBeqMBB ? solMBB : ReadBBFile(_rBB_file.c_str(),nbsolBB,lsolBB,typesolsBB,2,2);
            if (!solBB && !solbb) {
               cerr << "Fatal Error " << "solBB = " << solBB << " solbb = " << solbb << endl;
               exit(2);
            }
            ofstream *fbb = _set_wbb ? new ofstream(_wbb_file.c_str()) : NULL;
            ofstream *fBB = _set_wBB ? new ofstream(_wBB_file.c_str()) : NULL;
            long nbfieldBB=0, nbfieldbb=nbsolbb;
            if (fbb)
               *fbb << dim << " " << nbsolbb << " " << Th.nbv << " " << 2 << endl;
            if (fBB) {
               *fBB << dim << " " << nbsolBB;
               for (int i=0; i<nbsolBB; i++)
                  *fBB << " " << (typesolsBB ? typesolsBB[i]+1:1);
               *fBB << " " << Th.nbv << " " << 2 << endl;
               assert(dim==2);
               for (int i=0; i<nbsolBB; i++)
                  nbfieldBB += typesolsBB ? typesolsBB[i]+1 : 1;
            }
            if (_verb)
               cout << "nb of fields BB: " << nbfieldBB << endl;
            for (int i=0; i<Th.nbv; i++) {
               long i0, i1, i2;
               real_t a[3];
               Icoor2 dete[3];
               I2 I = Th.BTh.toI2(Th._vertices[i].r);
               Triangle &tb = *Th.BTh.FindTriangleContaining(I,dete);
               if (tb.det>0) {
                  a[0] = dete[0]/tb.det;
                  a[1] = dete[1]/tb.det;
                  a[2] = dete[2]/tb.det;
                  i0 = Th.BTh.Number(tb[0]);
                  i1 = Th.BTh.Number(tb[1]);
                  i2 = Th.BTh.Number(tb[2]);
               }
               else {
                  real_t aa, bb;
                  TriangleAdjacent ta = CloseBoundaryEdge(I,&tb,aa,bb).Adj();
                  int k = ta;
                  Triangle &tc = *(Triangle *)ta;
                  i0 = Th.BTh.Number(tc[0]);
                  i1 = Th.BTh.Number(tc[1]);
                  i2 = Th.BTh.Number(tc[2]);
                  a[VerticesOfTriangularEdge[k][1]] = aa;
                  a[VerticesOfTriangularEdge[k][0]] = bb;
                  a[OppositeVertex[k]] = 1 - aa - bb;
               }
               long ibb0=nbfieldbb*i0, ibb1=nbfieldbb*i1, ibb2=nbfieldbb*i2;
               long iBB0=nbfieldBB*i0, iBB1=nbfieldBB*i1, iBB2=nbfieldBB*i2;
               for (int j=0; j<nbfieldbb; j++)
                  *fbb << " " << (a[0]*solbb[ibb0++] + a[1]*solbb[ibb1++] + a[2]*solbb[ibb2++]);
               for (int j=0; j<nbfieldBB; j++)
                  *fBB << " " << (a[0]*solBB[iBB0++] + a[1]*solBB[iBB1++] + a[2]*solBB[iBB2++]);
               if (fbb)
                  *fbb << endl;
               if (fBB)
                  *fBB << endl;
            }
         }
         delete &Th;
      }
   }
   _ms.push_back(new Mesh);
   getBamg(_output_mesh_file,*_ms[_iter]);
   _theMesh = *_ms[_iter];
   _theMesh *= _scale_fact;
   _nb_nodes = _ms[_iter]->getNbNodes();
   _nb_elements = _ms[_iter]->getNbElements();
   _iter++;
   if (_iter>1) {
      for (size_t i=0; i<=_iter; i++)
          remove(string("mesh."+itos(i)+".bamg").c_str());
      remove("mesh.geo");
   }
   return 0;
}


void MeshAdapt::getSolution(Vect<real_t>& u,
                            int           is)
{
   ifstream wbbf(_wbb_file.c_str());
   int k, l, n, d;
   wbbf >> k >> l >> n >> d;
   real_t uu;
   u.setMesh(getMesh());
   for (int i=0; i<n; i++) {
      for (int s=1; s<=l; s++) {
         wbbf >> uu;
         if (is==s)
            u[i] = uu;
      }
   }
}


void MeshAdapt::Interpolate(const Vect<real_t>& u,
                            Vect<real_t>&       v)
{
   v.setMesh(_theMesh);
   MeshToMesh(u,v,_nb_subdiv,_nb_subdiv,0);
}


ostream& operator<<(ostream&         s,
                    const MeshAdapt& a)
{
   cout << "\nSUMMARY OF ADAPTION PROCESS:" << endl;
   cout << "Number of adaption iterations:\t" << setw(6) << a._iter << endl;
   cout << "Number of nodes:\t\t" << setw(6) << a._nb_nodes << endl;
   cout << "Number of elements:\t\t" << setw(6) << a._nb_elements << endl;
   return s;
}

} /* namespace OFELI */
