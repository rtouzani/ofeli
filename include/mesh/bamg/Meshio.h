// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
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

#ifndef __MESHIO_H
#define __MESHIO_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cctype>

using namespace std;
//  PB compilo HP aCC 
#if defined(__hpux) || defined(__SUNPRO_CC) 
#define IOS_OPENMODE int
#else
#define IOS_OPENMODE ios::openmode
#endif
extern  long verbosity;

namespace bamg {

extern void(*MeshIstreamErrorHandler)(ios &);

void WriteStr(ostream& out,
              char*    str);

double *ReadbbFile(const char* file,
                         long& nbsol,
                         long& lsol,
                   const int   dim=2,
                   const int   typesol=2);

double *ReadBBFile(const char* file,
                         long& nbsol,
                         long& lsol,
                         int*& typesols,
                   const int   dim=2,
                   const int typesol=2);

// Solution at vertex (P1)

union Char4orLong { char c[4]; long l; };

class MeshIstream {

 public:
   istream& in;
   int LineNumber, LineError, opened;
   const char* CurrentFile;

   istream& cm() // Eat white space and comments 
   {
      char c;
      int cmm=0;
      while (in.get(c) &&  
         (isspace(c) ? (((c=='\n'||c==char(12)||c==char(15)) && (LineNumber++,cmm=0)),1) : (cmm || (c=='#' && (cmm=1) ))))
         ;
      if (in.good())
         in.putback(c);
      return in;
  }

   void eol()
   {
      char c;
      while (in.get(c) && (c!='\n') && (c!='\r'))
         ;
   }
  
   void ShowIoErr(int);

   MeshIstream& err()
   {
      if (!in.good() ) ShowIoErr(in.rdstate());
         return *this;
   }

   MeshIstream(const char* file_name)
    : in(*new ifstream(file_name)), LineNumber(1), LineError(0) 
   {
      CurrentFile = file_name;
      if (!in) {
         cerr << " Error Opening file " << file_name, CurrentFile=0;
         ShowIoErr(1);
      }
      if (verbosity>4)
         cout << "    Openfile : " << file_name << endl;
         err();
   }

   int eof()
   {
      return in.eof();
   }

   ~MeshIstream()
   {
      delete &in;
   }

   int IsString(const char* s);
   char * ReadStr();
   MeshIstream& operator>>(short& i) { cm() >> i; return err(); }
   MeshIstream& operator>>(long& i) { cm() >> i; return err(); }
   MeshIstream& operator>>(int& i) { cm() >> i; return err(); }
   MeshIstream& operator>>(float& i) { cm() >> i; return err(); }
   MeshIstream& operator>>(double& i) { cm() >> i; return err(); }
   MeshIstream& operator>>(char*& i) { i=ReadStr(); return err(); }
};
// Fortran unformatted file  interface ----------

class IFortranUnFormattedFile
{
   long i, l, n, j, nb_rec;
   istream* f;
   int to_close;
   const char *file_name;
 
 public:

   IFortranUnFormattedFile(char *name)
    : i(0), l(0), n((long)-int(sizeof(long))),
      nb_rec(0), to_close(1)
   {
      f = new ifstream(name);
      if (!*f)
         Error(0);
      file_name = name;
   }

   IFortranUnFormattedFile(MeshIstream & ff)
      : i(0), l(0), n((long)-int(sizeof(long))), nb_rec(0), to_close(0)
   {
      f = &ff.in;
      if (! *f)
         Error(0);
      file_name = ff.CurrentFile;
   }
  
 ~IFortranUnFormattedFile();
  long Record();
  long where(){return j-i;}
  void read4(char *c, int);// for the fortran 77 char4
  void read(char *p, const size_t lg);
  void Error(int);
};


class OFortranUnFormattedFile
{
  long i, l, n, j, nb_rec;
  ostream* f;
  const static char *unknown;
  int to_close;
  const char* file_name;
 
 public:
  
  OFortranUnFormattedFile(const char*        name,
                                IOS_OPENMODE mm=ios::trunc)
    : i(0), l(0), n((long) -int(sizeof(long))), nb_rec(0), to_close(1)
    {
       f = new ofstream(name,mm);
       if (!*f)
          Error(0);
       file_name = name;
    }

  OFortranUnFormattedFile(ostream &ff)
    : i(0), l(0), n((long) -int(sizeof(long))), nb_rec(0), to_close(0)
    {
       f = &ff; 
       if (!*f)
          Error(0);
       file_name = unknown;
    }
  
  ~OFortranUnFormattedFile();

  long Record(long ll=0);
  long where() { return j-i; }
  void write4(const char *c, int);// for the fortran 77 char4
  void write(const char *p, const size_t lg);
  void Error(int );
};


inline void  IFortranUnFormattedFile::read(      char*  p,
                                           const size_t lg)
{  
   f->read(p,lg);
   j += long(lg);
   if (j>n)
      Error(1);
   else if (!f->good())
      Error(2);
}


inline void OFortranUnFormattedFile::write(const char*  p,
                                           const size_t lg)
{
   f->write(p,lg);
   j += long(lg);
   if (l && j>n)
      Error(1);
   else if (!f->good())
      Error(2);
}

template<class T> inline
IFortranUnFormattedFile & operator>>(IFortranUnFormattedFile& f,
                                     T&                       l)
{
   f.read((char *) &l,sizeof(l));
   return f;
}

inline 
OFortranUnFormattedFile & operator<<(      OFortranUnFormattedFile& f,
                                     const int&                     l)
{
   f.write((char *) &l,sizeof(l));return f;
}

inline 
OFortranUnFormattedFile & operator<<(      OFortranUnFormattedFile& f,
                                     const long&                    l)
{
   f.write((char *) &l,sizeof(l));
   return f;
}

inline 
OFortranUnFormattedFile& operator<<(      OFortranUnFormattedFile& f,
                                    const double&                  l)
{
   f.write((char *) &l,sizeof(l));
   return f;
}

inline 
OFortranUnFormattedFile & operator<<(      OFortranUnFormattedFile& f,
                                     const float&                   l)
{
   f.write((char *) &l,sizeof(l));
   return f;
}

inline void OFortranUnFormattedFile::write4(const char* c,
                                                  int   ll)
{
   int i, j;
   Char4orLong ch4;
   for (i=0; i<ll; i++) {
      ch4.l = 0;
      for (j=0; j<4; j++)
         ch4.c[j] = *c ? *c++ : ' ';
      *this << ch4.l;
   }
}

inline void IFortranUnFormattedFile::read4(char* c,
                                           int   ll)
{
   int i, j;
   Char4orLong ch4;
   for (i=0; i<ll; i++) {
      *this >> ch4.l;
      for (j=0; j<4; j++)
	 *c++ = ch4.c[j];
   }
   *c = 0;
}

}

#endif

