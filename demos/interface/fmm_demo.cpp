#include "OFELI.h"
using namespace OFELI;

int main(int argc, char* argv[])
{
   void initial(int option, Grid *g, Vect<double> &u);
   const int nx=60, ny=60, nz=60;
   int option = 0;
   if (argc>1) 
      option = atoi(argv[1]);
   else {
      cout << "Choice of demo examples\n";
      cout << "1. 2-D example of 1 point obstacle\n";
      cout << "2. 2-D example of a rectangular interface\n";
      cout << "3. 3-D example of a spherical interface\n";
      cout << "0. Exit" << endl;
      cout << "Enter your choice: ";
      cin >> option;
   }
   if (option==0)
      return EXIT_SUCCESS;
   if (option>3)
      return EXIT_FAILURE;

// Define a grid
   Grid *g;
   if (option==3)
      g = new Grid(-1.,1.,-1.,1.,-1.,1.,nx,ny,nz);
   else
      g = new Grid(-1.,1.,-1.,1.,nx,ny);

// Choose initialization
   Vect<double> u(*g);
   initial(option,g,u);

// Solve the problem
   FastMarching fm(*g,u);
   fm.run();

// Compute residual to test accuracy
   cout << "Consistency error: " << fm.getResidual() << endl;

// Save solution for plotting
   saveField(u,*g,"u.pos",GMSH);
   cout << "Solution is stored in file 'u.pos' for plotting with gmsh" << endl;

// Free pointer and exit
   delete g;
   return EXIT_SUCCESS;
}


// Initial solution
void initial(int option,  Grid *g, Vect<double> &u)
{
   u = INFINITY;
   int nx=g->getNx(), ny=g->getNy(), nz=g->getNz();
   double eps=1.e-5;

   switch (option) {

      case 1:
         for (int i=1; i<=nx+1; ++i) {
            for (int j=1; j<=ny+1; ++j) {
               double x=g->getX(i), y=g->getY(j);
               if (fabs(x+0.5)<eps && fabs(y+0.5)<eps)
                  u(i,j) = 0.;
            }
         }
         break;

      case 2:
         {
         double a = 0.4;
         for (int i=1; i<=nx+1; ++i) {
            for (int j=1; j<=ny+1; ++j) {
               double x=g->getX(i), y=g->getY(j);
               if (x>-a && x<a && fabs(fabs(y)-a)<eps)
                  u(i,j) = 0.;
               if (y>-a && y<a && fabs(fabs(x)-a)<eps)
                  u(i,j) = 0.;
               if (x>-a && x<a && y>-a && y<a)
                  u(i,j) = -u(i,j);
            }
         }
         }
         break;

      case 3:
         for (int i=1; i<=nx+1; ++i) {
            for (int j=1; j<=ny+1; ++j) {
               for (int k=1; k<=nz+1; ++k) {
                  double x=g->getX(i), y=g->getY(j), z=g->getZ(k);
                  if (fabs(x+0.5)<eps && fabs(y+0.5)<eps && fabs(z+0.5)<eps)
                     u(i,j,k) = 0.;
               }
            }
         }
         break;
   }
}
