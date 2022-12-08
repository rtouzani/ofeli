%==============================================================================
%
%                                    O  F  E  L  I
%
%                           Object  Finite  Element  Library
%
%==============================================================================
%
%   Copyright (C) 1998 - 2023 Rachid Touzani
%
%   This file is part of OFELI.
%
%   OFELI is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   OFELI is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License
%   along with OFELI. If not, see <http://www.gnu.org/licenses/>.
%
%==============================================================================
%
%                Matlab function to export mesh matrices to OFELI
%
%   Mesh properties include matrices p, e, t, and b
%   This function stores these data in a given file that can be read by using
%   the OFELI function 'Matlab2OFELI' and then acquire as a Mesh instance.
%   If the file name is 'myfile.dat', the mesh converter 'cmesh' will create
%   the mesh MDF file 'myfile.dat.m'
%
%==============================================================================*/

function Matlab2OFELI(p,e,t,b,file)
fm = fopen(file,'w');

% Matrix p
fprintf(fm,'%d\n',size(p,2));
for i=1:1:size(p,2)
   fprintf(fm,'%g %g\n',p(1,i),p(2,i));
end

% Matrix e
fprintf(fm,'%d\n',size(e,2));
for i=1:1:size(e,2)
   for j=1:1:7
      fprintf(fm,'%d  ',e(j,i));
   end
   fprintf(fm,'\n');
end

% Matrix t
fprintf(fm,'%d\n',size(t,2));
for i=1:1:size(t,2)
   for j=1:1:4
      fprintf(fm,'%d  ',t(j,i));
   end
   fprintf(fm,'\n');
end

% Matrix b
fprintf(fm,'%d\n',size(b,2));
for i=1:1:size(b,2)
   N = b(1,i);
   M = b(2,i);
   fprintf(fm,'%d  %d\n',N,M);
   for j=3:1:3+N*N-1
      fprintf(fm,'%d  ',b(j,i));
   end
   fprintf(fm,'\n');
   for j=3+N*N:1:3+N*N+N-1
      fprintf(fm,'%d  ',b(j,i));
   end
   fprintf(fm,'\n');
   for j=3+N*N+N:1:3+N*N+N+M*N-1
      fprintf(fm,'%d  ',b(j,i));
   end
   fprintf(fm,'\n');
   for j=3+N*N+N+M*N:1:3+N*N+N+M*N+M-1
      fprintf(fm,'%d  ',b(j,i));
   end

   fprintf(fm,'\n');
   for j=3+N*N+N+M*N+M:1:3+N*N+N+M*N+M+N*N-1
      s = sprintf('%s',b(j,i));
      fprintf(fm,'%s  ',s);
   end
   fprintf(fm,'\n');
   for j=3+N*N+N+M*N+M+N*N:1:3+N*N+N+M*N+M+N*N+N-1
      s = sprintf('%s',b(j,i));
      fprintf(fm,'%s  ',s);
   end
   fprintf(fm,'\n');
   for j=3+N*N+N+M*N+M+N*N+N:1:3+N*N+N+M*N+M+N*N+N+M*N-1
      s = sprintf('%s',b(j,i));
      fprintf(fm,'%s  ',s);
   end
   fprintf(fm,'\n');
   for j=3+N*N+N+M*N+M+N*N+N+M*N:1:3+N*N+N+M*N+M+N*N+N+M*N+M-1
      s = sprintf('%s',b(j,i));
      fprintf(fm,'%s  ',s);
   end
   fprintf(fm,'\n');
end

fclose(fm);

