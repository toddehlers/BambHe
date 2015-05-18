function convert
% CONVERT Matlab style into Fortran style

disp('3-D hexagon cylinder');

% six vertices of hexagon
n=6;
phi=(0:n)'/n*2*pi;
pv=[cos(phi),sin(phi)];

% 3-D tetrahedral mesh
[p,t]=distmeshnd(@fd,@huniform,0.3,[-1,-1,-1;1,1,1],[],pv);

% find number of nodes and elements
nnode=size(p,1);
nelem=size(t,1);

% find all the surface triangles e in tetrahedral mesh [p,t]
e=surftri(p,t);
nsfem=size(e,1);

% compute the volume of the simplex elements in mesh [p,t]
v=simpvol(p,t);
v=abs(v);

% number of node, element, and surface element
fid=fopen('nnes.dat','w');
fprintf(fid,'%12i %12i %12i',nnode,nelem,nsfem);
fclose(fid);

% node coordinate
fid=fopen('xyzn.dat','w');
fprintf(fid,'%18.15f %18.15f %18.15f\n',p');
fclose(fid);

% global node label 
fid=fopen('node.dat','w');
fprintf(fid,'%12i %12i %12i %12i\n',t');
fclose(fid);

% global surface triangles label
fid=fopen('sfnd.dat','w');
fprintf(fid,'%12i %12i %12i\n',e');
fclose(fid);

% volume of each tetrahedarl element
fid=fopen('volu.dat','w');
fprintf(fid,'%18.15f\n',v');
fclose(fid);
