function convert
% CONVERT Matlab style into Fortran style

disp('2-D hexagon');

% 6 vertices of hexagon
phi=(0:6)'/6*2*pi;
pfix=[cos(phi),sin(phi)];

% 2-D triangle mesh
[p,t]=distmesh2d(@dpoly,@huniform,0.2,[-1,-1;1,1],pfix,pfix);

% find number of nodes and elements
nnode=size(p,1);
nelem=size(t,1);

% find all the boundary edges e in triangular mesh [p,t]
e=boundedges(p,t);
nedge=size(e,1);

% compute the area of the triangular elements in mesh [p,t]
ds=simpvol(p,t);
ds=abs(ds);

% number of nodes, elements, and edges
fid=fopen('nnee.dat','w');
fprintf(fid,'%12i %12i %12i',nnode,nelem,nedge);
fclose(fid);

% node coordinates 
fid=fopen('xy.dat','w');
fprintf(fid,'%18.15f %18.15f\n',p');
fclose(fid);

% global node label
fid=fopen('node.dat','w');
fprintf(fid,'%12i %12i %12i\n',t');
fclose(fid);

% global edge label
fid=fopen('edge.dat','w');
fprintf(fid,'%12i %12i\n',e');
fclose(fid);

% area of triangular elements
fid=fopen('area.dat','w');
fprintf(fid,'%18.15f\n',ds');
fclose(fid);