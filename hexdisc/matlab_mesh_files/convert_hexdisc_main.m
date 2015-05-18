function convert_hexdisc_main
% CONVERT Matlab style into Fortran style

disp('2-D Hexagon Disc');

% Start of user input
% User can change the amplification variable, the amount of vertices of the disc,
% the edgel variable (to create fine or coarse mesh), and the xmin, xmax, ymin, ymax
% variables if nonuniform initial bounds are wanted

% 6 vertices of hexagon
a=1;                           % Amplification factor for the x coordinates
verts=6;                        % Number of vertices (5=pentagon,6=hexagon,8=octagon,etc)
phi=(0:verts)'/verts*2*pi;      % Finds the angles at which all the vertices are placed
pfix=[a*cos(phi),a*sin(phi)];   % X and y coordinates at the cooresponding angles with amplification included

edgel=.02;                        % Initial smallest edge length
xmin=-a;                        % Initial lower bound in the x direction
xmax=a;                         % Initial upper bound in the x direction
ymin=-a;                        % Initial lower bound in the y direction
ymax=a;                         % Initial upper bound in the y direction

% End of user input

% 2-D triangle mesh
[p,t]=distmesh2d(@dpoly,@huniform,edgel,[xmin,ymin;xmax,ymax],[],pfix);

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
fid=fopen('../nnee_hexdisc.dat','w');
fprintf(fid,'%12i %12i %12i',nnode,nelem,nedge);
fclose(fid);

% node coordinates 
fid=fopen('../xy_hexdisc.dat','w');
fprintf(fid,'%18.15e %18.15e\n',p');
fclose(fid);

% global node label
fid=fopen('../node_hexdisc.dat','w');
fprintf(fid,'%12i %12i %12i\n',t');
fclose(fid);

% global edge label
fid=fopen('../edge_hexdisc.dat','w');
fprintf(fid,'%12i %12i\n',e');
fclose(fid);

% area of triangular elements
fid=fopen('../area_hexdisc.dat','w');
fprintf(fid,'%18.15e\n',ds');
fclose(fid);

% Output area and perimeter of 2-D hexagonal disc
side=sqrt((a*cos(0.)-a*cos((2.*pi)/verts))^2. + (a*sin(0.)-a*sin((2.*pi)/verts))^2.);
area=3.*sqrt(3.)/2.*side^2.;
perim=verts*side;
perim_to_area=perim/area;

fid=fopen('../attribs_hexdisc.txt','w');
fprintf(fid,'2-D Hexagon Disc Geometry\n\n');
fprintf(fid,'Perimeter: %15.10e\n',perim);
fprintf(fid,'Area: %15.10e\n',area);
fprintf(fid,'Perimeter-area ratio: %15.10f\n',perim_to_area);
fclose(fid);