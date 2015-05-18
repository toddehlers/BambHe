function convert_hexcylinder_main
% CONVERT Matlab style into Fortran style

disp('3-D hexagon cylinder');

% global declaration of pheight and nheight (variables set in fd.m)
% pheight is the length of the cylinder above 0
% nheight is the length of the cylinder below 0
global pheight;
global nheight;

% Start of user input
% User can change the pheight and nheight variables to adjust the height of the cylinder,
% the amplification factor, the amount of vertices of the cylinder geometry,
% the edgel variable (to create fine or coarse mesh), and the xmin, xmax, ymin, ymax, zmin, zmax
% variables if nonuniform initial bounds are wanted

% The height of the entire cylinder can be found by adding pheight and nheight
pheight=1;                              % The length of the cylinder in the positive z direction
nheight=1;                              % The length of the cylinder in the negative z direction

% six vertices of hexagon
a=1;                                    % Amplification factor in the x and y directions for the hexagon
verts=6;                                % Number of vertices (5=pentagon,6=hexagon,8=octagon,etc)
phi=(0:verts)'/verts*2*pi;              % Angles at which the points will be placed
pv=[a*cos(phi),a*sin(phi)];             % X and Y coordinates of vertices based on angle and amplification

edgel=.08;                              % Initial smallest edge length
xmin=-a;                                % Initial lower bound in x direction
xmax=a;                                 % Initial upper bound in x direction
ymin=-a;                                % Initial lower bound in y direction
ymax=a;                                 % Initial upper bound in y direction
zmin=-nheight;                          % Initial lower bound in z direction
zmax=pheight;                           % Initial upper bound in z direction

% End of user input

% 3-D tetrahedral mesh
[p,t]=distmeshnd(@fd,@huniform,edgel,[xmin,ymin,zmin;xmax,ymax,zmax],[],pv);

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
fid=fopen('../nnes_hexcylinder.dat','w');
fprintf(fid,'%12i %12i %12i',nnode,nelem,nsfem);
fclose(fid);

% node coordinate
fid=fopen('../xyzn_hexcylinder.dat','w');
fprintf(fid,'%18.15e %18.15e %18.15e\n',p');
fclose(fid);

% global node label 
fid=fopen('../node_hexcylinder.dat','w');
fprintf(fid,'%12i %12i %12i %12i\n',t');
fclose(fid);

% global surface triangles label
fid=fopen('../sfnd_hexcylinder.dat','w');
fprintf(fid,'%12i %12i %12i\n',e');
fclose(fid);

% volume of each tetrahedarl element
fid=fopen('../volu_hexcylinder.dat','w');
fprintf(fid,'%18.15e\n',v');
fclose(fid);

% Output volume, surface area, surface area to volume ratio perimeter of hexagon cylinder
side=sqrt((a*cos(0.)-a*cos((2.*pi)/verts))^2.+(a*sin(0.)-a*sin((2.*pi)/verts))^2.);     % Simply using distance formula to find length of sides
volume=3.*sqrt(3.)/2.*side^2.*(pheight+nheight);
surf_area=3.*sqrt(3.)*side^2.+verts*side*(pheight+nheight);
sa_to_vol=surf_area/volume;
perim=verts*side;

fid=fopen('../attribs_hexcylinder.txt','w');
fprintf(fid,'3-D Hexagon Cylinder Geometry\n\n');
fprintf(fid,'Surface area: %15.10e\n',surf_area);
fprintf(fid,'Volume: %15.10e\n',volume);
fprintf(fid,'Surface area-volume ratio: %15.10f\n',sa_to_vol);
fprintf(fid,'Perimeter: %15.10e\n',perim);
fclose(fid);