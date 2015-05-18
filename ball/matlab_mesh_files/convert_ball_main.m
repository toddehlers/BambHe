function convert_ball_main
% CONVERT Matlab style into Fortran style

disp('3-D unit ball');

% Start of user input
% User can change the radius variable, the edgel variable (to create fine or coarse mesh),
% and the xmin, xmax, ymin, ymax, zmin, zmax variables if nonuniform initial bounds are wanted

rad=1;                                                  % Radius of the sphere
fd=inline(['sqrt(sum(p.^2,2))-' num2str(rad)],'p');     % distance function

edgel=.03;                                              % Initial smallest edge length
xmin=-rad;                                              % Initial lower bound in x direction
xmax=rad;                                               % Initial upper bound in x direction
ymin=-rad;                                              % Initial lower bound in y direction
ymax=rad;                                               % Initial upper bound in y direction
zmin=-rad;                                              % Initial lower bound in z direction
zmax=rad;                                               % Initial upper bound in z direction

% End of user input

% 3-D triangle mesh
[p,t]=distmeshnd(fd,@huniform,edgel,[xmin,ymin,zmin;xmax,ymax,zmax],[]);

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
fid=fopen('../nnes_ball.dat','w');
fprintf(fid,'%12i %12i %12i',nnode,nelem,nsfem);
fclose(fid);

% node coordinate
fid=fopen('../xyzn_ball.dat','w');
fprintf(fid,'%18.15e %18.15e %18.15e\n',p');
fclose(fid);

% global node label 
fid=fopen('../node_ball.dat','w');
fprintf(fid,'%12i %12i %12i %12i\n',t');
fclose(fid);

% global surface triangles label
fid=fopen('../sfnd_ball.dat','w');
fprintf(fid,'%12i %12i %12i\n',e');
fclose(fid);

% volume of each tetrahedarl element
fid=fopen('../volu_ball.dat','w');
fprintf(fid,'%18.15e\n',v');
fclose(fid);

% Output total volume, surface area, surface area to volume ratio and perimeter of the geometry
volume=(4./3.)*pi*rad^3.;
surf_area=4.*pi*rad^2.;
sa_to_vol=surf_area/volume;
perim=2.*pi*rad;

fid=fopen('../attribs_ball.txt','w');
fprintf(fid,'3-D Unit Ball Geometry\n\n');
fprintf(fid,'Surface area: %15.10e\n',surf_area);
fprintf(fid,'Volume: %15.10e\n',volume);
fprintf(fid,'Surface area-volume ratio: %15.10f\n',sa_to_vol);
fprintf(fid,'Perimeter: %15.10e\n',perim);
fclose(fid);