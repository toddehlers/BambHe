function convert_rounddisc_main
% CONVERT Matlab style into Fortran style

disp('2-D Unit Disc')

% Start of user input
% User can change the radius variable, the edgel variable (to create fine or coarse mesh),
% and the xmin, xmax, ymin, ymax variables if nonuniform initial bounds are wanted

rad=1;                                                  % Radius of the circular disc
fd=inline(['dcircle(p,0,0,' num2str(rad) ')'],'p');     % distance function

edgel=.03;                                              % Initial smallest edge length
xmin=-rad;                                              % Initial lower bound in the x direction
xmax=rad;                                               % Initial upper bound in the x direction
ymin=-rad;                                              % Initial lower bound in the y direction
ymax=rad;                                               % Initial upper bound in the y direction

% End of user input

% 2-D triangle mesh
[p,t]=distmesh2d(fd,@huniform,edgel,[xmin,ymin;xmax,ymax],[]);

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
fid=fopen('../nnee_rounddisc.dat','w');
fprintf(fid,'%12i %12i %12i',nnode,nelem,nedge);
fclose(fid);

% node coordinates 
fid=fopen('../xy_rounddisc.dat','w');
fprintf(fid,'%18.15e %18.15e\n',p');
fclose(fid);

% global node label
fid=fopen('../node_rounddisc.dat','w');
fprintf(fid,'%12i %12i %12i\n',t');
fclose(fid);

% global edge label
fid=fopen('../edge_rounddisc.dat','w');
fprintf(fid,'%12i %12i\n',e');
fclose(fid);

% area of triangular elements
fid=fopen('../area_rounddisc.dat','w');
fprintf(fid,'%18.15e\n',ds');
fclose(fid);

% Output area and perimeter of 2-D round disc
area=pi*rad^2.;
perim=2.*pi*rad;
perim_to_area=perim/area;

fid=fopen('../attribs_rounddisc.txt','w');
fprintf(fid,'2-D Round Disc Geometry\n\n');
fprintf(fid,'Perimeter: %15.10e\n',perim);
fprintf(fid,'Area: %15.10e\n',area);
fprintf(fid,'Perimeter-area ratio: %15.10f\n',perim_to_area);
fclose(fid);