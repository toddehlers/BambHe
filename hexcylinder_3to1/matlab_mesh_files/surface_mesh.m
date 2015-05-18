function surface_mesh
% MESH generate tetrahedarl surface mesh

load nnes.dat
load xyzn.dat
load sfnd.dat

nnode=nnes(1);
nelem=nnes(2);
nsfem=nnes(3);
p=xyzn';
e=sfnd';

for i=1:nsfem
    for j=1:3
        ip=e(j,i);
        x(j)=p(1,ip);
        y(j)=p(2,ip);
        z(j)=p(3,ip);
    end
    fill3(x,y,z,[0.4,0.0,1.0]);
    hold on;
end

camlight headlight;
axis equal;
axis off;


