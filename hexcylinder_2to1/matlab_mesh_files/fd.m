function d=fd(p,pv)
% distance function for hexagon cylinder 

global pheight;
global nheight;

np=size(p,1);
nvs=size(pv,1)-1;

ds=dsegment(p(:,1:2),pv);
d1=min(ds,[],2);
d1=(-1).^(inpolygon(p(:,1),p(:,2),pv(:,1),pv(:,2))).*d1;

z=p(:,3);
d2=z-pheight;
d3=-z-nheight;
d4=sqrt(d1.^2+d2.^2);
d5=sqrt(d1.^2+d3.^2);
d=dintersect(dintersect(d1,d2),d3);
ix=d1>0 & d2>0;
d(ix)=d4(ix);
ix=d1>0 & d3>0;
d(ix)=d5(ix);



