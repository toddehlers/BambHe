function he
% plot of solution u

load nnee.dat;
load node.dat;
load x.dat;
load y.dat;
load u.dat;

nelem=nnee(2);

for i=1:nelem
    for j=1:3
        ip=j+3*(i-1);
        xx(j)=x(ip);
        yy(j)=y(ip);
        uu(j)=u(ip);
    end
    fill3(xx,yy,uu,uu);
    hold on;
end

view(90,0)