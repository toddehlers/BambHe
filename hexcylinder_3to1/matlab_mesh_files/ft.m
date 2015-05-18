% plot of true f

load nnee.dat;
load x.dat;
load y.dat;

nelem=nnee(2);

axis([-1 1 -1 1 0.2 1.2]);
hold on;
grid on;

for i=1:nelem
    for j=1:3
        ip=j+3*(i-1);
        xx(j)=x(ip);
        yy(j)=y(ip);
        ff(j)=f(xx(j),yy(j));
    end
    fill3(xx,yy,ff,ff);
    hold on;
end

clear all;
