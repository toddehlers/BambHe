
% plot of numerical reconstruction of f

load nnee.dat;
load x.dat;
load y.dat;
load fn_ex2.dat;

nelem=nnee(2);

axis([-1 1 -1 1 0.2 1.2]);
hold on;
grid on;

for i=1:nelem
    for j=1:3
        ip=j+3*(i-1);
        xx(j)=x(ip);
        yy(j)=y(ip);
        ff(j)=fn_ex2(ip);
    end
    fill3(xx,yy,ff,ff);
    hold on;
end

xlabel('x_1');
ylabel('x_2');


clear all;
