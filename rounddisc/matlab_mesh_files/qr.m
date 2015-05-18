
% plot of numerical reconstruction of q

load qn_ex1.dat;

tmin=0.0;
tmax=1.0;
nt=20;
dt=(tmax-tmin)/nt;
ht=tmin:dt:tmax;

for i=1:nt+1
    t=tmin+(i-1)*dt;
    qt(i)=q(t);
end

axis([0 1 0.5 1.1]);
hold on;

plot(ht,qt,'r');
hold on;
plot(ht,qn_ex1,'bo');
hold on;
box on;

legend('true diffusivity','reconstructed diffusivity');
xlabel('Time (t)');

clear all;
