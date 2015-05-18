clear all;
load qn.dat;

tmin=0.0;
tmax=1.0;
nt=20;
dt=(tmax-tmin)/nt;
ht=tmin:dt:tmax;

for i=1:nt+1
    t=tmin+(i-1)*dt;
    qt(i)=q(t);
end

plot(ht,qt,'r-o');
hold on;
plot(ht,qn,'bs');
hold on;