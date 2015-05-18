function w=f(x,y)

r0=0.2;
r=sqrt(x*x+y*y);

if r <= 1.0-r0
    w=1.0;
else
    w=0.5+0.5*(1.0-r)/r0;
end
