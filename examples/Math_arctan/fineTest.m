function Rf = fineTest(xi,f)
% f = linspace(0,pi,101);
x1 = 0.4*xi(1);
x2 = 1.0*xi(2) - 4.0*pi;
% x1 = 1.0*xi(1);
% x2 = 1.0*xi(2) + 2.0*pi;
Rf = 1.0*atan( (x1*f') - x2 - 7.0*pi );
Rf = reshape(Rf,length(f),1);