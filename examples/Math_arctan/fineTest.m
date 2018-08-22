function Rf = fineTest(xi,f)
% f = linspace(0,pi,101);
x1 = 1.0*xi(1) + 0.0 + 0.0j; 
x2 = 1.0*xi(2) - 1.0*pi + 0.0j;
Rf = 1.0*atan(x1*(f' - x2));
Rf = reshape(Rf,length(f),1);