function Rc = coarseTest(xi,f)
% f = linspace(0,pi,101);
x1 = xi(1) + 0.0j;
x2 = xi(2) + 0.0j;
Rc = 1.0*atan(x1*(f' - x2));
Rc = reshape(Rc,length(f),1);
