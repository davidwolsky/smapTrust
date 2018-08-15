function Rc = coarseTest(xi,f)
% f = linspace(0,pi,101);
x1 = xi(1);
x2 = xi(2);
Rc = 1.0*atan( (x1*f') - x2 - 7.0*pi );
Rc = reshape(Rc,length(f),1);
