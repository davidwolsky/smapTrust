function [] = SM_MATLABtest(type)

if nargin < 1
    type = 'c';
end

opts.getA = 0;
opts.getB = 0;
opts.getc = 0;
opts.getG = 0;
opts.getxp = 0;
opts.getF = 0;
opts.getd = 0;
opts.getE = 0;

opts.optsFminS = optimset('display','none');


% x_init = [0.5,1.5,0.03,0.5,0.5]';
x_init = [0.5,1.5,0.03]';
% x_init = [0.5,1.5]';
xp_init = [2,1]';

Nn = length(x_init);
Nq = length(xp_init);

% xp_init = [];
fmin = 2;
fmax = 8;  
Nf = 101;
f = linspace(fmin,fmax,Nf)';

Nm = length(f);

% opts.ximin = 0.5*x_init;
% opts.ximax = 1.8*x_init;

% Set up fine model
Mf.path = [pwd,'\'];
Mf.name = @modelF;
Mf.solver = 'MATLAB';
Mf.params = {'x','f'};
% Mf.ximin = zeros(Nn,1);
Mf.ximin = -ones(Nn,1).*2;
Mf.ximax = ones(Nn,1).*2;
Mf.freq = f;

% Set up coarse model 
Mc.path = [pwd,'\'];
Mc.name = @modelC;  
Mc.solver = 'MATLAB';
Mc.params = {'x','xp','f'}; % Must be this order (and names) for MATLAB sims.  Can omit xp or f though...
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.xpmin = xp_init.*0.5;
Mc.xpmax = xp_init.*3;
Mc.freq = f;

% Set up the SM
% The initial SM structure
Sinit = [];
Sinit.xp = xp_init;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
SMopts.optsFminS = optimset('display','none','TolX',1e-6,'TolFun',1e-6);
SMopts.optsPBIL.Nfeval = 5000;
SMopts.wk = 1.0;

% Set up the optimization
OPTopts.ximin = Mf.ximin;
OPTopts.ximax = Mf.ximax;
OPTopts.Ni = 5;
OPTopts.TRNi = OPTopts.Ni*2;
% OPTopts.TRNi = 1;
OPTopts.Rtype = {'Gen'};
% OPTopts.globOpt = 1;
OPTopts.globOpt = 0;
OPTopts.goalType = {'minimax'};
OPTopts.goalResType = {'Gen'};
% OPTopts.goalVal = {-20};
OPTopts.goalWeight = {1};
OPTopts.goalStart = {3};
OPTopts.goalStop = {7};
OPTopts.errNorm = {2};
OPTopts.optsPBIL.display =  'none'; 
OPTopts.optsPBIL.Nfeval = 5000;
OPTopts.optsPBIL.Nbest = 10; % DOM
OPTopts.M_PBIL = 4;
OPTopts.optsFminS = optimset('display','none');
OPTopts.TolX = 10e-4;
OPTopts.eta1 = 0.05;
OPTopts.eta2 = 0.9;
OPTopts.alp1 = 2.5;
OPTopts.alp2 = 0.25;
OPTopts.testEnabled = 0;
% OPTopts.optsFminS = optimset('MaxFunEvals',10,'display','iter');



%% Build the models


modPar = [];
if any(type == 'A')
    A = rand(1);
    modPar.A = A;
    SMopts.getA = 1;
end
if any(type == 'B')
    B = rand(Nn)./500;
    Bd = eye(Nn) + eye(Nn).*(rand(Nn) - 0.5) ;
    B = B + Bd;
    modPar.B = B;
    SMopts.getB = 1;
end
if any(type == 'c')
    c = (rand(Nn,1) - 0.5)./5;
%     c = x_init./10;
    modPar.c = c;
    SMopts.getc = 1;
end
if any(type == 'G')
    G = rand(Nq,Nn);
    modPar.G = G;
    SMopts.getG = [1,0,1;0,1,0];
end
if any(type == 'p')
    xp = rand(Nq,1);
    modPar.xp = xp;
    SMopts.getxp = 1;
end
if any(type == 'F')
    F = rand(2,1)./10;
    modPar.F = F;
    SMopts.getF = 1;
end

% Sinit.c = c.*1.01;

% Rc = modelC(x_init,xp_init,f);
% Rfi = model(x_init,xp_init,f,modPar);
% Rfit = modelF(x_init,f);
% S.coarse = @modelC;
% S.xp = xp_init;
% S.f = f;
% Si = buildSurr(x_init,Rfi,S,opts);
% Rs = evalSurr(x_init,Si);
% 
% 
% plot(f,Rc,'k'), grid on, hold on
% plot(f,Rfi,'r')
% plot(f,Rfit,'r.')
% plot(f,Rs,'b--')


%% Run the main loop
[Ri,Si,Pi,Ci,Oi,Li,Ti] = SMmain(x_init,Sinit,SMopts,Mf,Mc,OPTopts);

keyboard;


function Rf = modelF(xi,f)
Rf = model(xi,xp_init,f,modPar);
end

end


function Rc = modelC(xi,xpi,f)
Rc = model(xi,xpi,f);
end


function R = model(xi,xpi,fi,modPar)

if isempty(xpi)
    xpi = 0;
end
if isempty(fi)
    fi = linspace(-5,5,1001);
end

Nn = length(xi);
Nq = length(xpi);
Nm = length(fi);


A = 1;
B = diag(ones(Nn,1));
c = zeros(Nn,1);
G = zeros(Nq,Nn);
F = [1,0]';

if nargin > 3
    if isfield(modPar,'A'), A = modPar.A; end
    if isfield(modPar,'B'), B = modPar.B; end
    if isfield(modPar,'c'), c = modPar.c; end
    if isfield(modPar,'G'), G = modPar.G; end
    if isfield(modPar,'F'), F = modPar.F; end
end

x = B*xi + c;
xp = G*x + xpi;
f = F(1).*fi + F(2);

% R = sin(f);
g = 0;
for nn = 1:Nn-1
%     R = R + x(nn).*sin((2*nn-1).*f) + (x(nn) - 2).^(nn./3);
%     R = R + x(nn).*sin((2*nn-1).*f) + x(nn).*(f./100).^(Nn-nn);
% R = R + (f - x(nn)).^(2);
g = g + 100.*(x(nn+1) - x(nn).^2).^2 + (x(nn) - 1).^2;
% g = g + x(nn).^2;
end
R = abs((f - 5).^(2+log(g+1)));

for qq = 1:Nq
    R = R + 0*cos(f + xp(qq));
end


R = A.*reshape(R,Nm,1);
end