function [Z0,lam_g,lam_cutoff] = rectWG(freq,width,height,m,n,eps_r,mu_r)

% function [Z0] = ridgeWG(freq,Wg,Hg,Wr,Hr,type,eps_r,mu_r)
% Calculates the wave impedance of rectangular waveguide
% Much expansion possible!
% Inputs:
% freq - frequency in Hz
% width - guide width in m 
% height - guide height in m 
% m - width mode number {1}
% n - height mode number {0}
% eps_r - relative permittivity
% mu_r - relative permeability

% Equations from Balanis 

if nargin < 4
    m = 1;
    n = 0;
    eps_r = 1;
    mu_r= 1;
elseif nargin < 5
    n = 0;
    eps_r = 1;
    mu_r= 1;
elseif nargin < 6
    eps_r = 1;
    mu_r= 1;
elseif nargin < 7
    mu_r = 1;
end

assert(~(m~=0 && n~=0),'m and n cannot both be zero');

%TODO_DWW: DDDV: Where do these come from?
% load EMconstants
eps0 = 8.85418781761e-12;
mu0 = pi*4e-7;

ep = eps0*eps_r;
mu = mu0*mu_r;

% Rename variables to fit paper
a = width;
b = height;

c = sqrt(1/(ep*mu));
lam = c./freq;
bet = 2*pi./lam;
eta = sqrt(mu/ep);

bet_x = m*pi./a;
bet_y = n*pi./b;
bet_c = sqrt(bet_x.^2 + bet_y.^2);

lam_x = 2*pi./bet_x;
lam_y = 2*pi./bet_y;
lam_cutoff = 2.*pi./bet_c;

fc = bet_c./(2*pi*sqrt(mu*ep));

fsqrt = sqrt(1 - (fc./freq).^2);
Z0 = eta./fsqrt;
lam_g = lam./fsqrt;



