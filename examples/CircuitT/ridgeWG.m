function [Z0,lam_g,lam_cutoff] = ridgeWG(freq,Wg,Hg,Wr,Hr,type,eps_r,mu_r)

% function [Z0] = ridgeWG(freq,Wg,Hg,Wr,Hr,type,eps_r,mu_r)
% Calculates the characteristic impedance of rectangular ridge waveguide
% Much expansion possible!
% Inputs:
% freq - frequency in Hz
% Wg - guide width in m 
% Hg - guide height in m 
% Wr - ridge width in m 
% Hr - ridge height in m 
% type - 1 for single, 2 for double ridge
% eps_r - relative permittivity
% mu_r - relative permeability

% Equations from Hoefer 1982


if nargin < 6
    type = 1;
    eps_r = 1;
    mu_r= 1;
elseif nargin < 7
    type = 1;
    mu_r = 1;
elseif nargin < 8
    mu_r = 1;
end

load constants

ep = eps0*eps_r;
mu = mu0*mu_r;

% Rename variables to fit paper
a = Wg;
s = Wr;
b = 2*Hg/type;
d = 2*(Hg-Hr)/type;

% Accuracy checks
if d/b < 0.01 || d/b > 1
    warning('Decreased accuracy for d/b outside [0.01,1]');
end
if b/a < 0 || b/a > 1
    warning('Decreased accuracy for b/a outside [0,1]');
end
if s/a < 0.01 || s/a > 0.45
    warning('Decreased accuracy for s/a outside [0,0.45]');
end

lam = c0./freq;

% Cutoff of open WG
b_lamcr = b./(2.*(a-s))./sqrt(1 + 4./pi.*(1 + 0.2.*sqrt(b./(a-s))).*b./(a-s).*log(csc(pi.*d./(2.*b))) + (2.45 + 0.2.*s./a).*s.*b./(d.*(a-s)));
lam_cutoff = b./b_lamcr;

lam_g = lam./sqrt(1 - (lam./lam_cutoff).^2);

% Impedance (voltage/current based)
th1 = pi.*s./lam_cutoff;
th2 = pi.*(a-s)./lam_cutoff;
B0_Y0 = 2.*b./lam_cutoff.*log(csc(pi.*d./(2.*b)));
Z0inf = 120*pi.^2.*(b./lam_cutoff)./(b./d.*sin(th1) + (B0_Y0 + tan(th2./2)).*cos(th1));
Z0 = Z0inf./sqrt(1 - (lam./lam_cutoff).^2)./2.*type;


