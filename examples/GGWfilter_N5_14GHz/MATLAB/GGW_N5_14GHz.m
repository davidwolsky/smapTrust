function R = GGW_N5_14GHz(x,xp,f)
% input parameters:
%   x       - vector of parameters [Nn,1]
%   xp      - vector of implicit parameters [Nq,1]
%   f       - frequency vector [Nm,1]
%
% output parameters:
%   R       - response vector only (no frequency), S1,1 in dB in this case

%% test inputs
% f = linspace(13.8e9,14.2e9,301);
% x = [0, 0, 0, 0, 0, 0];
% xp = [0];

%% Scale inputs
f = f.*1e0;

%% Decode inputs
L1_term = x(1);
L2_term = x(2);
L3_term = x(3);
s12_term = x(4);
s23_term = x(5);
x_term = x(6);

P = xp(1);
% as_fact = xp(2);
% bs_fact = xp(3);
% cs_fact = xp(4);
% ax_fact = xp(5);
% bx_fact = xp(6);
% cx_fact = xp(7);
% L1c = xp(8);
% L2c = xp(9);

%% Fixed parameters
as_fact = 1;
bs_fact = 1;
cs_fact = 1;
ax_fact = 1;
bx_fact = 1;
cx_fact = 1;

Pdeembed = 100*P;
% L1c = -1.04;
L1c = -4.8927;
L2c = L1c;
L3c = L1c;
 
L1 = (15.6 + L1c + L1_term).*1e-3;
L2 = (15.6 + L2c + L2_term).*1e-3;
L3 = (15.6 + L3c + L3_term).*1e-3;

s12_0 = 0.5119;
s23_0 = 0.7792;
s12 = s12_0 + s12_term;
s23 = s23_0 + s23_term;
as = 8.03.*as_fact;
bs = 2.364.*bs_fact;
cs = 1.682.*cs_fact;
k12 = as.*(s12 + cs).^-bs./100;
k23 = as.*(s23 + cs).^-bs./100;

x0 = 3.7620;
x = x0 + x_term;
ax = 0.3803.*ax_fact;
bx = 1.698.*bx_fact;
cx = -1.118.*cx_fact;
Qe = ax.*(x + cx).^(-bx).*1000;

c0 = 299792458;
lam = c0./f;

width = 15.8e-3;
d = 6.25e-3;
h = 1e-3;
height = d+h;
Lf = 4e-3;

%% Get circuit elements

[Z0_RWG,lam_g_RWG,lam_cutoff_RWG] = rectWG(f,width,height);
TXlineType = 'TEM';
% TXlineType = 'WG';
switch TXlineType
    case 'TEM'
        bet = 2.*pi./lam;
        Zline = ones(size(f)).*50;
    case 'WG'
        bet = 2.*pi./lam_g_RWG;
        Zline = Z0_RWG;
end

% FeedLine
% TODO_DWW: DDV: d2r not defined.
% Lphase = d2r(Pdeembed)./bet;
Lphase = deg2rad(Pdeembed)./bet;
Lfeed = Lf + Lphase;
Tfeed = TlineABCD(Zline,bet,Lfeed);

% Input coupler
TQe = KinvABCD(Zline./sqrt(Qe));

% Resonator 1
TR1 = TlineABCD(Zline,bet,L1);

% Coupler 12
TC12 = KinvABCD(Zline.*k12);

% Resonator 2
TR2 = TlineABCD(Zline,bet,L2);

% Coupler 23
TC23 = KinvABCD(Zline.*k23);

% Resonator 3
TR3 = TlineABCD(Zline,bet,L3);

%% Build circuit
[Ts,T] = deal(zeros(2,2,length(f)));
for ff = 1:length(f)
    T(:,:,ff) = Tfeed(:,:,ff)*TQe(:,:,ff)*TR1(:,:,ff)*TC12(:,:,ff)*TR2(:,:,ff)*TC23(:,:,ff)*TR3(:,:,ff)*TC23(:,:,ff)*TR2(:,:,ff)*TC12(:,:,ff)*TR1(:,:,ff)*TQe(:,:,ff)*Tfeed(:,:,ff);
end
Smat = ABCD2S(T,Zline,Zline);
R = squeeze(Smat(1,1,:));
% R = squeeze(Smat(1,1,:));

% plot(f./1e9,dB20(R)), grid on




