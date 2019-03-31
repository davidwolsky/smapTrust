function [T] = JinvABCD(J)

% function [T] = JinvABCD(J)
% Shunt resistor ABCD matrix of size [2,2,Nf], where Nf is the number of
% frequencies as reflected in the length of J.
%
% Inputs:
% J - Inverter parameter (can, in general, be function of frequency)

Nf = length(J);
T(1,1,:) = zeros(Nf,1);
T(1,2,:) = 1./(1i.*J);
T(2,1,:) = -1i.*K;
T(2,2,:) = zeros(Nf,1);

