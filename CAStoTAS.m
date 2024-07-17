
%
%   Flight Management and Procedure Design
%
%   Conversion from CAS to TAS
%
%   Copyright(c) Chair for Flight Guidance and Air Traffic
%   Institute of Aeronautics and Astronautics
%   TU Berlin
%
function [Vtas] = CAStoTAS(Vcas, rho, p)

% Constants
R = 287.05287;
kappa = 1.4;
mu = (kappa-1)/kappa;

% Sea-level
[rho0, p0, temp0, soundSpeed0] = Atmos(0.0);

% Compute TAS from CAS
temp1 = (1+(mu/2)*(rho0/p0)*Vcas.^2);
temp2 = temp1.^(1/mu) - 1;
temp3 = 1 + (p0./p).*temp2;
temp4 = temp3.^mu - 1;
temp5 = (2/mu)*(p./rho).*temp4;

% Compute airspeed
Vtas = sqrt(temp5);
