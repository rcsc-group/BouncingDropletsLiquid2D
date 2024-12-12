function [nt, phi_t, Zt,Wt,ct,ct_dot] = Func(n_hat, phi_hat, Z, W,c,c_dot,t, Params)
%% Right hand side of ODE in pseudospectral solving of the Lubrication-mediated Model

%% Initialising 
c_modes = 1:length(c);
[P,Grad,P_n,P_f,h,h_t,rho,rho_t] = Pressure_Solve(n_hat,phi_hat,Z,W,c,c_dot,Params);
P_hat = fft(P);

% Unpacking Params
Dz = Params.Dz;
mu_liq = Params.mu_liq;
rho_liq = Params.rho_liq;
Dxx = Params.Dxx;
sig = Params.sig;
g = Params.g;
lambda_impact = Params.lambda_impact;
lambda_liq = Params.lambda_liq;
womega_sq = Params.womega_sq;
R = Params.R;

%% ODEs
%Free Surface Waves
nt = Dz.*phi_hat + (2*mu_liq/rho_liq)*Dxx.*n_hat;
phi_t = (-1/rho_liq)*P_hat + (sig/rho_liq)*Dxx.*n_hat  + (2*mu_liq/rho_liq)*Dxx.*phi_hat - g*n_hat;
% Centre of Mass Trajectory
Zt = W;
Wt = P_f - g;
%Droplet Deformation
ct = c_dot;

% If resolving solid cases, multiply RHS of both statements below by zero to not
% permit oscillation. 
if min(h)<0.01*R 
    ct_dot = P_n - 2*lambda_impact'.*c_dot - womega_sq'.*c;
else
    ct_dot = P_n - 2*lambda_liq'.*c_dot - womega_sq'.*c;
end

end