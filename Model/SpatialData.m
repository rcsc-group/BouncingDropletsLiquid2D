%Spatial Data

%% Physical Parameters

R = 0.2e-3; % Radius  Choose from: [0.05, 0.1,0.2,0.3,0.4,0.5,0.6,0.7]*1e-3;

% AIR 
rho_air = 1.225;      % air     % density air    kg/m^3
mu_air  = 1.825e-5;   % air     % viscosity air  kg/ms

%WATER 
rho_liq = 1000;      % water  % density liquids     kg/m^3 
sig     = 7.2e-2;    % water  % surface tension     kg/s^2
mu_liq  = 9.78e-4;   % water  % dynamic viscosity   kg/ms

We = 0.2;  % Choose from: [0.05, 0.1111,0.3403,0.2222,0.5,0.68,0.8]; 

% Compute Initial velocity 
W0 = -sqrt((We_fixed*sig)/(rho_liq*R));
%% Numerical Discretization 

%Discretization Choices
N_x = 2^13;                     % # Gridpoints along bath 
N_theta = 2^8;                  % # Gridpoints along droplet
dtheta = 2*pi/N_theta;          % Spacing of drop
dx = R*dtheta/4;                % Spacing of bath
theta = -pi:dtheta:pi;          % Grid on drop
L = dx*N_x;                     % Length of bath
x = -L/2:dx:L/2;                % Grid on bath
dk = (2*pi)/L;                  % Fourier modes spacing on bath
k = dk*[0:N_x/2 -N_x/2+1:-1];   % Fourier modes on bath

%Code Tidy Coefficient Names
Origin = N_x/2+1;               % Index of x=0
O_theta = N_theta/2+1;          % Index of theta=0

%Fourier mode Operators
Int = 1./(1i*k);        Int(Origin) = 0; Int(1)= 0;    % integrate wrt x
Dxx = (1i*k).^2;        Dxx(Origin) = 0;            % second order x deriv
Dz = abs(k);            Dz(Origin) = 0;             % First order z deriv
Dx = (1i*k);

%% Droplet Parameters

% Droplet Values
rho_drop =  rho_liq;   
sig_drop = sig;   
mu_drop = mu_liq; 


n_modes = (1:32);
damping_coeff = 150;

lambda_impact(1:32)= 2*n_modes(2)*damping_coeff*mu_drop/(rho_drop*R^2); %For During Impact 
lambda_liq = 2*n_modes.*((n_modes)-1).*mu_liq/(rho_drop*R^2); % for Free Flight
womega_sq = n_modes.*((n_modes).^2-1)*sig_drop/(rho_drop*R^3);


%% Non-dim parameters 
vare_liq = (mu_air/(rho_liq*R*abs(W0)))^(1/3);  %h/R < vare
We_liq = (rho_liq*R*W0^2)/sig; % Weber Number
Re_liq = rho_liq*R*abs(W0)/mu_liq; %Reynolds Number
Oh_liq = mu_liq/sqrt(rho_liq*R*sig); %Ohnesorge Number

vare_drop = (mu_air/(rho_drop*R*abs(W0)))^(1/3);  %h/R < vare
We_drop = (rho_drop*R*W0^2)/sig_drop; % Weber Number
Re_drop = rho_drop*R*abs(W0)/mu_drop(end); %Reynolds Number
Oh_drop = mu_drop/sqrt(rho_drop*R*sig_drop); %Ohnesorge Number


