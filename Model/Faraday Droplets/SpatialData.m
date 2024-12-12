%SpatilData

R = 0.83e-3;
W0 = -0.1;

%% Numerical discretization 

%Discretization Choices
N_x = 2^7*3;                    % # Gridpoints along bath 
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

%Fourier mode stuff
Int = 1./(1i*k);        Int(Origin) = 0; Int(1)= 0;    % integrate wrt x
Dxx = (1i*k).^2;        Dxx(Origin) = 0;            % second order x deriv
Dz = abs(k);            Dz(Origin) = 0;             % First order z deriv
Dx = (1i*k);

%% Physical Parameters
% AIR 
rho_air = 1.225;      % air     % density air kg/m^3
mu_air  = 1.825e-5;   % air     % viscosity air  kg/ms

%Silicone oil
rho_liq = 949;            % density liquids     kg/m^3 
sig     = 2.06e-2;        % surface tension     kg/s^2
mu_liq  = 0.8205*1.9e-2;  % dynamic viscosity   kg/ms

% Droplet Values
rho_drop =  rho_liq;   
sig_drop = sig;   
mu_drop = mu_liq; 

n_modes = (1:32);
lambda_impact= 2*n_modes.*((n_modes)-1).*150*mu_drop/(rho_drop*R^2);
%No damping required as increased viscosity is sufficient
lambda_liq = 2*n_modes.*((n_modes)-1).*mu_liq/(rho_drop*R^2);
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


