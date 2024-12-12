function [P,Grad,P_n,P_f,h,h_t,rho,rho_t] = Pressure_Solve(n_hat,phi_hat,Z,W,c,c_dot,Params)
%% Resolving the thin film elliptic problem to determine pressure within the air layer

%Unpacking Params
N_x = Params.N_x;
Dx = Params.Dx;
R = Params.R;
theta = Params.theta;
O_theta = Params.O_theta;
x = Params.x;
Origin = Params.Origin;
Dz = Params.Dz;
mu_liq = Params.mu_liq;
rho_liq = Params.rho_liq; 
mu_air = Params.mu_air;
Dxx = Params.Dxx;
Q = Params.Q;
rho_drop = Params.rho_drop;

%Initialise Empty Vectors
P = zeros(1,N_x);
Grad = zeros(1,N_x);
P_n = zeros(size(c_dot));
P_f = 0;

n = real(ifft(n_hat,'symmetric')); %CHANGE HERE
phi_x = real(ifft(Dx.*phi_hat,'symmetric'));

%% Calculate h

% Compute r(theta,t)
rho = R + sum(c(2:end).*cos((2:length(c))'.*theta)); %r(theta,t)
[rho_x,rho_z] = pol2cart(theta-pi/2,rho);

%Widest points
rho_max = find(diff(rho_x(O_theta:end))<0,1) -1; %Widest point on sphere
[var, x_max] = min(abs(rho_x(O_theta+rho_max) - x));

%find region of defined height
x_defined = x(Origin:x_max-1);
S_defined = interp1(rho_x(O_theta:O_theta+rho_max),rho_z(O_theta:O_theta+rho_max),x_defined,'spline');

%Find height along points
h_defined = S_defined+Z-n(Origin:x_max-1);

%find x^* (Last point that satisfies h(x)/R< epsilon, out from origin!
x_star = find(h_defined >= 0.01*R ,1)-1; 

%% Pressure Computation
if x_star > 4
    %% Calculaute h_t
    x_Lub = x_defined(1:x_star);
    h_Lub = h_defined(1:x_star);

    %Find corresponding theta_star
    [var,theta_star] = min(abs(x_defined(x_star)-rho_x(O_theta:O_theta + rho_max)));

    %Find ht
    rho_t =  sum(c_dot(2:end).*cos((2:length(c))'.*theta)); 
    [rho_t_x,rho_dt_z] = pol2cart(theta-pi/2,rho_t);
    S_dt= interp1(rho_x(O_theta:O_theta+rho_max),rho_dt_z(O_theta:O_theta+rho_max),x_Lub,'spline');
    n_dt = ifft(real( Dz.*phi_hat + (2*mu_liq/rho_liq)*Dxx.*n_hat),'symmetric');
    h_t = S_dt + W-n_dt(Origin:Origin+x_star-1);

    %% Calculate P(x,t)
    % Find Flux = int ht dx
    Flux = cumtrapz(x_Lub,h_t);

    % Find Slip term 
    Phi_x_Lub = phi_x(Origin:Origin+x_star-1);
    Slip = 0.5*Q*mu_air*Phi_x_Lub.*h_Lub.^(-2);

    % Find Grad = 12 mu *h^(-3) Flux + Slip
    Grad_Lub = Q*mu_air*Flux.*h_Lub.^(-3) + Slip;
    
    % Find P = int Grad dt + C
    P_Lub = flip(cumtrapz(flip(x_Lub),flip(Grad_Lub)));

    %Assign Full Vectors 
    Grad(Origin - x_star+1:Origin) = -flip(Grad_Lub);
    Grad(Origin:Origin+x_star-1) = Grad_Lub;
    P(Origin - x_star+1:Origin) = flip(P_Lub);
    P(Origin:Origin+x_star-1) = P_Lub;
    
    %% Obtain P(theta,t) from P(x,t) 
    %Polar convert back
    rho_Lub = rho_x(O_theta:O_theta+theta_star-1);
    P_theta = interp1(x_Lub,P_Lub,rho_Lub,'spline');
    P_Polar = zeros(size(theta));
    P_Polar(O_theta:O_theta+theta_star-1) = P_theta;
    P_Polar(O_theta - theta_star+1:O_theta) = flip(P_theta);
    
    %% Calculate P_n
    P_n = -(1:size(c_dot))'.*Pmode(P_Polar,(1:size(c_dot))')/(rho_drop*R);
    
    %% Calculate P_f
    Force = trapz([-flip(x_Lub(2:end)),x_Lub],[flip(P_Lub(2:end)),P_Lub]);
    mass = pi*rho_drop*R^2;
    P_f = Force/mass;
    
    %% Send h and h_t to plotting
    h  = h_Lub; 
else
    P = zeros(1,N_x);
    Grad = zeros(1,N_x);
    P_n = zeros(size(c_dot));
    P_f = 0;
    h = h_defined;
    h_t = zeros(size(h_defined));
    rho_t = zeros(size(rho));
end


function p_n = Pmode(P,n)
    % Sphere surface Forcing [-pi,pi]
    p_n = trapz(theta(1:end-1), P(1:end-1).*cos(n.*theta(1:end-1)),2)/pi;       
end


end