%% Main body of code to simulate a 2D bouncing droplet. 

clear;
close all;

fullrun = tic;

Q = 12;

g= 9.81;
Gamma = 3;
omega = 80;


nprint = 4000;
t_final = 1.2;

SpatialData;
save('Params.mat') %Allow all parameters to be loaded into 
Params = load('Params.mat'); % Converts the mat into a struct.
File_suffix = input("Suffix for results: ",'s');


%% Intial Values
n = zeros(1,N_x);
phi = zeros(1,N_x);
Z = 1.1*R;
W = W0;
c = zeros(32,1);
c_dot = zeros(32,1);
i=1; % loop index

n_hat = fft(n);
phi_hat = fft(phi);
[P,Grad,P_n,P_f,h,h_t,rho,rho_t] = Pressure_Solve(n_hat,phi_hat,Z,W,c,c_dot,Params);

%% TimeStepping 
%Initial Conditions 
y0 = [n_hat,phi_hat,Z,W,c',c_dot']';
tspan = 0:t_final/nprint:t_final;

RHS = @(t,y) FuncM(t,y,Params);
%Loop
opts = odeset('RelTol',1e-5,'AbsTol',1e-9); %1e-1, 1e-3 causes 23s to throw out complex values for rho
[t,y] = ode45(@(t,y)RHS(t,y),tspan,y0,opts);


%% reconstructing Pressure

for i = 1:length(t)
    n_hat = y(i,1:N_x);
    phi_hat = y(i,N_x+1:2*N_x);
    Z = y(i,2*N_x+1);
    W = y(i,2*N_x+2);
    c = y(i,2*N_x+3:2*N_x+34)';
    c_dot = y(i,2*N_x+35:end)';
    [P,Grad,P_n,P_f,h,h_t,rho,rho_t] = Pressure_Solve(n_hat,phi_hat,Z,W,c,c_dot,Params);
    Pressure(i,:) = P;
    PressureGrad(i,:) = Grad;
    Force(i) = P_f;
    height{i} = h;
    Height_t{i} = h_t;
    Rho(i,:) = rho;
    Rho_t(i,:) = rho_t;
    Eta(i,:)= real(ifft(n_hat));
    Phi(i,:) = real(ifft(phi_hat));
    Position(i) = Z;
    Speed(i) =W;
    C(:,i) = c;
    C_dot(:,i) = c_dot;
end


File.('T') = tspan;
File.('x') = x;
File.('eta') = Eta;
File.('phi') = Phi;
File.('Z') = Position;
File.('W') = Speed;
File.('P') = Pressure;
File.('Grad') = PressureGrad; 
File.('h') = height;
File.('h_t') = Height_t;
File.('rho') = Rho;
File.('rho_t') = Rho_t;  
File.('Force') = Force;
File.('Params') = Params;
File.('C') = C;
File.('C_dot') = C_dot;
save(['Run',File_suffix,'.mat'],'-struct','File')

toc(fullrun)

%% Functions

function [y] =  FuncM(t,y,Params)
    N_x = Params.N_x;
    %Separating out Values
    n_hat = y(1:N_x)';
    phi_hat = y(N_x+1:2*N_x)';
    Z = y(2*N_x+1);
    W = y(2*N_x+2);
    c = y(2*N_x+3:2*N_x+34);
    c_dot = y(2*N_x+35:end);
    %Compute evolution
    [nt, phit, Zt,Wt,ct,c_dot_t] = Func(n_hat, phi_hat, Z, W,c,c_dot,t, Params);
    %Reconstruct Vetor
    y = [nt,phit,Zt,Wt,ct',c_dot_t']';
    
    disp(["t="+num2str(t,'%.10f'),"Speed="+num2str(W,'%.6f')])
end

