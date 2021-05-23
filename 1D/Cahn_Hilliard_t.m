% Cahn-Hilliard equation 1D
% $\partial_t u - D \partial_x (u(1-u)^\alpha\partial_x \mu) = 0$
% $\mu = -\gamma \Delta n + \psi^\prime(n)$
close all; clear all; format long;

lim = 'MC'; % choose flux limiter

%% Model parameters
D = 1;
gamma = 1e-4; % squared of length of diffuse interface
alpha=2; % exponent for mobility
nstar = 0.6; % equilibrium density for attraction/repulsion
k = 10; % constant to have a positive energy

b = @(n) n.*(1-n).^alpha; % mobility function
psip = @(n) (1-nstar)./(1-n) - n.^2 - (1-nstar).*n -(1-nstar); % derivative of potential
psi = @(n) -(1-nstar).*log(1-n) - n.^3./3 - (1-nstar).*n.^2./2 -(1-nstar).*n +k;

%% Space grid
nbx = 200;                          % Number of points in the domain
L = 1;     
xm=0;
xM=1;% Length of the domain
h = (xM-xm)/(nbx-1);                     % grid space
x = xm:h:xM;                            % discretization

%% Time discretization
T = 2;
dt = 1e-6;
t = 0:dt:T;  
CFL=0.0001;

%% Initial condition
% initial condition
r = -1 + (2)*rand(nbx,1);
n = abs(0.3+0.05*r);  % Random numbers
mu = zeros(nbx,1);
mu(2:nbx-1) = -gamma *(n(3:nbx)-2*n(2:nbx-1)+n(1:nbx-2))./h.^2 + psip(n(2:nbx-1));
mu(1) = -gamma*(n(2)-n(1))./h^2 + psip(n(1));
mu(nbx) = -gamma*(-n(nbx)+n(nbx-1))./h^2 + psip(n(nbx));

% initialisation flux at boundary of cells
fluxm_h = zeros(nbx+1,1); 

figure(1)
plot(x,n, '-*');
hold on 
plot(x,mu);
hold off;

tt = dt; inc = 1; inc2=1;
%%
while(tt < T)
    tt
    % the velocity here is minus the gradient of mu
    v=zeros(nbx+1,1);
    v(2:nbx) = -(mu(2:nbx)-mu(1:nbx-1))./h;
    
    % compute flux
    fluxm_h = upwind_flux_degenerate_mobility(n,v,nbx,h,dt,alpha);
    
    % CFL condition: adaptive time step
    dt = 1e-5; 
    dt = min(dt, CFL.*h./max(abs(v)));
    % compute solution
    % density
    n = n - dt*(fluxm_h(2:nbx+1)-fluxm_h(1:nbx))./h;
    % chemical potential
    mu(2:nbx-1) = -gamma* (n(3:nbx)-2*n(2:nbx-1)+n(1:nbx-2))./(h.^2) + psip(n(2:nbx-1));
    mu(1) = -gamma*(n(2)-n(1))./(h^2) + psip(n(1));
    mu(nbx) = -gamma*(-n(nbx)+n(nbx-1))./(h^2) + psip(n(nbx));
    
    % plot the solution during simulation 
    if mod(inc,1000)==0 
        figure(1)
        plot(x,n, '-*');
        
        
        figure(2)
        plot(x,mu, '-*')
        
        mass(inc2) = h*sum(n);
        figure(3)
        plot(mass)
        inc2=inc2+1;
    end
    inc = inc +1;
    tt = tt+dt;
end