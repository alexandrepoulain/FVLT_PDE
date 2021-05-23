% Keller-Segel equation 1D
% $\partial_t u - D \Delta_x u + \chi \partial_x (u\partial_x c) = 0$
% $\partial_t c - D_2 \Delta_x c = \alpha u - \beta c$
close all; clear all; format long;

lim = 'MC'; % choose flux limiter

%% Model parameters
D = 1;
D2 = 1;
alpha =1; 
beta = 0.1;

chi = 100; % chemoattractant strength

%% Space grid
nbx = 200;                          % Number of points in the domain
L = 10;     
xm=0;
xM=10;% Length of the domain
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
u = abs(0.3+0.05*r);  % Random numbers
c = u;

% initialisation flux at boundary of cells
fluxm_h = zeros(nbx+1,1); 

figure(1)
plot(x,u, '-*');
hold on 
plot(x,c);
hold off;

u_up = u;
tt = dt; inc = 1; inc2=1;
%%
while(tt < T)
    tt
    % the velocity here is the gradient of c
    v=zeros(nbx+1,1);
    v(2:nbx) = (c(2:nbx)-c(1:nbx-1))./h;
    
    % compute flux
    fluxm_h = ospre_convection_flux(u,v,nbx,h,dt);
    flux_diff = zeros(nbx+1,1);
    flux_diff(2:nbx) = (u(2:nbx)-u(1:nbx-1))./h;
    % CFL condition: adaptive time step
    dt = 1e-5; 
    dt = min(dt, CFL.*h./max(abs(v)));
    % compute solution
    u = u - chi*dt*(fluxm_h(2:nbx+1)-fluxm_h(1:nbx))./h + dt*D*(flux_diff(2:nbx+1)-flux_diff(1:nbx))./h;
    c = c + D2*dt*(v(2:nbx+1)-v(1:nbx))./h +dt*alpha.*u - dt*beta .*c;
    
    
    % plot the solution during simulation 
    if mod(inc,1000)==0 
        figure(1)
        plot(x,u, '-*');
        
        
        figure(2)
        plot(x,c, '-*')
        
        mass(inc2) = h*sum(u);
        figure(3)
        plot(mass)
        inc2=inc2+1;
    end
    inc = inc +1;
    tt = tt+dt;
end