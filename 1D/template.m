% Template transient equation 1D

close all; clear all; format long;

lim = 'MC'; % choose flux limiter if you use TVD_convection_flux

%% Model parameters

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

%% Initial conditions for unknows here
% initial condition

% initialisation flux at boundary of cells (need to initialize array of size (nbx+1,1))
fluxm_h = zeros(nbx+1,1); 


tt = dt; inc = 1; inc2=1;
%% start temporal loop
while(tt < T)
    tt
    % Compute the velocity at each interface here 
    v=zeros(nbx+1,1);
    v(2:nbx) = ;
    
    % compute flux
    fluxm_h = upwind_convection_flux(u,v,nbx,h,dt); % Or use TVD_convection_flux for high order
    
    % CFL condition: adaptive time step
    dt = 1e-5; 
    dt = min(dt, CFL.*h./max(abs(v)));
    
    % compute solution (solve your equation here)
    
    
    
    % plot the solution during simulation 
    if mod(inc,1000)==0 
        
        inc2=inc2+1;
    end
    inc = inc +1;
    tt = tt+dt;
end