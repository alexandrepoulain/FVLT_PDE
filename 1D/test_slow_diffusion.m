% Unsteady slow diffusion equation 1D
% The flux is projected on a pieciwise constant basis
% Low-order; upwind + High-Order: Lax-Wendroff
% $\partial_t u - \partial_x (u\partial_x u) = 0$
close all; clear all; format long;

lim = 'MC'; % choose flux limiter
v=1;
%% Space grid
nbx = 200;                          % Number of points in the domain
L = 20;     
xm=-10;
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
u = zeros(nbx,1); % numerical solution
u = u_exact(x,0);
uex = u_exact(x,0);

% initialisation flux at boundary of cells
fluxm_h = zeros(nbx+1,1); 

figure(1)
plot(x,u, '-*');
hold on 
plot(x,uex);
hold off;

u_up = u;
tt = dt; inc = 1; inc2=1;
%%
while(tt < T)
    tt
    % the velocity here is minus the gradient of u
    v=zeros(nbx+1,1);
    v(2:nbx) = -(u(2:nbx)-u(1:nbx-1))./h;
    
    % compute flux
%     fluxm_h = TVD_convection_flux(u,v,nbx,h,dt,lim);
    fluxm_h = ospre_convection_flux(u,v,nbx,h,dt);
    flux_up = upwind_convection_flux(u_up,v,nbx,h,dt);
    % CFL condition: adaptive time step
    dt = 1e-5; 
    dt = min(dt, CFL.*h./max(abs(v)));
    % compute solution
    u = u - 2*dt*(fluxm_h(2:nbx+1)-fluxm_h(1:nbx))./h;
    u_up = u_up - 2*dt*(flux_up(2:nbx+1)-flux_up(1:nbx))./h;
    
    
    % plot the solution during simulation 
    if mod(inc,1000)==0 
        uex = u_exact(x,tt);

        figure(1)
        plot(x,u, '-*');
        hold on 
        plot(x,u_up, '-*')
        plot(x,uex);
        hold off;
        
        mass(inc2) = h*sum(uex);
        figure(2)
        plot(mass)
        inc2=inc2+1;
    end
    inc = inc +1;
    tt = tt+dt;
end

% exact solution
function [u_ex]=u_exact(x,t)
    m=2;
    k = (m+1)^(-1);
    u_ex = zeros(numel(x),1);
    % Barenblatt solution 
    for i=1:numel(x)
        u_ex(i) = (t+1)^(-k)*max(0,(1-k*(m-1)/(2*m)*abs(x(i))^2/(t+1)^(2*k))^(1/(m-1)));      
    end
%    % heavyside
%     for i=1:numel(x)
%         if (abs(x(i)-t-0.2) < 0.1)
%             u_ex(i) = 1;
%         end
%     end
end

function [L] = TVD_limiter(a,b)
    S = (sign(a)+sign(b))*0.5;
    %L = S*(min(abs(a),abs(b))); % minmod
    %L = S.*2*abs(a).*abs(b)./( abs(a)+abs(b) ) ; % Van Leer
    %L = S.*min(0.5.*abs(a+b) ,min(2*abs(a) ,2*abs(b))); % MC
    L = S.*max( min(2*abs(a),abs(b) ),min(abs(a),2*abs(b))); % Superbee
end