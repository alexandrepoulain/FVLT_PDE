% Unsteady Convection equation 1D
% The flux is projected on a pieciwise constant basis
% Low-order; upwind + High-Order: Lax-Wendroff
% $\partial_t u + \partial_x f(u) = 0$
close all; clear all; format long;

lim = 'minmod';
% two populations Scheme
gam = 130;
Kgam = (gam+1)/gam;

am =0.1;
an =0.8;
bm =0.5;
bn =5e-7;
Pm =  30;
rate=20;

nu1 = 2; % mobility 1
nu2 = 1; % mobility 2
%Pe=1;
%% Space grid
nbx = 100;                          % Number of points in the domain
L = 10;     
xm=0;
xM=10;% Length of the domain
h = (xM-xm)/(nbx-1);                     % grid space
x = xm:h:xM;                            % discretization

%% Time discretization
T = 0.5;
dt = 1e-6;
t = 0:dt:T;  
CFL=0.0001;

%% Initialisation
n = zeros(nbx,1); % non-growing population
m = zeros(nbx,1); % growing population
p = zeros(nbx,1); % pressure
% initial condition
sigx = 0.3;% spread of gaussian
x0=0; %sigx =1;
am =0.5;
x0 = 1;% pos of gaussian
m(:,1) = am.*exp(-bm.*(x-x0).^2/(2*sigx^2));
x0 = 4;% pos of gaussian
n(:,1) = am.*exp(-bm.*(x-x0).^2/(2*sigx^2));
% 
% intn=2*exp(-0.1*(x+5).^2);
% intn(x<-5)=0;
% intn(intn<0.5)=1e-30;
% 
% intm=1*exp(-1*(x+10).^2);
% intm(intm<0.008)=0;
% m=intm'; %(initial condition)
% 
% n=intn'; %(initial condition)


% initialisation flux at boundary of cells
fluxm_h = zeros(nbx+1,1); 
fluxn_h=zeros(nbx+1,1); 
    
figure(1)
plot(x,m);
hold on 
plot(x,n);
hold off;

m_up=m; n_up=n;

tt = dt; inc = 1; inc2=1;tt_plot=0;
%%
while(tt < T)
    tt
    fmip12 = zeros(nbx+1,1); fluxm_h = zeros(nbx+1,1);
    fnip12 = zeros(nbx+1,1); fluxn_h = zeros(nbx+1,1);

    % pressure
    p = Kgam*(m+n).^gam;
    % compute the velocity at each interfaces
    v = zeros(nbx+1,1);
    v(2:nbx) = -(p(2:nbx)-p(1:nbx-1))./h;
    
%     fluxm_h = TVD_convection_flux(m,v,nbx,h,dt,lim);
%     fluxn_h = TVD_convection_flux(n,v,nbx,h,dt,lim);
    fluxm_h = ospre_convection_flux(m,v,nbx,h,dt);
    fluxn_h = ospre_convection_flux(n,v,nbx,h,dt);
    
    
    % CFL condition: adaptive time step
    dt = 1e-6;
    grad_pres = (p(2:nbx)-p(1:nbx-1))./h;
    grad_m = abs(m(2:nbx)-m(1:nbx-1))./h;
    grad_n = abs(n(2:nbx)-n(1:nbx-1))./h;
    dt = min(dt, CFL.*h./(max(nu1,nu2)*max(max(abs( grad_pres )),max( max(grad_n), max(grad_m))) ) );
    
    
    Gp = rate*(1/(pi)*atan(4*max(0,Pm-p))).*m; 
    m = m - dt*nu1.*(fluxm_h(2:nbx+1)-fluxm_h(1:nbx))./h + dt*Gp;
    n = n - dt*nu2.*(fluxn_h(2:nbx+1)-fluxn_h(1:nbx))./h;
    
    % pressure
    p = Kgam*(m_up+n_up).^gam;
    % compute the velocity at each interfaces
    v = zeros(nbx+1,1);
    v(2:nbx) = -(p(2:nbx)-p(1:nbx-1))./h;
    
    fluxm_h = upwind_convection_flux(m_up,v,nbx,h,dt);
    fluxn_h = upwind_convection_flux(n_up,v,nbx,h,dt);
    
    
    % CFL condition: adaptive time step
%     dt = 1e-6;
%     grad_pres = (p(2:nbx)-p(1:nbx-1))./h;
% 
%     dt = min(dt, CFL.*h./(max(nu1,nu2)*max(max(abs( grad_pres )) ) ));
    
    
    Gp = rate*(1/(pi)*atan(4*max(0,Pm-p))).*m_up; 
    m_up = m_up - dt*nu1.*(fluxm_h(2:nbx+1)-fluxm_h(1:nbx))./h + dt*Gp;
    n_up = n_up - dt*nu2.*(fluxn_h(2:nbx+1)-fluxn_h(1:nbx))./h;
    
    
    
    
    if tt_plot >0.01
        tt_plot=0;
        figure(1)
        plot(x,m);
        hold on 
        plot(x,n);
        plot(x,m_up, '-*');
        plot(x,n_up, '-*');
        hold off;
        
        figure(3)
        plot(x,p);
        
        figure(2)
        mass(inc2) = h*sum(n);
        plot(mass)
        inc2=inc2+1;
    end
    inc = inc +1;
    tt = tt+dt;
    tt_plot=tt_plot+dt;
end

function [L] = TVD_limiter(a,b)
    S = (sign(a)+sign(b))*0.5;
    %L = S*(min(abs(a),abs(b))); % minmod
    %L = S.*2*abs(a).*abs(b)./( abs(a)+abs(b) ) ; % Van Leer
    %L = S.*min(0.5.*abs(a+b) ,min(2*abs(a) ,2*abs(b))); % MC
    L = S.*max( min(2*abs(a),abs(b) ),min(abs(a),2*abs(b))); % Superbee
end