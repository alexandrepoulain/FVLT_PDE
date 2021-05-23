% 2D diffusion Finite volume scheme
% solve d_t u - D(d_x u + d_y u) = 0
clear all; close all; format long;

% flux limiter for TVD scheme
lim = 'VL'


%% Parameters of the model 
mu = 1; % mobilities
nu = 1;
CFL = 0.0005;
gam = 30; % gamma for pressure
Kg = (gam+1)/gam; % pressure coeff
Pm = 30; % max pressure
pres = @(m,n) Kg*(m+n).^gam; % pressure fun
Gp = @(p) 200/pi*atan(4*max(0,Pm-p)); % growth fun

% Space discretization
Lx = 10;
Ly = 10;

dx = 0.1; dy = 0.1;

x = dx/2:dx:Lx;
y = dy/2:dy:Ly;

[X,Y] = meshgrid(x,y);
left_boundary = 0;
right_boundary = Lx;
bottom_boundary =0;
top_boundary=Ly;
% Time discretization 
dt = 1e-6;
T = 1;
tt=0;

Z = zeros(size(X));


% initial condition
am = 1;
an = 0.8;
bm =5e-1;
bn =5e-8;

sig = 0.5;
x0=2; y0=2;
init_f = @(x,y,t) ((sig/sqrt(4*t+sig^2))^2*(exp(-((x-x0)/sqrt(4*t+sig^2)).^2)).*(exp(-((y-y0)/sqrt(4*t+sig^2)).^2))) ;
u_init = init_f(X,Y,0);


m = u_init;

h1 = figure(1)
surf(X,Y,m);
view(2);
colorbar
shading interp


tt_plot=0;
nbx= size(x,2);
nby=size(y,2);
inc=1;inc2=1;
%%
while tt<T
    tt
    
    v_x = zeros(nby,nbx+1); % velocity in x-direction
    v_x(:,2:nbx) = -(m(:,2:nbx)-m(:,1:nbx-1))./dx; 
    v_y = zeros(nby+1,nbx); % velocity in x-direction
    v_y(2:nby,:) = -(m(2:nby,:)-m(1:nby-1,:))./dy; 
    
%     [Fmxi12] = TVD_convection_flux_x(ones(nbx,nby),v_x,nbx,nby,dx,dt,lim);
%     [Fmyi12] = TVD_convection_flux_y(ones(nbx,nby),v_y,nbx,nby,dx,dt,lim);
    
    [Fmxi12] = ospre_convection_flux_x(ones(nbx,nby),v_x,nbx,nby,dx,dt);
    [Fmyi12] = ospre_convection_flux_y(ones(nbx,nby),v_y,nbx,nby,dx,dt);
    
    dt=1e-4;
    ui12= v_x;
    dt= min(dt,CFL*min(dx)/(max(nu,mu)*(max(max(ui12, [], 'all'),abs(min(ui12, [], 'all'))))));
    ui12= v_y;
    dt= min(dt,CFL*min(dx)/(max(nu,mu)*(max(max(ui12, [], 'all'),abs(min(ui12, [], 'all'))))));
    

        m = -mu*dt*(Fmxi12(:,2:nbx+1)-Fmxi12(:,1:nbx))./(dy) ...
        -mu*dt*(Fmyi12(2:nby+1,:)-Fmyi12(1:nby,:))./(dy) + m ;

    
        u_ex = init_f(X,Y,tt);

    if tt_plot>0.01
        tt_plot=0;
        h1 = figure(1)
        surf(X,Y,m);
        view(2);
        colorbar
        shading interp
        
        figure(2)
        surf(X,Y,u_ex);
        view(2);
        colorbar
        shading interp

        error(inc2) = dx.*dy.*sqrt(sum(sum((m-u_ex).^2)));
        
        
        massn(inc2) = dx.*dy.*sum(sum(m(:,:)));
        figure(3)
        plot(massn)
      
        figure(4)
        plot(error)
        
        inc2=inc2+1;
        
        
    end
    tt = dt+tt;
    tt_plot = tt_plot+dt;
    inc=inc+1;
end