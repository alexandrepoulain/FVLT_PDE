% 2D pure convection Finite volume scheme
% $\partial_t u + \partial_x f(u) = 0$
clear all; close all; format long;

% flux limiter for TVD scheme
lim = 'MC'


%% Parameters of the model 
mu = 2; % mobilities
nu = 1;
CFL = 0.0005;
gam = 30; % gamma for pressure
Kg = (gam+1)/gam; % pressure coeff
Pm = 30; % max pressure
pres = @(m,n) Kg*(m+n).^gam; % pressure fun
Gp = @(p) 200/pi*atan(4*max(0,Pm-p)); % growth fun

% Space discretization
Lx = 1;
Ly = 1;

dx = 0.01; dy = 0.01;

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

x0=0; y0=0;

u = u_exact(x,y,0);

h1 = figure(1)
surf(X,Y,u);
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
    v_x = 1; % velocity in x-direction
    v_y = 1; % velocity in x-direction
    
    [Fmxi12] = TVD_convection_flux_x(u,v_x,nbx,nby,dx,dt,lim);
    [Fmyi12] = TVD_convection_flux_y(u,v_y,nbx,nby,dx,dt,lim);
    
    dt=1e-4;
    dt= min(dt,CFL*min(dx)/((max(max(v_x, [], 'all')))));

  
    u = -dt*(Fmxi12(:,2:nbx+1)-Fmxi12(:,1:nbx))./(dy) ...
        -dt*(Fmyi12(2:nby+1,:)-Fmyi12(1:nby,:))./(dy) + u ;

    u_ex = u_exact(x,y,tt);
    
    if tt_plot>0.01
        tt_plot=0;
        h1 = figure(1)
        surf(X,Y,u);
        view(2);
        colorbar
        shading interp
        
        h2 = figure(2)
        surf(X,Y,u_ex);
        view(2);
        colorbar
        shading interp
        
        massn(inc2) = dx.*dy.*sum(sum(u(:,:)));
        figure(3)
        plot(massn)
        

        
        
        inc2=inc2+1;
        
        
    end
    tt = dt+tt;
    tt_plot = tt_plot+dt;
    inc=inc+1;
end

% exact solution
function [u_ex]=u_exact(x,y,t)
    u_ex = zeros(numel(x),numel(x));
%     % cosinus
%     for j = 1:numel(y)
%         for i=1:numel(x)
%             if (abs(x(i)+y(j)-t-0.2) < 0.12)
%                 u_ex(j,i) = 0.5*(1+cos(pi.*(x(i)+y(j)-t-0.2)/0.12));
%             end
%         end
%     end
   % heavyside Travel to the right
%    for j = 1:numel(y)
%         for i=1:numel(x)
%             if (abs(x(i)-t-0.2) < 0.1)
%                 u_ex(j,i) = 1;
%             end
%         end
%    end
%    heavyside Travel to the right and top
%    for j = 1:numel(y)
%         for i=1:numel(x)
%             if (abs(x(i)+y(j)-t-0.2) < 0.1)
%                 u_ex(j,i) = 1;
%             end
%         end
%    end
end