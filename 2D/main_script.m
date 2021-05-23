% 2D diffusion Finite volume scheme
% solve d_t u - D(d_x u + d_y u) = 0
clear all; close all; format long;

% flux limiter for TVD scheme
lim = 'VL'


%% Parameters of the model 
mu = 2; % mobilities
nu = 1;
CFL = 0.0001;
gam = 30; % gamma for pressure
Kg = (gam+1)/gam; % pressure coeff
Pm = 30; % max pressure
pres = @(m,n) Kg*(m+n).^gam; % pressure fun
Gp = @(p) 200/pi*atan(4*max(0,Pm-p)); % growth fun

% Space discretization
Lx = 45;
Ly = 45;

dx = 0.2; dy = 0.2;

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
am = 0.1;
an = 0.8;
bm =5e-1;
bn =5e-8;

x0=0; y0=0;

m = am*exp(-bm*((X-x0).^2+(Y-y0).^2));
n = an*exp(-bn*((X-x0).^2+(Y-y0).^2));

% h1 = figure(1)
% surf(X,Y,m);
% view(2);
% colorbar
% shading interp
% filename = ['Images2/m-', num2str(0), '.png'];
% saveas(h1,filename)

% h2 = figure(2)
% surf(X,Y,n);
% view(2);
% colorbar
% shading interp
filename = ['Images3/m-000', num2str(0), '.vtk'];
vtkwrite(filename, 'structured_grid', X, Y, Z,'scalars', 'm', m);


tt_plot=0;
nbx= size(x,2);
nby=size(y,2);
inc=1;inc2=1;
%%
while tt<T
    tt
    p = Kg*(m+n).^gam;
    
    v_x = zeros(nby,nbx+1); % velocity in x-direction
    v_x(:,2:nbx) = -(p(:,2:nbx)-p(:,1:nbx-1))./dx; 
    v_y = zeros(nby+1,nbx); % velocity in x-direction
    v_y(2:nby,:) = -(p(2:nby,:)-p(1:nby-1,:))./dy; 
    
    [Fmxi12] = ospre_convection_flux_x(m,v_x,nbx,nby,dx,dt);
    [Fnxi12] = ospre_convection_flux_x(n,v_x,nbx,nby,dx,dt);
    [Fmyi12] = ospre_convection_flux_y(m,v_y,nbx,nby,dx,dt);
    [Fnyi12] = ospre_convection_flux_y(n,v_y,nbx,nby,dx,dt);
    
    dt=1e-4;
    ui12= v_x;
    dt= min(dt,CFL*min(dx)/(max(nu,mu)*(max(max(ui12, [], 'all'),abs(min(ui12, [], 'all'))))));
    ui12= v_y;
    dt= min(dt,CFL*min(dx)/(max(nu,mu)*(max(max(ui12, [], 'all'),abs(min(ui12, [], 'all'))))));
    

    
    gp = Gp(p);
    m = -mu*dt*(Fmxi12(:,2:nbx+1)-Fmxi12(:,1:nbx))./(dy) ...
        -mu*dt*(Fmyi12(2:nby+1,:)-Fmyi12(1:nby,:))./(dy) + m ;
    m = m + dt*gp.*m;
    n = -nu*dt*(Fnxi12(:,2:nbx+1)-Fnxi12(:,1:nbx))./(dy) ...
        -nu*dt*(Fnyi12(2:nby+1,:)-Fnyi12(1:nby,:))./(dy) + n;
    
    
    if tt_plot>0.01
        tt_plot=0;
        h1 = figure(1)
        surf(X,Y,m);
        view(2);
        colorbar
        shading interp
%         filename = ['Images2/m-', num2str(inc2), '.png'];
%         saveas(h1,filename)

        h2 = figure(2)
        surf(X,Y,n);
        view(2);
        colorbar
        shading interp
%         filename = ['Images2/n-', num2str(inc2), '.png'];
%         saveas(h2,filename)
        
        massn(inc2) = dx.*dy.*sum(sum(n(:,:)));
        figure(3)
        plot(massn)
        
        if inc2 < 10
            filename = ['Images3/m-000', num2str(inc2), '.vtk'];
        else
            if inc2 <100
                filename = ['Images3/m-00', num2str(inc2), '.vtk'];
            else
                if inc2 < 1000
                    filename = ['Images3/m-0', num2str(inc2), '.vtk'];
                else
                    filename = ['Images3/m-', num2str(inc2), '.vtk'];
                end
            end
        end
        vtkwrite(filename, 'structured_grid', X, Y, Z,'scalars', 'm', m);
        
        if inc2 < 10
            filename = ['Images3/n-000', num2str(inc2), '.vtk'];
        else
            if inc2 <100
                filename = ['Images3/n-00', num2str(inc2), '.vtk'];
            else
                if inc2 < 1000
                    filename = ['Images3/n-0', num2str(inc2), '.vtk'];
                else
                    filename = ['Images3/n-', num2str(inc2), '.vtk'];
                end
            end
        end
        
        vtkwrite(filename, 'structured_grid', X, Y, Z,'scalars', 'n', n);
        
        
        inc2=inc2+1;
        
        
    end
    tt = dt+tt;
    tt_plot = tt_plot+dt;
    inc=inc+1;
end