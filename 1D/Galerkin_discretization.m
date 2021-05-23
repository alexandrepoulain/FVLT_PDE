% Unsteady Convection equation 1D
% The flux is projected on a pieciwise constant basis
% Low-order; upwind + High-Order: Lax-Wendroff
% $\partial_t u + \partial_x f(u) = 0$
close all; clear all; format long;

% two populations Scheme
gam = 30;
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
nb_pts = 101;                          % Number of points in the domain
L = 10;     
xm=0;
xM=10;% Length of the domain
h = (xM-xm)/(nb_pts-1);                     % grid space
nbx = (xM-xm)/h;                            % number of elements;
x = xm:h:xM;                            % discretization

%% Time discretization
T = 0.5;
dt = 1e-6;
t = 0:dt:T;  

 theta = 0; % Forward Euler, explicit first order
% % theta = 0.5; % Crank-Nicolson, implicit second order
%theta = 1; % Backward Euler, explicit first order

%% Initialisation
n = zeros(nbx+1,1); % non-growing population
m = zeros(nbx+1,1); % growing population
p = zeros(nbx+1,1); % pressure
% initial condition
sigx = 0.3;% spread of gaussian
x0=0; %sigx =1;
am =0.5;
x0 = 1;% pos of gaussian
m(:,1) = am.*exp(-bm.*(x-x0).^2/(2*sigx^2));
x0 = 4;% pos of gaussian
n(:,1) = am.*exp(-bm.*(x-x0).^2/(2*sigx^2));


% intn=2*exp(-0.1*(x+8).^2);
% intn(x<-8)=0;
% intn(intn<0.5)=1e-30;
% 
% intm=1*exp(-1*(x+10).^2);
% intm(intm<0.008)=0;
% m=intm'; %(initial condition)
% 
% n=intn'; %(initial condition)

%% Assembling
% Shape functions in reference space
psi = 1; % constant basis function


N1 = @(xi) 1/2*(1-xi);
N2 = @(xi) 1/2*(1+xi);
B=@(e) (1/h(e))*[-1. 1.]; % derivative of shape function in physical space
dN = 1/2.*[-1. 1.];       % derivative of shape fn in ref space
% Gauss points xi1 = -1/sqrt(3) xi2 =1/sqrt(3)
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
% weight for gauss points
w1 = 1;
w2 = 1;

%% Constant matrices
% assembling convec matrix
C = zeros(nbx+1, nbx+1);
M = zeros(nbx+1, nbx+1);
Mtilde = zeros(nbx+1, nbx+1);


celem = zeros(2,2);
melem = zeros(2,2);
melemtilde = zeros(2,2);


for e=1:nbx
    % Position of x in physical space
    e1 = x(e);
    e2 = x(e+1);
    Xe = [e1 e2];
    J = dN*Xe';
    global_dN = J\dN;       
    % submatrix
    melem = h/2*(w1*[N1(xi1)^2 N1(xi1)*N2(xi1) ; N1(xi1)*N2(xi1) N2(xi1)^2] + ...
        w2*[N1(xi2)^2 N1(xi2)*N2(xi2); N1(xi2)*N2(xi2) N2(xi2)^2]);
    % return to matrix M
    M(e:e+1,e:e+1) = M(e:e+1,e:e+1) + melem(:,:);
    
end

Mlump = zeros(nbx+1, nbx+1);
% Mass lump
for i=1:nbx+1
    Mlump(i,i) = sum(M(i,:));
end
phi=1;
Cm = zeros(nbx+1, nbx+1);
Cn = zeros(nbx+1, nbx+1);
for e=1:nbx
    % Position of x in physical space
    e1 = x(e);
    e2 = x(e+1);
    Xe = [e1 e2];
    J = dN*Xe';
    global_dN = J\dN;
    % submatrix
    NN1 = [phi,phi]; NN2 = [phi,phi];
    cmelem =abs(det(J))*(w1*NN1'*global_dN + w2*NN2'*global_dN);

    cnelem =abs(det(J))*(w1*NN1'*global_dN + w2*NN2'*global_dN);
    Cm(e:e+1,e:e+1) = Cm(e:e+1,e:e+1) + cmelem(:,:);
    Cn(e:e+1,e:e+1) = Cn(e:e+1,e:e+1) + cnelem(:,:);
end
fluxm_h = zeros(nbx+1,1); fluxn_h=zeros(nbx+1,1);
    
figure(1)
plot(x,m);
hold on 
plot(x,n);
hold off;
tt = dt; inc = 1; inc2=1;
while(tt < T)
    tt
    fm = zeros(numel(x),1);
    fn = zeros(numel(x),1);
    p = Kgam*(m+n).^gam;
    dnp = Kgam*gam*(m+n).^(gam-1);
    dmp=dnp;
    for ii= 1: numel(x)-1
        fm(ii) = m(ii)*(p(ii+1)-p(ii))/h; 
        fn(ii) = n(ii)*(p(ii+1)-p(ii))/h; 
    end
    
    % calculating flux
    for ii= 2:numel(x)-1
        
        % derivative of fluxes
        % First order upwind scheme
        
            vnmip12 = (n(ii)*(dmp(ii))+n(ii+1)*(dmp(ii+1)))*0.5*(m(ii+1)-m(ii))/h;
            vmmip12 = (p(ii+1)-p(ii))/h + (m(ii)*(dmp(ii))+m(ii+1)*(dmp(ii+1)))*0.5*(m(ii+1)-m(ii))/h;
        
            vnmim12 = (n(ii-1)*(dmp(ii-1))+n(ii)*(dmp(ii)))*0.5*(m(ii)-m(ii-1))/h ;
            vmmim12 = (p(ii)-p(ii-1))/h + (m(ii-1)*(dmp(ii-1))+m(ii)*(dmp(ii)))*0.5*(m(ii)-m(ii-1))/h;
        
            vnnip12 = (p(ii+1)-p(ii))/h + (n(ii)*dnp(ii)+n(ii+1)*dnp(ii+1))*0.5*(n(ii+1)-n(ii))/h;
            vmnip12 = (m(ii)*dnp(ii)+m(ii+1)*dnp(ii+1))*0.5*(n(ii+1)-n(ii))/h;
        
            vnnim12 = (p(ii)-p(ii-1))/h + (n(ii-1)*dnp(ii-1)+n(ii)*dnp(ii))*0.5*(n(ii)-n(ii-1))/h;
            vmnim12 = (m(ii-1)*dnp(ii-1)+m(ii)*dnp(ii))*0.5*(n(ii)-n(ii-1))/h;
        

        % fluxes at interface: first-order upwind
        if ((p(ii+1)-p(ii)) <= 0)
            noeud_up_p = ii-1;
            fmip12 = (m(ii))*((p(ii+1)-p(ii))/h) ;
            fnip12 = (n(ii))*((p(ii+1)-p(ii))/h);
        else
            noeud_up_p = ii+1;
            fmip12 = (m(ii+1))*((p(ii+1)-p(ii))/h) ;
            fnip12 = (n(ii+1))*((p(ii+1)-p(ii))/h);
        end
        if (p(ii)-p(ii-1)) <= 0
            noeud_up_m = ii-1;
            fmim12 = (m(ii-1))*((p(ii)-p(ii-1))/h) ;
            fnim12 = (n(ii-1))*((p(ii)-p(ii-1))/h);
        else
            noeud_up_m = ii+1;
            fmim12 = (m(ii))*((p(ii)-p(ii-1))/h) ;
            fnim12 = (n(ii))*((p(ii)-p(ii-1))/h);
        end
        if (noeud_up_p == 1)
            noeud_up_p=2;
        end
        if (noeud_up_m == 1)
            noeud_up_m=2;
        end
        % high-order reconstruction
        % local Courant Number
        nu_mip12 = (p(ii+1)-p(ii))/h*dt/h;
        nu_mim12 = fmim12*dt/h;
        nu_nip12 = fnip12*dt/h;
        nu_nim12 = fnim12*dt/h;
        
        fmip12 = fmip12 + (1-abs(nu_mip12))*abs((p(ii+1)-p(ii))/h)*0.5...
            *(TVD_limiter(m(ii+1)-m(ii),m(noeud_up_p)-m(noeud_up_p-1)));
        fmim12 = fmim12 + (1-abs(nu_mim12))*abs((p(ii)-p(ii-1))/h)*0.5...
            *(TVD_limiter(m(ii)-m(ii-1),m(noeud_up_m)-m(noeud_up_m-1)));
        
        fnip12 = fnip12 + (1-abs(nu_nip12))*abs((p(ii+1)-p(ii))/h)*0.5...
            *(TVD_limiter(n(ii+1)-n(ii),n(noeud_up_p)-n(noeud_up_p-1)));
        fnim12 = fnim12 + (1-abs(nu_nim12))*abs((p(ii)-p(ii-1))/h)*0.5...
            *(TVD_limiter(n(ii)-n(ii-1),n(noeud_up_m)-n(noeud_up_m-1)));
        
        % difference of fluxes at interface
        fluxm_h(ii) = fmip12 - fmim12;
        fluxn_h(ii) = fnip12 - fnim12;
    end
    
    % boundary conditions: zeros Neumann
    % left boundary
    ii=1;
    if ((p(ii+1)-p(ii)) <= 0)
        noeud_up_p = ii-1;
        fmip12 = (m(ii))*((p(ii+1)-p(ii))/h) ;
        fnip12 = (n(ii))*((p(ii+1)-p(ii))/h);
    else
        noeud_up_p = ii+1;
        fmip12 = (m(ii+1))*((p(ii+1)-p(ii))/h) ;
        fnip12 = (n(ii+1))*((p(ii+1)-p(ii))/h);
    end

    if (noeud_up_p == 1)
        noeud_up_p=2;
    end
    if (noeud_up_m == 1)
        noeud_up_m=2;
    end
    % high-order reconstruction
    % local Courant Number
    nu_mip12 = (p(ii+1)-p(ii))/h*dt/h;
    nu_nip12 = fnip12*dt/h;

    fmip12 = fmip12 + (1-abs(nu_mip12))*abs((p(ii+1)-p(ii))/h)*0.5...
        *(TVD_limiter(m(ii+1)-m(ii),m(noeud_up_p)-m(noeud_up_p-1)));
    fnip12 = fnip12 + (1-abs(nu_nip12))*abs((p(ii+1)-p(ii))/h)*0.5...
        *(TVD_limiter(n(ii+1)-n(ii),n(noeud_up_p)-n(noeud_up_p-1)));
        
    fluxm_h(ii) = fmip12;
    fluxn_h(ii) =  fnip12;
    % Right boundary
    ii = nbx+1;
    
    % fluxes at interface: first-order upwind

    if (p(ii)-p(ii-1)) <= 0
        noeud_up_m = ii-1;
        fmim12 = (m(ii-1))*((p(ii)-p(ii-1))/h) ;
        fnim12 = (n(ii-1))*((p(ii)-p(ii-1))/h);
    else
        noeud_up_m = ii+1;
        fmim12 = (m(ii))*((p(ii)-p(ii-1))/h) ;
        fnim12 = (n(ii))*((p(ii)-p(ii-1))/h);
    end
    if (noeud_up_p == 1)
        noeud_up_p=2;
    end
    if (noeud_up_m == 1)
        noeud_up_m=2;
    end
    % high-order reconstruction
    % local Courant Number
    nu_mim12 = fmim12*dt/h;
    nu_nim12 = fnim12*dt/h;


    fmim12 = fmim12 + (1-abs(nu_mim12))*abs((p(ii)-p(ii-1))/h)*0.5...
        *(TVD_limiter(m(ii)-m(ii-1),m(noeud_up_m)-m(noeud_up_m-1)));


    fnim12 = fnim12 + (1-abs(nu_nim12))*abs((p(ii)-p(ii-1))/h)*0.5...
        *(TVD_limiter(n(ii)-n(ii-1),n(noeud_up_m)-n(noeud_up_m-1)));
    
    fluxm_h(ii) = -fmim12;
    fluxn_h(ii) = -fnim12;
    
    
    Gp = rate*(1/(pi)*atan(4*max(0,Pm-p))).*m; 
    m = (Mlump./dt)^(-1)*(Mlump./dt*m + nu1.*fluxm_h + Mlump*(Gp));
    n = (Mlump./dt)^(-1)*(Mlump./dt*n + nu2.*fluxn_h);
    if mod(inc,1000)==0 
        figure(1)
        plot(x,m);
        hold on 
        plot(x,n);
        hold off;
        
        figure(2)
        mass(inc2) = h*sum(n);
        plot(mass)
        inc2=inc2+1;
    end
    inc = inc +1;
    tt = tt+dt;
end

function [L] = TVD_limiter(a,b)
    S = (sign(a)+sign(b))*0.5;
    L = S*(min(abs(a),abs(b))); % minmod
end