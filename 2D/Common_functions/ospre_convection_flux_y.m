function [fluxm_h] = ospre_convection_flux(u,v,nbx,nby,dx,dt)
%CONVECTION_FLUX This function computes the convective fluxes at the
% interface of the cells for a transient equation
% An example of handled flux is: du/dt + div(u v) = 0
% Inputs:
% - u: array of scalars
% - v: velocity (real or array)
% - nbx: number of points (corresponding to the number of cells)
% - h: constant grid size.
% - lim: name of the flux limiter (see the end of this file)
% Outputs:
% - fluxm_h : array of fluxes at i+1/2

fmip12 = zeros(nbx+1,nby); fluxm_h = zeros(nbx+1,nby);
% test if velocity is a real or an array
if size(v) == 1 % if constant scalar velocity
    % then expand the coeff on each interfaces
    v = ones(nby+1,nbx).*v;
else
    if size(v) == size(fmip12) % velocity is given at each interfaces
        ; % nothing to do
    else
        msg = 'Incorrect size for velocity array: I need an array with velocities at each interfaces of cells';
        error(msg)
    end
end

% Construct first order approximation of u inside the cells.
r = zeros(nbx,nby); u_l12 = zeros(nbx+1,nby); u_r12 = zeros(nbx+1,nby);
r(2:nbx-1,:) = (u(2:nbx-1,:)-u(1:nbx-2,:))./(u(3:nbx,:)-u(2:nbx-1,:));
r(1,:)=1./(u(2,:)-u(1,:)); r(nbx,:)= u(nbx,:)-u(nbx-1,:);
u_l12(2:nbx,:) = u(1:nbx-1,:) +0.5*phi(r(1:nbx-1,:)).*(u(2:nbx,:)-u(1:nbx-1,:)); 
u_r12(2:nbx-1,:) = u(2:nbx-1,:) -0.5*phi(r(2:nbx-1,:)).*(u(3:nbx,:)-u(2:nbx-1,:)); 
u_r12(1,:) = u(1,:) -0.5*phi(r(1,:)).*(u(2,:)-u(1,:));


% calculating flux at interior borders of cells
% first order upwind
fmip12(2:nbx,:) = u_l12(2:nbx,:).*max(0,v(2:nbx,:)) + u_r12(2:nbx,:).*min(0,v(2:nbx,:)); 

fluxm_h(2:nbx,:) = fmip12(2:nbx,:);


function [L] = phi(r)
    % lim is a string that gives the flux limiter
    L = 1.5.*((r).^2+r)./(r.^2+r+1);
    L(isfinite(L)~=1)=1.5;
end


end

