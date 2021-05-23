function [fluxm_h] = Upwind_convection_flux(u,v,nbx,h,dt)
%UPWIND_CONVECTION_FLUX This function computes the convective fluxes at the
% interface of the cells for a transient equation using an upwind approach
% An example of handled flux is: du/dt + div(u v) = 0
% Inputs:
% - u: array of scalars
% - v: velocity (real or array)
% - nbx: number of points (corresponding to the number of cells)
% - h: constant grid size.
% Outputs:
% - fluxm_h : array of fluxes at i+1/2

fmip12 = zeros(nbx+1,1); fluxm_h = zeros(nbx+1,1);
% test if velocity is a real or an array
if size(v) == 1 % if constant scalar velocity
    % then expand the coeff on each interfaces
    v = ones(nbx+1,1).*v;
else
    if size(v) == size(fmip12) % velocity is given at each interfaces
        ; % nothing to do
    else
        msg = 'Incorrect size for velocity array: I need an array with velocities at each interfaces of cells';
        error(msg)
    end
end
        

% calculating flux at interior borders of cells
% first order upwind
fmip12(3:nbx-1) = u(2:nbx-2).*max(0,v(3:nbx-1)) + u(3:nbx-1).*min(0,v(3:nbx-1)); 
fluxm_h(3:nbx-1) = fmip12(3:nbx-1);
% boundary conditions: zeros Neumann. Nothing to add because the flux
% is zero on the boundaries
ii =2;
fmip12(ii) = u(ii-1).*max(0,v(ii)) + u(ii).*min(0,v(ii)); 
fluxm_h(ii) = fmip12(ii);

ii =nbx;
fmip12(ii) = u(ii-1).*max(0,v(ii)) + u(ii).*min(0,v(ii)); 
fluxm_h(ii) = fmip12(ii);



end

