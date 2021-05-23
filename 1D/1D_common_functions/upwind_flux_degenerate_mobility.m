function [fluxm_h] = upwind_flux_degenerate_mobility(u,v,nbx,h,dt,deg_1,deg_2,BC)
%UPWIND_FLUX_DEGENERATE_MOBILITY This function computes the convective fluxes at the
% interface of the cells for a transient equation using an upwind approach.
% The mobility degenerates at two different points.
% Inputs:
% - u: array of scalars
% - v: velocity (real or array)
% - nbx: number of points (corresponding to the number of cells)
% - h: constant grid size.
% - deg_1: first part of degenerate mobility (degeneracy in 0)
% - deg_2: second part: degenerates at another point
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
fmip12(3:nbx-1) = deg_1(u(2:nbx-2)).*deg_2(u(3:nbx-1)).*max(0,v(3:nbx-1))...
    + deg_1(u(3:nbx-1)).*deg_2(u(2:nbx-2)).*min(0,v(3:nbx-1)); 
fluxm_h(3:nbx-1) = fmip12(3:nbx-1);

% boundary conditions:
switch BC
    case 'Homo-Neumann'
        %zeros Neumann. Nothing to add because the flux
        % is zero on the boundaries
        ii =2;
        fmip12(ii) = deg_1(u(ii-1)).*deg_2(u(ii)).*max(0,v(ii)) ...
            + deg_1(u(ii)).*deg_2(u(ii-1)).*min(0,v(ii)); 
        fluxm_h(ii) = fmip12(ii);

        ii =nbx;
        fmip12(ii) = deg_1(u(ii-1)).*deg_2(u(ii)).*max(0,v(ii))...
            + deg_1(u(ii)).*deg_2(u(ii-1)).*min(0,v(ii)); 
        fluxm_h(ii) = fmip12(ii);
    
    otherwise
        error('Boundary conditions not implemented yet');
end


end

