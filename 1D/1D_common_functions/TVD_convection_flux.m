function [fluxm_h] = TVD_convection_flux(u,v,nbx,h,dt,lim)
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
%index_j(3:nbx-1) = (3:nbx-1)'-sign(v(3:nbx-1));
index_j(3:nbx-1) = (3:nbx-1)'-1;
% high-order reconstruction
% local Courant Number
nu_mip12 = zeros(nbx+1,1);
nu_mip12(3:nbx-1) = v(3:nbx-1)*dt/h;

fmip12(3:nbx-1) =fmip12(3:nbx-1) + (1-abs(nu_mip12(3:nbx-1))).*abs(v(3:nbx-1))*0.5...
        .*(TVD_limiter(u(3:nbx-1)-u(2:nbx-2),u(index_j(3:nbx-1))-u(index_j(3:nbx-1)-1),lim));

fluxm_h(3:nbx-1) = fmip12(3:nbx-1);
% boundary conditions: zeros Neumann. Nothing to add because the flux
% is zero on the boundaries
ii =2;

fmip12(ii) = u(ii-1).*max(0,v(ii)) + u(ii).*min(0,v(ii)); 
jj  = ii-sign(v(ii));
nu_mip12(ii) = v(ii)*dt/h;
if jj~= 1
    fmip12(ii) = fmip12(ii) + (1-abs(nu_mip12(ii)))*abs(v(ii))*0.5...
        *(TVD_limiter(u(ii)-u(ii-1),u(jj)-u(jj-1),lim));
else
    fmip12(ii) = fmip12(ii) + (1-abs(nu_mip12(ii)))*abs(v(ii))*0.5...
        *(TVD_limiter(u(ii)-u(ii-1),u(ii)-u(ii-1),lim));
end

fluxm_h(ii) = fmip12(ii);

ii =nbx;

fmip12(ii) = u(ii-1).*max(0,v(ii)) + u(ii).*min(0,v(ii)); 
jj  = ii-sign(v(ii));
nu_mip12(ii) = v(ii)*dt/h;
if jj~=nbx+1
    fmip12(ii) = fmip12(ii) + (1-abs(nu_mip12(ii)))*abs(v(ii))*0.5...
        *(TVD_limiter(u(ii)-u(ii-1),u(jj)-u(jj-1),lim));
else
    fmip12(ii) = fmip12(ii) + (1-abs(nu_mip12(ii)))*abs(v(ii))*0.5...
        *(TVD_limiter(u(ii)-u(ii-1),u(ii)-u(ii-1),lim));
end


fluxm_h(ii) = fmip12(ii);


function [L] = TVD_limiter(a,b,lim)
    % lim is a string that gives the flux limiter
    S = (sign(a)+sign(b))*0.5;
    
    switch lim
        case 'minmod'
            L = S.*(min(abs(a),abs(b))); % minmod
        case 'VL'
            L = S.*2.*abs(a).*abs(b)./( abs(a)+abs(b) ) ; % Van Leer
            L(isfinite(L)~=1)=2;
            L(a==b)=0;
        case 'MC'
            L = S.*min(0.5.*abs(a+b) ,min(2*abs(a) ,2*abs(b))); % MC
        case 'SB'
            L = S.*max( min(2*abs(a),abs(b) ),min(abs(a),2*abs(b))); % Superbee
        case 'ospre'
            r = b./a;
            L = 1.5.*((r).^2+r)./(r.^2+r+1);
            L(isfinite(L)~=1)=1.5;
        otherwise
            error('Wrong Flux limiter name')
    end
end


end

