function [fluxm_h] = TVD_convection_flux_y(u,v,nbx,nby,h,dt,lim)
%CONVECTION_FLUX This function computes the convective fluxes in the y direction
% at the interface of the cells for a transient equation
% An example of handled flux is: du/dt + div(u v) = 0
% Inputs:
% - u: matrix of scalars nbx*nby
% - v: y-component of velocity (real or matrix of size nby*nbx+1)
% - nbx: number of cells in x-direction
% - nby: number of cells in y-direction
% - h: constant grid size.
% - lim: name of the flux limiter (see the end of this file)
% Outputs:
% - fluxm_h : array of fluxes at the North boundary (top interface of each cell)

fmip12 = zeros(nby+1,nbx); fluxm_h = zeros(nby+1,nbx);
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
        

% calculating flux at interior borders of cells
% first order upwind
fmip12(3:nby-1,:) = u(2:nby-2,:).*max(0,v(3:nby-1,:)) + u(3:nby-1,:).*min(0,v(3:nby-1,:)); 
% matrix that contains at each interface which node is upwind in the x-direction
big_one = ones(nby+1,nbx);
index_j(3:nby-1,:) = big_one(3:nby-1,:).*(3:nby-1)'-sign(v(3:nby-1,:)); 
% high-order reconstruction
% local Courant Number
nu_mip12 = zeros(nby+1,nbx);
nu_mip12(3:nby-1,:) = v(3:nby-1,:)*dt/h;

fmip12(3:nby-1,:) = fmip12(3:nby-1,:) + (1-abs(nu_mip12(3:nby-1,:))).*abs(v(3:nby-1,:))*0.5...
        .*(TVD_limiter(u(3:nby-1,:)-u(2:nby-2,:),u(index_j(3:nby-1,:))-u(index_j(3:nby-1,:)-1),lim));

fluxm_h(3:nby-1,:) = fmip12(3:nby-1,:);
% boundary conditions: zeros Neumann. Nothing to add because the flux
% is zero on the boundaries
ii =2;

fmip12(ii,:) = u(ii-1,:).*max(0,v(ii,:)) + u(ii,:).*min(0,v(ii,:)); 
index_j(ii,:)  = ii-sign(v(ii,:));
nu_mip12(ii,:) = v(ii,:)*dt/h;


index_good = find(index_j(ii,:)~=1);
index_bad = find(index_j(ii,:)==1);

index_vec = index_j(ii,index_good);
global_index_good = (index_good-1).*nby+index_vec;
index_vec = index_j(ii,index_bad);
global_index_bad = (index_bad-1).*nby+index_vec;

m_diff_test = zeros(1,nbx);
m_diff_test(index_good) =  u(global_index_good)-u(global_index_good-1);
m_diff_test(index_bad) = u(global_index_bad+1)-u(global_index_bad);

fmip12(ii,index_good) = fmip12(ii,index_good) + (1-abs(nu_mip12(ii,index_good))).*abs(v(ii,index_good))*0.5...
        .*(TVD_limiter(u(ii,index_good)-u(ii-1,index_good),...
        m_diff_test(index_good),lim));
    
fmip12(ii,index_bad) = fmip12(ii,index_bad) + (1-abs(nu_mip12(ii,index_bad))).*abs(v(ii,index_bad))*0.5...
        .*(TVD_limiter(u(ii,index_bad)-u(ii-1,index_bad),...
         m_diff_test(index_bad),lim));

fluxm_h(ii,:) = fmip12(ii,:);

ii =nbx;

fmip12(ii,:) = u(ii-1,:).*max(0,v(ii,:)) + u(ii,:).*min(0,v(ii,:)); 
index_j(ii,:)  = ii-sign(v(ii,:));
nu_mip12(ii,:) = v(ii,:)*dt/h;


index_good = find(index_j(ii,:)~=nbx+1);
index_bad = find(index_j(ii,:)==nbx+1);

index_vec = index_j(ii,index_good);
global_index_good = (index_good-1).*nby+index_vec;
index_vec = index_j(ii,index_bad);
global_index_bad = (index_bad-1).*nby+index_vec;

m_diff_test = zeros(1,nbx);
m_diff_test(index_good) =  u(global_index_good)-u(global_index_good-1);
m_diff_test(index_bad) = u(global_index_bad-1)-u(global_index_bad);

fmip12(ii,index_good) = fmip12(ii,index_good) + (1-abs(nu_mip12(ii,index_good))).*abs(v(ii,index_good))*0.5...
        .*(TVD_limiter(u(ii,index_good)-u(ii-1,index_good),...
        m_diff_test(index_good),lim));
    
fmip12(ii,index_bad) = fmip12(ii,index_bad) + (1-abs(nu_mip12(ii,index_bad))).*abs(v(ii,index_bad))*0.5...
        .*(TVD_limiter(u(ii,index_bad)-u(ii-1,index_bad),...
         m_diff_test(index_bad),lim));

fluxm_h(ii,:) = fmip12(ii,:);


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
            L = ones(size(r))*1.5;
            L(isfinite(r)) = 1.5.*(r(isfinite(r)).^2+r(isfinite(r)))./(r(isfinite(r)).^2+r(isfinite(r))+1);
            L(a==b)=0;
        otherwise
            error('Wrong Flux limiter name')
    end
end

end

