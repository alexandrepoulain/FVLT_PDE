function [fluxm_h,fluxn_h] = flux_x(m,n,p,dt,h,nbx,nby,Kgam,gam)
%FLUX_X Compute fluxes across the x boundary

fmip12 = zeros(nby,nbx+1); fnip12 = zeros(nby,nbx+1);
fluxm_h = zeros(nby,nbx+1); fluxn_h=zeros(nby,nbx+1);

fmip12(:,3:nbx-1) = (m(:,3:nbx-1)).*max((p(:,3:nbx-1)-p(:,2:nbx-2))/h,0)+...
         (m(:,2:nbx-2)).*min((p(:,3:nbx-1)-p(:,2:nbx-2))/h,0);
fnip12(:,3:nbx-1) = (n(:,3:nbx-1)).*max((p(:,3:nbx-1)-p(:,2:nbx-2))/h,0)+...
     (n(:,2:nbx-2)).*min((p(:,3:nbx-1)-p(:,2:nbx-2))/h,0);

nu_mip12 = zeros(nby,nbx+1); nu_nip12= zeros(nby,nbx+1);
nu_mip12(:,3:nbx-1) = fmip12(:,3:nbx-1)*dt/h;
nu_nip12(:,3:nbx-1) = fnip12(:,3:nbx-1)*dt/h;

% Coumpute index of upwind nodes
index_j = zeros(nby,nbx+1);
index_j(:,3:nbx-1) = ones(nby,nbx-3).*(3:nbx-1)-sign((p(:,3:nbx-1)-p(:,2:nbx-2))/h);

fluxm_h(:,3:nbx-1) = fmip12(:,3:nbx-1) + (1-abs(nu_mip12(:,3:nbx-1))).*abs((p(:,3:nbx-1)-p(:,2:nbx-2))/h)*0.5...
        .*(TVD_limiter(m(:,3:nbx-1)-m(:,2:nbx-2),m(index_j(:,3:nbx-1))-m(index_j(:,3:nbx-1)-1)));

fluxn_h(:,3:nbx-1) = fnip12(:,3:nbx-1) + (1-abs(nu_nip12(:,3:nbx-1))).*abs((p(:,3:nbx-1)-p(:,2:nbx-2))/h)*0.5...
    .*(TVD_limiter(n(:,3:nbx-1)-n(:,2:nbx-2),n(index_j(:,3:nbx-1))-n(index_j(:,3:nbx-1)-1)));



% boundary conditions: zeros Neumann. Nothing to add because the flux
% is zero on the boundaries
ii =2;


fmip12(:,ii) = (m(:,ii)).*max((p(:,ii)-p(:,ii-1))/h,0)+...
         (m(:,ii-1)).*min((p(:,ii)-p(:,ii-1))/h,0);
fnip12(:,ii) = (n(:,ii)).*max((p(:,ii)-p(:,ii-1))/h,0)+...
     (n(:,ii-1)).*min((p(:,ii)-p(:,ii-1))/h,0);


nu_mip12(:,ii) = fmip12(:,ii)*dt/h;
nu_nip12(:,ii) = fnip12(:,ii)*dt/h;

index_j(:,ii) = ii.*ones(nbx,1)-sign((p(:,ii)-p(:,ii-1))/h);

test = zeros(nby,1);
    
index_correct = index_j( index_j(:,ii)~=1 ,ii);
index_not_correct = index_j( index_j(:,ii)==1 ,ii);


    fmip12(:,ii)  = fmip12(:,ii) + (1-abs(nu_mip12(:,ii))).*abs((p(:,ii)-p(:,ii-1))/h)*0.5...
            .*(TVD_limiter(m(:,ii)-m(:,ii-1),  ...
            (m(:, )-m(:,index_correct-1) ) ...
            + (m(:,index_not_correct+1)-m(:,index_not_correct ) )        ));
    fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(:,ii)-p(:,ii-1))/h)*0.5...
        *(TVD_limiter(n(ii)-n(ii-1),n(jj)-n(:,jj-1)));
    fluxm_h(ii) = fmip12;
    fluxn_h(ii) = fnip12;

    fmip12 = fmip12 + (1-abs(nu_mip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
            *(TVD_limiter(m(ii)-m(ii-1),m(ii)-m(ii-1)));
    fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
        *(TVD_limiter(n(ii)-n(ii-1),n(ii)-n(ii-1)));
    fluxm_h(ii) = fmip12;
    fluxn_h(ii) = fnip12;




ii =nbx;
if ((p(ii)-p(ii-1)) <= 0)
    jj = ii-1;
    fmip12 = (m(ii-1))*((p(ii)-p(ii-1))/h) ;
    fnip12 = (n(ii-1))*((p(ii)-p(ii-1))/h);
else
    jj= ii+1;
    fmip12 = (m(ii))*((p(ii)-p(ii-1))/h) ;
    fnip12 = (n(ii))*((p(ii)-p(ii-1))/h);
end
nu_mip12(ii) = fmip12*dt/h;
nu_nip12(ii) = fnip12*dt/h;
if jj~=nbx+1
    fmip12 = fmip12 + (1-abs(nu_mip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
            *(TVD_limiter(m(ii)-m(ii-1),m(jj)-m(jj-1)));
    fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
        *(TVD_limiter(n(ii)-n(ii-1),n(jj)-n(jj-1)));
    fluxm_h(ii) = fmip12;
    fluxn_h(ii) = fnip12;
else
    fmip12 = fmip12 + (1-abs(nu_mip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
            *(TVD_limiter(m(ii)-m(ii-1),m(ii)-m(ii-1)));
    fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
        *(TVD_limiter(n(ii)-n(ii-1),n(ii)-n(ii-1)));
    fluxm_h(ii) = fmip12;
    fluxn_h(ii) = fnip12;
end


function [L] = TVD_limiter(a,b)
    S = (sign(a)+sign(b))*0.5;
    %L = S*(min(abs(a),abs(b))); % minmod
    %L = S.*2*abs(a).*abs(b)./( abs(a)+abs(b) ) ; % Van Leer
    %L = S.*min(0.5.*abs(a+b) ,min(2*abs(a) ,2*abs(b))); % MC
    L = S.*max( min(2*abs(a),abs(b) ),min(abs(a),2*abs(b))); % Superbee
end

end

