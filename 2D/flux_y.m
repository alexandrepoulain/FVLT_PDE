function [fluxm_h,fluxn_h] = flux_y(m,n,p,dt,h,nbx,nby,Kgam,gam)
%FLUX_X Compute fluxes across the x boundary


fluxm_h = zeros(nby,nbx); fluxn_h=zeros(nby,nbx);

% derivatives of p with respect to m and n
dnp = Kgam*gam*(m+n).^(gam-1);
dmp=dnp;

% calculating flux
for ii=1:nbx
    for jj= 2:nbx-1

        % derivative of fluxes
        % First order upwind scheme

            vnmip12 = (n(jj,ii)*(dmp(jj,ii))+n(jj+1,ii)*(dmp(jj+1,ii)))*0.5*(m(jj+1,ii)-m(jj,ii))/h;
            vmmip12 = (p(jj+1,ii)-p(jj,ii))/h + (m(jj,ii)*(dmp(jj,ii))+m(jj+1,ii)*(dmp(jj+1,ii)))*0.5*(m(jj+1,ii)-m(jj,ii))/h;

            vnmim12 = (n(jj-1,ii)*(dmp(jj-1,ii))+n(jj,ii)*(dmp(jj,ii)))*0.5*(m(jj,ii)-m(jj-1,ii))/h ;
            vmmim12 = (p(jj,ii)-p(jj-1,ii))/h + (m(jj-1,ii)*(dmp(jj-1,ii))+m(jj,ii)*(dmp(jj,ii)))*0.5*(m(jj,ii)-m(jj-1,ii))/h;

            vnnip12 = (p(jj+1,ii)-p(jj,ii))/h + (n(jj,ii)*dnp(jj,ii)+n(jj+1,ii)*dnp(jj+1,ii))*0.5*(n(jj+1,ii)-n(jj,ii))/h;
            vmnip12 = (m(jj,ii)*dnp(jj,ii)+m(jj+1,ii)*dnp(jj+1,ii))*0.5*(n(jj+1,ii)-n(jj,ii))/h;

            vnnim12 = (p(jj,ii)-p(jj-1,ii))/h + (n(jj-1,ii)*dnp(jj-1,ii)+n(jj,ii)*dnp(jj,ii))*0.5*(n(jj,ii)-n(jj-1,ii))/h;
            vmnim12 = (m(jj-1,ii)*dnp(jj-1,ii)+m(jj,ii)*dnp(jj,ii))*0.5*(n(jj,ii)-n(jj-1,ii))/h;


        % fluxes at interface: first-order upwind
        if ((p(jj+1,ii)-p(jj,ii)) <= 0)
            noeud_up_p = jj-1;
            fmip12 = (m(jj,ii))*((p(jj+1,ii)-p(jj,ii))/h) ;
            fnip12 = (n(jj,ii))*((p(jj+1,ii)-p(jj,ii))/h);
        else
            noeud_up_p = jj+1;
            fmip12 = (m(jj+1,ii))*((p(jj+1,ii)-p(jj,ii))/h) ;
            fnip12 = (n(jj+1,ii))*((p(jj+1,ii)-p(jj,ii))/h);
        end
        if (p(jj,ii)-p(jj-1,ii)) <= 0
            noeud_up_m = jj-1;
            fmim12 = (m(jj-1,ii))*((p(jj,ii)-p(jj-1,ii))/h) ;
            fnim12 = (n(jj-1,ii))*((p(jj,ii)-p(jj-1,ii))/h);
        else
            noeud_up_m = jj+1;
            fmim12 = (m(jj,ii))*((p(jj,ii)-p(jj-1,ii))/h) ;
            fnim12 = (n(jj,ii))*((p(jj,ii)-p(jj-1,ii))/h);
        end
        if (noeud_up_p == 1)
            noeud_up_p=2;
        end
        if (noeud_up_m == 1)
            noeud_up_m=2;
        end
        % high-order reconstruction
        % local Courant Number
        nu_mip12 = (p(jj+1,ii)-p(jj,ii))/h*dt/h;
        nu_mim12 = fmim12*dt/h;
        nu_nip12 = fnip12*dt/h;
        nu_nim12 = fnim12*dt/h;

        fmip12 = fmip12 + (1-abs(nu_mip12))*abs((p(jj+1,ii)-p(jj,ii))/h)*0.5...
            *(TVD_limiter(m(jj+1,ii)-m(jj,ii),m(noeud_up_p,ii)-m(noeud_up_p-1,ii)));
        fmim12 = fmim12 + (1-abs(nu_mim12))*abs((p(jj,ii)-p(jj-1,ii))/h)*0.5...
            *(TVD_limiter(m(jj,ii)-m(jj-1,ii),m(noeud_up_m,ii)-m(noeud_up_m-1,ii)));

        fnip12 = fnip12 + (1-abs(nu_nip12))*abs((p(jj+1,ii)-p(jj,ii))/h)*0.5...
            *(TVD_limiter(n(jj+1,ii)-n(jj,ii),n(noeud_up_p,ii)-n(noeud_up_p-1,ii)));
        fnim12 = fnim12 + (1-abs(nu_nim12))*abs((p(jj,ii)-p(jj-1,ii))/h)*0.5...
            *(TVD_limiter(n(jj,ii)-n(jj-1,ii),n(noeud_up_m,ii)-n(noeud_up_m-1,ii)));

        % difference of fluxes across interface
        fluxm_h(jj,ii) = fmip12 - fmim12;
        fluxn_h(jj,ii) = fnip12 - fnim12;
    end


    % boundary conditions: zeros Neumann
    % left boundary
    jj=1;
    if ((p(jj+1,ii)-p(jj,ii)) <= 0)
        fmip12 = (m(jj,ii))*((p(jj+1,ii)-p(jj,ii))/h) ;
        fnip12 = (n(jj,ii))*((p(jj+1,ii)-p(jj,ii))/h);
    else
        fmip12 = (m(jj+1,ii))*((p(jj+1,ii)-p(jj,ii))/h) ;
        fnip12 = (n(jj+1,ii))*((p(jj+1,ii)-p(jj,ii))/h);
    end
    fluxm_h(jj,ii) = fmip12;
    fluxn_h(jj,ii) =  fnip12;
    % Right boundary
    jj = nbx;

    if (p(jj,ii)-p(jj-1,ii)) <= 0

        fmim12 = (m(jj-1,ii))*((p(jj,ii)-p(jj-1,ii))/h) ;
        fnim12 = (n(jj-1,ii))*((p(jj,ii)-p(jj-1,ii))/h);
    else
        fmim12 = (m(jj,ii))*((p(jj,ii)-p(jj-1,ii))/h) ;
        fnim12 = (n(jj,ii))*((p(jj,ii)-p(jj-1,ii))/h);
    end

    fluxm_h(jj,ii) = -fmim12;
    fluxn_h(jj,ii) = -fnim12;
end


function [L] = TVD_limiter(a,b)
    S = (sign(a)+sign(b))*0.5;
    L = S*(min(abs(a),abs(b))); % minmod
end


end

