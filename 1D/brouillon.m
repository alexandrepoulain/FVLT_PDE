   
    % fluxes at interface: first-order upwind
%     fmip12(3:nbx-1) = (m(3:nbx-1)).*max((p(3:nbx-1)-p(2:nbx-2))/h,0)+...
%          (m(2:nbx-2)).*min((p(3:nbx-1)-p(2:nbx-2))/h,0);
%     fnip12(3:nbx-1) = (n(3:nbx-1)).*max((p(3:nbx-1)-p(2:nbx-2))/h,0)+...
%          (n(2:nbx-2)).*min((p(3:nbx-1)-p(2:nbx-2))/h,0);
%     
%     nu_mip12 = zeros(nbx+1,1); nu_nip12= zeros(nbx+1,1);
%     nu_mip12(3:nbx-1) = (p(3:nbx-1)-p(2:nbx-2))/h*dt/h;
%     nu_nip12(3:nbx-1) = (p(3:nbx-1)-p(2:nbx-2))/h*dt/h;
%      
%     % Coumpute index of upwind nodes
%     index_j = zeros(nbx+1,1);
%     index_j(3:nbx-1) = (3:nbx-1)'-sign((p(3:nbx-1)-p(2:nbx-2))/h);
%     
%     fluxm_h(3:nbx-1) = fmip12(3:nbx-1) + (1-abs(nu_mip12(3:nbx-1))).*abs((p(3:nbx-1)-p(2:nbx-2))/h)*0.5...
%             .*(TVD_limiter(m(3:nbx-1)-m(2:nbx-2),m(index_j(3:nbx-1))-m(index_j(3:nbx-1)-1)));
%         
%     fluxn_h(3:nbx-1) = fnip12(3:nbx-1) + (1-abs(nu_nip12(3:nbx-1))).*abs((p(3:nbx-1)-p(2:nbx-2))/h)*0.5...
%         .*(TVD_limiter(n(3:nbx-1)-n(2:nbx-2),n(index_j(3:nbx-1))-n(index_j(3:nbx-1)-1)));
%     
%     
%    
%     % boundary conditions: zeros Neumann. Nothing to add because the flux
%     % is zero on the boundaries
%     ii =2;
%     if ((p(ii)-p(ii-1)) <= 0)
%         jj = ii-1;
%         fmip12 = (m(ii-1))*((p(ii)-p(ii-1))/h) ;
%         fnip12 = (n(ii-1))*((p(ii)-p(ii-1))/h);
%     else
%         jj= ii+1;
%         fmip12 = (m(ii))*((p(ii)-p(ii-1))/h) ;
%         fnip12 = (n(ii))*((p(ii)-p(ii-1))/h);
%     end
%     nu_mip12(ii) = (p(ii)-p(ii-1))/h*dt/h;
%     nu_nip12(ii) = (p(ii)-p(ii-1))/h*dt/h;
%     if jj~= 1
%         fmip12 = fmip12 + (1-abs(nu_mip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%                 *(TVD_limiter(m(ii)-m(ii-1),m(jj)-m(jj-1)));
%         fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%             *(TVD_limiter(n(ii)-n(ii-1),n(jj)-n(jj-1)));
%         fluxm_h(ii) = fmip12;
%         fluxn_h(ii) = fnip12;
%     else
%         fmip12 = fmip12 + (1-abs(nu_mip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%                 *(TVD_limiter(m(ii)-m(ii-1),m(ii)-m(ii-1)));
%         fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%             *(TVD_limiter(n(ii)-n(ii-1),n(ii)-n(ii-1)));
%         fluxm_h(ii) = fmip12;
%         fluxn_h(ii) = fnip12;
%     end
%     
% 
%     
%     ii =nbx;
%     if ((p(ii)-p(ii-1)) <= 0)
%         jj = ii-1;
%         fmip12 = (m(ii-1))*((p(ii)-p(ii-1))/h) ;
%         fnip12 = (n(ii-1))*((p(ii)-p(ii-1))/h);
%     else
%         jj= ii+1;
%         fmip12 = (m(ii))*((p(ii)-p(ii-1))/h) ;
%         fnip12 = (n(ii))*((p(ii)-p(ii-1))/h);
%     end
%     nu_mip12(ii) = (p(ii)-p(ii-1))/h*dt/h;
%     nu_nip12(ii) = (p(ii)-p(ii-1))/h*dt/h;
%     if jj~=nbx+1
%         fmip12 = fmip12 + (1-abs(nu_mip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%                 *(TVD_limiter(m(ii)-m(ii-1),m(jj)-m(jj-1)));
%         fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%             *(TVD_limiter(n(ii)-n(ii-1),n(jj)-n(jj-1)));
%         fluxm_h(ii) = fmip12;
%         fluxn_h(ii) = fnip12;
%     else
%         fmip12 = fmip12 + (1-abs(nu_mip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%                 *(TVD_limiter(m(ii)-m(ii-1),m(ii)-m(ii-1)));
%         fnip12 = fnip12 + (1-abs(nu_nip12(ii)))*abs((p(ii)-p(ii-1))/h)*0.5...
%             *(TVD_limiter(n(ii)-n(ii-1),n(ii)-n(ii-1)));
%         fluxm_h(ii) = fmip12;
%         fluxn_h(ii) = fnip12;
%     end
    