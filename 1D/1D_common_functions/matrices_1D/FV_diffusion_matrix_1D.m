function [K] = FV_diffusion_matrix_1D(dx,nbx,BC)
%FV_diffusion_matrix_1D This function compute the matrix associated to the
%laplacian 

K = zeros(nbx,nbx);

if BC == 'Homo-Neumann'
    for ii=2:nbx-1
        K(ii,ii) = -2/dx.^2;
        K(ii,ii+1) = 1/dx.^2;
        K(ii,ii-1) = 1/dx.^2;
    end
    ii=1;
    K(ii,ii) = -1/dx.^2;
    K(ii,ii+1) = 1/dx.^2;
    ii=nbx;
    K(ii,ii-1) = 1/dx.^2;
     K(ii,ii) = -1/dx.^2;

end



end

