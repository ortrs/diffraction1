function [v_out,V_out,V_in,Kx,Ky,Prop,prop]=propagation(v_in,x,y,z,lambda,sw)    
% Propagation over a predefined distance    
% Solution to the Helmholtz equation    
% All spatial coordinates in the unit of the wavelength    
% function call [v_out,V_out,V_in,Kx,Ky]=propagation(v_in,x,y,z);    
% v_in: input field (2D-matrix)    
% x,y: Spatial coordinates (2D-matrices)    
% z: propagation distance    
% v_out: output field in the spatial domain (2D-matrix)    
% V_out: output field in the Fourier domain (2D-matrix)    
% V_in:  input field in Fourier domain (2D-matrix)    
% Kx,Ky: Coordinates of the Fourier space (2D-matrix)
% sw: Use one or another method available for computation:
    % "manDFT" = Manual DFT. Slow and probably would crash
    % "redDFT" = reduced DFT. Using vectorization
    % "biDFT" = built-in DFT. 
    
% 
%tic;

    k = 2*pi/lambda;
    %Kx = 2*pi./x;
    %Ky = 2*pi./y;
    [Nx,Ny] = size(v_in);
    const = (1./(Nx*Ny));
    Kx = linspace(Nx/4/min(min(x)),Ny/4/max(max(x)),Nx);
    Ky = linspace(Ny/4/min(min(y)),Ny/4/max(max(y)),Ny);
    %Kx = (Nx.*lambda)./(2)./x;
    %Ky = (Ny.*lambda)./(2)./y;
    [Kx Ky] = meshgrid(Kx,Ky);
    
obj_size = size(v_in);
phy_x = max(max(x)); 
phy_y = max(max(y)); 
    
Fs_x = obj_size(2)/phy_x; 
Fs_y = obj_size(1)/phy_y; 
dx2 = Fs_x^(-1); 
dy2 = Fs_y^(-1); 
%x2 = dx2*(0:(obj_size(2) - 1))'; 
%y2 = dy2*(0:(obj_size(1) - 1))'; 
dFx = Fs_x/obj_size(2); 
dFy = Fs_y/obj_size(1); 
Fx = (-Fs_x/2:dFx:(Fs_x/2)); 
Fy = (-Fs_y/2:dFy:(Fs_y/2)); 
%Fx = (-Fs_x/2:dFx:(Fs_x/2 - dFx)); %this applies for non-symmetric things

%Fy = (-Fs_y/2:dFy:(Fs_y/2 - dFy)); 

alpha = lambda.*Fx; 
beta = lambda.*Fy; % gamma 

gama = zeros(length(beta), length(alpha)); 
for j = 1:length(beta);     
   for i = 1:length(alpha);         
       if (alpha(i)^2 + beta(j)^2) > 1;             
           gama(j, i) = 0;         
       else
           gama(j, i) = sqrt(1 - alpha(i)^2 - beta(j)^2);        
       end;     
   end; 
end; 

    
    if (strcmp(sw,'biDFT'))
       V_in =  fftshift(fft2(v_in));
       %Prop = exp ((1i).*sqrt(k^2 - Kx.^2 -Ky.^2).* z);
       Prop = sqrt(k^2 - Kx.^2 -Ky.^2).* z;
       prop = ifft2(fftshift(Prop));
       V_out = V_in .*  exp (1i.*sqrt(k^2 - Kx.^2 -Ky.^2).*z);
       %V_out = fftshift(V_out);
       v_out = ifft2(ifftshift(V_out));
    end
       
    if (strcmp(sw,'redDFT'))
        v_in = reshape (v_in,1,Nx*Ny);
        V_in = zeros (1,Nx.*Ny);
        m=0:Nx*Ny-1;
        m=repmat(m,Nx*Ny,1);
        m=m.*m';
        m = exp(-1i.*(2*pi/(Nx*Ny).*(m)));
        mi = exp(+1i.*(2*pi/(Nx*Ny).*(m)));
        V_out = (m*v_in')';
        V_out = reshape(V_out,Nx,Ny);
        %fft-shift re-shaping
        V_out = [circshift(V_out(1:Ny/2,1:end),[0 Nx/2-1]);
                 circshift(V_out(Ny/2+1:end,1:end),[0 Nx/2])];
        V_out = circshift(V_out,[Nx/2 0]);
        %V_out =  [V_out(Nx/2+1:end,Ny/2+1:end), ...
        %          V_out(Nx/2+1:end,1:Ny/2);
        %          V_out(1:Nx/2,Ny/2+1:end), ...
        %          V_out(1:Nx/2,1:Ny/2)];  
        %wave propagation
        V_out = V_out .* ...
             exp (1i .*sqrt(k^2 - Kx.^2 -Ky.^2).* z);
        trans = exp (1i .*sqrt(k^2 - Kx.^2 -Ky.^2).* z);
        %fft-shift
        %V_cout = circshift(V_out,[Ny/2 Nx/2]);
        V_cout = [circshift(V_out(1:Ny/2,1:end),[0 Nx/2-1]);
                 circshift(V_out(Ny/2+1:end,1:end),[0 Nx/2])];
        V_cout = circshift(V_out,[Nx/2 0]);
        %V_cout =  [V_out(Nx/2+1:end,Ny/2+1:end), ...
        %          V_out(Nx/2+1:end,1:Ny/2);
        %          V_out(1:Nx/2,Ny/2+1:end), ...
        %          V_out(1:Nx/2,1:Ny/2)];
        v_out=reshape (V_cout,1,Nx*Ny);
        v_out = (1/(Nx*Ny)).*(m*v_out')';
        v_out = reshape(v_out,Nx,Ny);
        %mesh(x,y,real(v_out));
        %figure;
        %mesh(x,y,imag(v_out));
    end
 
    if (strcmp(sw,'manDFT'))

        for m=1:Nx*Ny
            for n=1:Nx*Ny
                V_in(m) =  v_in(n) .* exp(-1i .* (2*pi*const).*m.*n); 
            end
        end
    end
    %V_in = 1/const .* reshape(Nx,Ny);
    %tac = 1;
end

