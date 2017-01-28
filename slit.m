function	[v,x,y]=slit(xa,ya,xb,yb,Nx,Ny)
%	Function	generates	around	the	centre	of	the	coordinate	
%	a	rectangular	amplitude	distribu4on	with	the	value	1	
%	Function call [v,x,y]=rectangle(xa,ya,xb,yb,Nx,Ny);	
%	v: final field	
%	x: matrix of the x-coordinates	
%	y: matrix of the y-coordinates	
%	xa:	size of	the	entire	domain	in	x-direction	
%	ya:	size of	the	entire	domain	in	y-direction	
%	xb:	size of	the	rectangle	in	x-direction	
%	yb:	size of	the	rectangle	in	y-direction	
%	Nx:	Number of points	in	x	direction	
%	Ny:	Number of points	in	y	direction
% Function call: [v,x,y]=slit(40,40,4,8,512,512);

% Initialize field

%%%%% Explicit %%%%%%%

%Make sure the index is odd (symmetric wrt point count)
ax = 2*round((xb/xa*Nx)/2);
ay = 2*round((yb/ya*Ny)/2);
%Create aperture
v = ones (ax,ay);
%Pad with zeroes around aperture
%Using floor to ensure it always rounds to 512
v = padarray(v,[floor((Nx-ax)/2) floor((Ny-ay)/2)]);
%v = padarray(v,[230 205]); %debug line
x = linspace(-xa/2,xa/2,Nx);
y = linspace(-xa/2,xa/2,Nx);
[x,y] = meshgrid(x,y);

%%% Non-explicit %%%
%ax = 2*round((xb/xa*Nx)/2);
%ay = 2*round((yb/ya*Ny)/2);
%v = padarray(ones(ax,ay),[floor((Nx-ax)/2) floor((Ny-ay)/2)]);
%[x,y] = meshgrid(linspace(-xa/2,xa/2,Nx),linspace(-xa/2,xa/2,Nx));
end

