function Coils = Simulate_Coils(Sz,NCoil)
  
% SENSITIVITY MAP creates aritficial sensitivity maps
%
%   The map is created such as the Sum of Squares is 1
%
%                      sum(Coils.^2,3)=ones(Sz)
%
%  INPUTS:
%               Sz=[Mx,My,Mz]   Size of the map
%               NCoil        Number of coils   (even number)
%
%  OUTPUT
%               Coils     [Mx,My,Mz,NCoils] Sensitivity Map
%
% EXAMPLE
%              Coils=Simulate_Coils([256,256],4);
%
% Santiago Aja-Fernandez
% Parallel Imaging Toolbox
% Feb. 2012
% www.lpi.tel.uva.es/~santi
  
  
Mx=Sz(1);
My=Sz(2);
Mz=Sz(3);
ejex=1:My;
ejey=1:Mx;
ejez=1:Mz;
[vY,vX,vZ] = meshgrid(ejex,ejey,ejez);

vX=pi/2.*vX./max(vX(:));
vY=pi/2.*vY./max(vY(:));
vZ=pi/2.*vZ./max(vZ(:));
  
Coils=zeros([Sz,NCoil]);
Theta=0:(2*pi/NCoil):(2*pi-(2*pi/NCoil));
Theta=Theta(1:end-1);
Phi=0:(2*pi/NCoil):(2*pi-(2*pi/NCoil));
Phi=Phi(1:end-1);

for ii=1:NCoil/2
  if (Theta(ii)<=pi/2)
      N1=vX.*cos(Theta(ii)).*sin(Phi(ii))+...
          vY.*sin(Theta(ii)).*sin(Phi(ii))+...
          vZ.*cos(Phi(ii));
      N1=N1./max(N1(:)).*pi/2; 
      Coils(:,:,:,ii)=cos(N1);
      Coils(:,:,:,NCoil./2+ii)=sin(N1);
  else %>pi/2
      N1=vX.*abs(cos(Theta(ii))).*sin(Phi(ii))+...
          (pi/2-vY).*sin(Theta(ii)).*sin(Phi(ii))+...
          vZ.*abs(cos(Phi(ii)));

      N1=N1./max(N1(:)).*pi/2; 
      Coils(:,:,:,ii)=sin(N1);
      Coils(:,:,:,NCoil./2+ii)=cos(N1);
        
  end
  
end
  
Coils=Coils./sqrt(NCoil/2);