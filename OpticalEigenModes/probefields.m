function [x,y,z,t,nx,ny,nz,nt,ff]=probefields(type,num,theta)

%%%
%type=1 bessel
%   param(1)=omega   %optical frequency :wavelength=2 pi c/(omega)
%   param(2)=gamma  %cone anlge
%   param(3)=pol    %A polaraisation 1:vert; 2:horiz; 3:longi;
%   param(4)=L      %L number
%   param(5)=n      %index of refraction

pr=probemodes(type,num,theta);
x1=-2e-6:20e-9:2e-6;
x1=x1;
y1=x1;
[x,y]=meshgrid(x1,y1);
z=0*x;t=0*x;
 
for ipar=1:size(pr,1)
[efx,efy,efz,hfx,hfy,hfz]=field(type,x,y,z,t,pr(ipar,:));
% pf(ipar).efx=efx;
% pf(ipar).efy=efy;
% pf(ipar).efz=efz;
% pf(ipar).hfx=hfx;
% pf(ipar).hfy=hfy;
% pf(ipar).hfz=hfz;
ff(1,ipar,:)=reshape(efx,1,[]);
ff(2,ipar,:)=reshape(efy,1,[]);
ff(3,ipar,:)=reshape(efz,1,[]);
ff(4,ipar,:)=reshape(hfx,1,[]);
ff(5,ipar,:)=reshape(hfy,1,[]);
ff(6,ipar,:)=reshape(hfz,1,[]);
end

x=reshape(x,1,[]);
y=reshape(y,1,[]);
z=reshape(z,1,[]);
t=reshape(t,1,[]);
nx=size(x1,2);
ny=size(y1,2);
nz=1;
nt=1;

