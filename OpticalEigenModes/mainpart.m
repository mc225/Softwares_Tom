%%
function mainpart()
    lgct=[];
    %% new 
    %for na=namin:0.1:namax
    na=0.65;
    theta=asin(na);
    prModes=probemodes(5,16,theta);
    %% new end

    lambda=0.8e-6;
    q1=1.22*lambda/(2*na);%airy disk radius
    qb=2.4*lambda/2/pi/sin(theta); %BB radius first zero
    [x,y,z,t,nx,ny,nz,nt,ff]=probefields(5,16,theta);
    dx=((max(x)-min(x))/nx);
    %%

    efx=squeeze(ff(1,:,:));
    efy=squeeze(ff(2,:,:));
    efz=squeeze(ff(3,:,:));
    hfx=squeeze(ff(4,:,:));
    hfy=squeeze(ff(5,:,:));
    hfz=squeeze(ff(6,:,:));

    %figure(1)
    %ii=23;
    %imagesc(reshape(abs(efx(ii,:)).^2+abs(efz(ii,:)).^2+abs(efy(ii,:)).^2+...
    %    abs(hfx(ii,:)).^2+abs(hfy(ii,:)).^2+abs(hfz(ii,:)).^2,nx,ny))

    %%
    ii=0;
    rmax=3.65e-6;
    %rmax=.4e-6:.01e-6:.47e-6;

    dr=[];deff=[];
    for roi=rmax
        rr=x.^2+y.^2;
        rmask=rr<(roi)^2;
        op=repmat(rmask,size(efx,1),1);

        maif=real((efx.*op)*hfy'-(efy.*op)*hfx');

        [vec lam]=eig(maif);
        dlam=real(diag(lam));
        idx=(dlam>1e-2*max(dlam));
        iev=vec(:,idx);
        ilam=diag(1./sqrt(dlam(idx)));
        av=iev*ilam;
        
        op=repmat(rr.*rmask,size(efx,1),1);
        maif2=av'*real((efx.*op)*hfy'-(efy.*op)*hfx')*av;

        [vec2 lam2]=eig((maif2+maif2')/2);
        dr=[dr sqrt(lam2(1,1))/q1];
        deff=[deff 1/max(dlam)/((av*vec2(:,1))'*(av*vec2(:,1)))];
        [roi sqrt(lam2(1,1)) ]

        %new
        figure(11)
        qm=1;ii=ii+1;
        lgc=av*vec2(:,qm);
        lgc=lgc*sign(lgc(1));
        lgct(ii,:)=lgc;
        nlgc=lgc/abs(lgc(end));
        plot(abs(nlgc).^2);drawnow;

        thp0x=nlgc'*efx;thp0y=nlgc'*efy;thp0z=nlgc'*efz;
        thp1x=nlgc'*hfx;thp1y=nlgc'*hfy;thp1z=nlgc'*hfz;
        thp0=real(conj(thp0x).*thp1y-conj(thp0y).*thp1x);
        thp=reshape(thp0,nx,ny);
        ti2=sum(rmask.*thp0);

        figure(3)
        imagesc(thp);%title(['oei:' num2str(ti0/ti2)])
        figure(4);
        subplot(4,3,1);imagesc(reshape(real(thp0x),nx,ny));
        subplot(4,3,2);imagesc(reshape(real(thp0y),nx,ny));
        subplot(4,3,3);imagesc(reshape(real(thp0z),nx,ny));
        subplot(4,3,4);imagesc(reshape(imag(thp0x),nx,ny));
        subplot(4,3,5);imagesc(reshape(imag(thp0y),nx,ny));
        subplot(4,3,6);imagesc(reshape(imag(thp0z),nx,ny));
        subplot(4,3,7);imagesc(reshape(real(thp1x),nx,ny));
        subplot(4,3,8);imagesc(reshape(real(thp1y),nx,ny));
        subplot(4,3,9);imagesc(reshape(real(thp1z),nx,ny));
        subplot(4,3,10);imagesc(reshape(imag(thp1x),nx,ny));
        subplot(4,3,11);imagesc(reshape(imag(thp1y),nx,ny));
        subplot(4,3,12);imagesc(reshape(imag(thp1z),nx,ny));
        %%
        return;
        
        [xslm,yslm]=meshgrid(1:800,1:600);
        rslm=sqrt((xslm-400).^2+(yslm-300).^2);
        irs=rslm<260;
        figure(4)
        ma1=0*xslm;
        ma1(irs)=interp1(0:260/(16-1):260,nlgc,rslm(irs));
        imagesc(ma1)
        ma(1,:,:)=ma1;
        tname(1)={'OEi'};
        efi(1)=sum(sum(abs(ma1).^2));
        %save(['/Users/mm17/Desktop/oeimask.mat'],'ma','title');

        %%
        %%new
        nroi=rmax/wave;
        kv=sin(prModes(:,2)); % [omega, gamma]
    end
    %%new end
    %additional masks
    ma1=0*xslm;
    ma1(irs)=1;
    ma(2,:,:)=ma1;
    tname(2)={'Dif. Limited'};
    efi(2)=sum(sum(abs(ma1).^2));

    irs=rslm<260*sqrt(.8);
    ma1=0*xslm;
    ma1(irs)=1;
    ma(3,:,:)=ma1;
    tname(3)={'Truncated R=80%'};
    efi(3)=sum(sum(abs(ma1).^2));
    efi/efi(2)

    irs=rslm<260*sqrt(.6);
    ma1=0*xslm;
    ma1(irs)=1;
    ma(4,:,:)=ma1;
    tname(4)={'Truncated R=60%'};
    efi(4)=sum(sum(abs(ma1).^2));
    efi/efi(2)

    irs=rslm<260*sqrt(.4);
    ma1=0*xslm;
    ma1(irs)=1;
    ma(5,:,:)=ma1;
    tname(5)={'Truncated R=40%'};
    efi(5)=sum(sum(abs(ma1).^2));
    efi/efi(2)

    stop

    %%
    for ii=1:16
        lgc0=0*lgc;lgc0(ii)=1;
        thp0x=lgc0'*efx;thp0y=lgc0'*efy;thp0z=lgc0'*efz;
        thp1x=lgc0'*hfx;thp1y=lgc0'*hfy;thp1z=lgc0'*hfz;
        thp0=real(conj(thp0x).*thp1y-conj(thp0y).*thp1x);
        thp=reshape(thp0,nx,ny);

        figure(6);subplot(4,4,ii)
        imagesc(thp);axis equal
    end
    %

    %%
    qm=1;
    lgc=av*vec2(:,qm);
    thp0x=lgc'*efx;thp0y=lgc'*efy;thp0z=lgc'*efz;
    thp1x=lgc'*hfx;thp1y=lgc'*hfy;thp1z=lgc'*hfz;
    thp0=abs(thp0x.^2)+abs(thp0y.^2)+abs(thp0z.^2)+abs(thp1x.^2)+abs(thp1y.^2)+abs(thp1z.^2);

    ImageHandle=figure(2);
    thp=reshape(thp0,nx,ny);
    ma=reshape(rmask,nx,ny);
    xx=reshape(x,nx,ny);
    yy=reshape(y,nx,ny);

    zpos=thp;zpos(~ma)=NaN;
    zrgb=sc(zpos,'jet');
    zpos2=thp;zpos2(ma)=NaN;
    zrgb2=sc(zpos2,'gray');


    sh(1)=imagesc(zrgb2);
    %%%%%colorbar('location','east')
    hold on
    sh(2)=imagesc(zrgb);
    %%%%%set(sh,'edgecolor','none')
    set(sh(2),'alphadata',ma) 
    axis equal;
    %colorbar('location','west')

    max((1-rmask).*thp0)/max((rmask).*thp0);
    set(ImageHandle,'units','centimeters','position',[5 5 15 15]) % set the screen size and position
    set(ImageHandle,'paperunits','centimeters','paperposition',[6 6 14 14]) % set size and position for printing
    set(gca,'units','normalized','position',[0 0 1 1]) % make sure axis fills entire figure

    %%
    rad=sqrt(max(rr(rmask)))/((max(x)-min(x))/nx);

    ImageHandle=figure(4);
    thp=reshape(thp0,nx,ny)/max(thp0);
    hold on;imshow(thp*256,jet(256));axis equal; 

    plot((nx)/2+1+rad*cos(0:pi/30:2*pi),(ny)/2+1+rad*sin((0:pi/30:2*pi)),'--y','LineWidth',3);drawnow
    sqrt(2*sum(sum((abs(rmask.*(x.^2+y.^2).*thp0))))/sum(sum((abs(rmask.*thp0)))))
    
    set(ImageHandle,'units','centimeters','position',[5 5 15 15]) % set the screen size and position
    set(ImageHandle,'paperunits','centimeters','paperposition',[6 6 14 14]) % set size and position for printing
    set(gca,'units','normalized','position',[0 0 1 1]) % make sure axis fills entire figure

    stop
    %%
    figure(1)
    %quiver(reshape(real(thp0x),nx,ny),reshape(real(thp0y),nx,ny));axis equal; drawnow
    quiver(reshape(real(rmask.*thp0x),nx,ny),reshape(real(rmask.*thp0y),nx,ny));axis equal; drawnow
    figure(3)
    %quiver(reshape(real(thp0x),nx,ny),reshape(real(thp0y),nx,ny));axis equal; drawnow
    subplot(2,3,1);imagesc(reshape(real(thp0x),nx,ny));axis equal;
    subplot(2,3,2);imagesc(reshape(real(thp0y),nx,ny));axis equal;
    subplot(2,3,3);imagesc(reshape(real(thp0z),nx,ny));axis equal;
    subplot(2,3,4);imagesc(reshape(real(thp1x),nx,ny));axis equal;
    subplot(2,3,5);imagesc(reshape(real(thp1y),nx,ny));axis equal;
    subplot(2,3,6);imagesc(reshape(real(thp1z),nx,ny));axis equal;
    %% first intenisty eigemode
    figure(5)
    lgc=av(:,1);
    thp0x=lgc'*efx;thp0y=lgc'*efy;thp0z=lgc'*efz;
    thp1x=lgc'*hfx;thp1y=lgc'*hfy;thp1z=lgc'*hfz;
    thp0=abs(thp0x.^2)+abs(thp0y.^2)+abs(thp0z.^2)+abs(thp1x.^2)+abs(thp1y.^2)+abs(thp1z.^2);
    thp=reshape(rmask.*thp0,nx,ny);
    imagesc(thp);axis equal; drawnow
    
    sqrt(sum(sum((abs(rmask.*(x.^2+y.^2).*thp0))))/sum(sum((abs(rmask.*thp0)))))
    %% max Bessel
    ImageHandle=figure(6);
    lgc=zeros(size(av,1),1);lgc(end)=1;
    thp0x=lgc'*efx;thp0y=lgc'*efy;thp0z=lgc'*efz;
    thp1x=lgc'*hfx;thp1y=lgc'*hfy;thp1z=lgc'*hfz;
    thp0=abs(thp0x.^2)+abs(thp0y.^2)+abs(thp0z.^2)+abs(thp1x.^2)+abs(thp1y.^2)+abs(thp1z.^2);
    thp=reshape(thp0,nx,ny)/max(thp0);
    imshow(thp*256,jet(256));axis equal; hold on;
    rb=qb/dx;
    plot((nx)/2+1+rb*cos(0:pi/30:2*pi),(ny)/2+1+rb*sin((0:pi/30:2*pi)),'--y','LineWidth',3);drawnow
    hold off
    rmask=rr<(qb)^2;
    sqrt(2*sum(sum((abs(rmask.*(x.^2+y.^2).*thp0))))/sum(sum((abs(rmask.*thp0)))))

    set(ImageHandle,'units','centimeters','position',[6 6 15 15]) % set the screen size and position
    set(ImageHandle,'paperunits','centimeters','paperposition',[5 5 16 16]) % set size and position for printing
    set(gca,'units','normalized','position',[0 0 1 1]) % make sure axis fills entire figure
    Iseg = getframe(ImageHandle);
    imwrite(Iseg.cdata,'/Users/mm17/Documents/OLDmac/Desktop/QME-OPEX/maxbessel.1.png')

    %% Airy disk Bessel
    close(6)
    dx=((max(x)-min(x))/nx);
    ImageHandle=figure(6);
    thp0=(2*besselj(1,sqrt(rr)/q1*3.83)./sqrt(rr)).^2;
    thp=reshape(thp0,nx,ny)/max(thp0);
    imshow(thp*256,jet(256));axis equal; 
    hold on
    rb=q1/dx;
    plot((nx)/2+1+rb*cos(0:pi/30:2*pi),(ny)/2+1+rb*sin((0:pi/30:2*pi)),'--y','LineWidth',3);drawnow
    hold off;
    set(ImageHandle,'units','centimeters','position',[6 6 15 15]) % set the screen size and position
    set(ImageHandle,'paperunits','centimeters','paperposition',[5 5 16 16]) % set size and position for printing
    set(gca,'units','normalized','position',[0 0 1 1]) % make sure axis fills entire figure
    Iseg = getframe(ImageHandle);
    imwrite(Iseg.cdata,'/Users/mm17/Documents/OLDmac/Desktop/QME-OPEX/airy.1.png')

    rmask=rr<(q1)^2;
    sqrt(2*sum(sum((abs(rmask.*(x.^2+y.^2).*thp0))))/sum(sum((abs(rmask.*thp0)))))
end