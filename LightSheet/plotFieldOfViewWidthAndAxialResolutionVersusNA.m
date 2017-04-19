% [numericalApertures,FOVWidths,axialResolutions]=plotFieldOfViewWidthAndAxialResolutionVersusNA(refractiveIndex,wavelength)
%
% Calculates the field of view and axial resolution of a light sheet microscope as a function of the illumination NA.
%
% All SI units.
%
function [numericalApertures,FOVWidths,axialResolutions]=plotFieldOfViewWidthAndAxialResolutionVersusNA(refractiveIndex,wavelength)
    close all;
    
    if (nargin<1 || isempty(refractiveIndex))
        refractiveIndex=1.4;
    end
    if (nargin<2 || isempty(wavelength))
        wavelength=532e-9;
    end    
    
    numericalApertures=.001:.001:(refractiveIndex/1);
    
    alpha=7;
    betas=[0.10 0.05];
    
    fontSize=24;
    lineColors=[.33 .33 .33; 0.2 0.2 0.8; 0.8 0 0; 0 .67 0];
    lineStyles={':','-','--','-'};
    lineWidths=[2 3 2 3];
    
    % Calculate
    FOVGaussian=4*wavelength*refractiveIndex*numericalApertures.^-2;
    axialResolutionGaussian=wavelength./(2*numericalApertures*0.88);
    FOVBessel10=(wavelength/refractiveIndex)./(2*(1-sqrt(1-(numericalApertures/refractiveIndex).^2))*betas(1));
    axialResolutionBessel10=wavelength*pi*.05./(numericalApertures*betas(1));
    FOVBessel5=(wavelength/refractiveIndex)./(2*(1-sqrt(1-(numericalApertures/refractiveIndex).^2))*betas(2));
    axialResolutionBessel5=wavelength*pi*.05./(numericalApertures*betas(2));
    FOVAiry=(6*alpha*wavelength/refractiveIndex)./(1-sqrt(1-(numericalApertures/refractiveIndex).^2));
    maxSpFreqAiry=min(0.88,1/(0.05^2*48*alpha)).*numericalApertures/(wavelength/2);
    axialResolutionAiry=1./maxSpFreqAiry;
    
    maxDetectionNAs=refractiveIndex*sin(pi/2-asin(numericalApertures/refractiveIndex));
    maxDetectionNAs=maxDetectionNAs*0+0.40;
    spotLengths=4*refractiveIndex*wavelength./(maxDetectionNAs.^2);
    pixelSize=[1 1]*7.4e-6;
    realMagnification=20*0.176/0.160;
    depthOfFocus=refractiveIndex*(wavelength./(maxDetectionNAs.^2)+pixelSize(1)./(realMagnification.*maxDetectionNAs));
    
%     FOVWidths_old=2*wavelength*sqrt(1-(numericalApertures./refractiveIndex).^2)*refractiveIndex./(2*sqrt(2)*numericalApertures.^2);
    FOVWidths=4*wavelength*refractiveIndex./numericalApertures.^2;
    axialResolutions=wavelength/2./numericalApertures;
        
    if (nargout==0)
        % Display results
        figs(1)=figure('Position',[100 100 1024 768]);
        axs(1)=subplot(2,2,1,'Parent',figs(1));
        axs(2)=subplot(2,2,2,'Parent',figs(1));
        axs(3)=subplot(2,2,[3 4],'Parent',figs(1));
        labels=[];
        semilogy(numericalApertures,axialResolutionGaussian*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionBessel10*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',lineColors(2,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionBessel5*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionAiry*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(1)); hold(axs(1),'on');
        semilogy(numericalApertures,axialResolutionGaussian*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(1)); hold(axs(1),'on');
%         semilogy(numericalApertures,spotLengths*1e6,':','LineWidth',2,'Color',[1 1 1]*.33,'Parent',axs(1)); hold(axs(1),'on');
        xlim(axs(1),[0 numericalApertures(end)]); ylim(axs(1),[.1 500]);
        set(axs(1),'YTick',10.^[0:1:100]);
        labels(end+1)=xlabel(axs(1),'NA');
        labels(end+1)=ylabel(axs(1),'Axial Resolution [\mum]');
        legend(axs(1),{'Gaussian','Bessel10','Bessel5','Airy'});
        semilogy(numericalApertures,FOVGaussian*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(2)); hold(axs(2),'on');
        semilogy(numericalApertures,FOVBessel10*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',lineColors(2,:),'Parent',axs(2)); hold(axs(2),'on');
        semilogy(numericalApertures,FOVBessel5*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(2)); hold(axs(2),'on');
        semilogy(numericalApertures,FOVAiry*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(2)); hold(axs(2),'on');
        xlim(axs(2),[0 numericalApertures(end)]); ylim(axs(2),[.1 500]);
        set(axs(2),'YTick',10.^[0:1:100]);
        labels(end+1)=xlabel(axs(2),'NA');
        labels(end+1)=ylabel(axs(2),'FOV Width [\mum]');
        legend(axs(2),{'Gaussian','Bessel10','Bessel5','Airy'});
%         loglog(axialResolutionGaussian*1e6,FOVGaussian*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(3)); hold(axs(3),'on');
%         loglog(axialResolutionBessel10*1e6,FOVBessel10*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',[0.8 0.2 0.2],'Parent',axs(3)); hold(axs(3),'on');
%         loglog(axialResolutionBessel5*1e6,FOVBessel5*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(3)); hold(axs(3),'on');
%         loglog(axialResolutionAiry*1e6,FOVAiry*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVGaussian*1e6,axialResolutionGaussian*1e6,lineStyles{1},'LineWidth',lineWidths(1),'Color',lineColors(1,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVBessel10*1e6,axialResolutionBessel10*1e6,lineStyles{2},'LineWidth',lineWidths(2),'Color',lineColors(2,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVBessel5*1e6,axialResolutionBessel5*1e6,lineStyles{3},'LineWidth',lineWidths(3),'Color',lineColors(3,:),'Parent',axs(3)); hold(axs(3),'on');
        plot(FOVAiry*1e6,axialResolutionAiry*1e6,lineStyles{4},'LineWidth',lineWidths(4),'Color',lineColors(4,:),'Parent',axs(3)); hold(axs(3),'on');
        xlim(axs(3),[0 1000]); ylim(axs(3),[0 10]);
        set(axs(3),'YTick',[0:1:100]);
        labels(end+1)=xlabel(axs(3),'FOV Width [\mum]');
        labels(end+1)=ylabel(axs(3),'Axial Resolution [\mum]');
        legend(axs(3),{'Gaussian','Bessel10','Bessel5','Airy'},'Parent',figs(1));
        set(axs,'LineWidth',3,'FontSize',fontSize,'FontWeight','bold');
        set(labels,'FontSize',fontSize,'FontWeight','bold');
        
%         figs(1)=figure('Position',[100 100 1024 768]);
%         axs(1)=axes('Parent',figs(1));
%         % Add horizontal gray lines
%         horizontalLinePositions=[.5 1 5 10 50 100];
%         semilogy(repmat([0 refractiveIndex].',size(horizontalLinePositions)),repmat(horizontalLinePositions,[2 1]),'LineWidth',lineWidths(1)1,'Color',[0.8 0.8 0.8],'Parent',axs(1));
%         hold(axs(1),'on');
%         semilogy(numericalApertures,FOVWidths*1e6,'LineWidth',lineWidths(1)3,'Color',lineColors(1,:),'Parent',axs(1));
%         semilogy(numericalApertures,axialResolutions*1e6,'LineWidth',lineWidths(1)2,'Color',[0 .33 0],'Parent',axs(1));
%         xlim([0 refractiveIndex]);
%         ylim([.1 500]);
%         xlabel('NA','FontSize',fontSize,'FontWeight','bold');
%         ylabel('FOV width and axial resolution [\mum]','FontSize',fontSize,'FontWeight','bold');
%         legend({'Gaussian FOV Width','Gaussian axial resolution'});
%         hold off;
%         set(axs(1),'LineWidth',lineWidths(1)3,'FontSize',fontSize,'FontWeight','bold');
%         set(axs(1),'XTick',[0:.2:10]);
%                 
%         figs(2)=figure('Position',[100 100 1024 768]);
%         axs(2)=axes('Parent',figs(2));
%         hold(axs(1),'on');
%         plot(FOVWidths*1e6,axialResolutions*1e6,'LineWidth',lineWidths(1)3,'Color',lineColors(1,:),'Parent',axs(2));
%         xlim([0 5]);
% %         ylim([0 500]);
%         xlabel('FOV width and axial resolution [\mum]','FontSize',fontSize,'FontWeight','bold');
%         ylabel('approx. light sheet width [\mum]','FontSize',fontSize,'FontWeight','bold');
%         legend({'Gaussian'});
%         hold off;
%         set(axs(2),'LineWidth',lineWidths(1)3,'FontSize',fontSize,'FontWeight','bold');
%         set(axs(2),'XTick',[0:1:10]);
%     
%         print(figs(1),'plotFieldOfViewWidthAndAxialResolutionVersusNA.eps','-depsc')
%         print(figs(2),'plotFieldOfViewWidthVersusAxialResolution.eps','-depsc')
        
        clear numericalApertures;
    end
end