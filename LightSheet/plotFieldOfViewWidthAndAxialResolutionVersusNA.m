% [numericalApertures,FOVWidths,axialResolutions]=plotFieldOfViewWidthAndAxialResolutionVersusNA(refractiveIndex,wavelength)
%
% Calculates the field of view and axial resolution of a light sheet microscope as a function of the illumination NA.
%
% All SI units.
%
function [numericalApertures,FOVWidths,axialResolutions]=plotFieldOfViewWidthAndAxialResolutionVersusNA(refractiveIndex,wavelength)
    close all;
    
    if (nargin<1 || isempty(refractiveIndex))
        refractiveIndex=1.33;
    end
    if (nargin<2 || isempty(wavelength))
        wavelength=.5;
    end
    
    
    numericalApertures=.001:.001:refractiveIndex;
    
    % Calculate
    
%     FOVWidths_old=2*wavelength*sqrt(1-(numericalApertures./refractiveIndex).^2)*refractiveIndex./(2*sqrt(2)*numericalApertures.^2);
    FOVWidths=4*wavelength*refractiveIndex./numericalApertures.^2;
    axialResolutions=wavelength/2./numericalApertures;
        
    if (nargout==0)
        % Display results
        semilogy(numericalApertures,FOVWidths,'LineWidth',3,'Color',[.33 .33 1]);
        hold on;
        semilogy(numericalApertures,axialResolutions,'LineWidth',2,'Color',[0 .33 0]);
        % Add horizontal gray lines
        semilogy(repmat([0 refractiveIndex],[4 1]).',repmat([.5 1 5 10].',[1 2]).','LineWidth',1,'Color',[0.8 0.8 0.8]);
        xlim([0 refractiveIndex]);
        ylim([.1 500]);
        xlabel('NA');
        ylabel('FOV width and axial resolution [\mum]');
        legend({'Gaussian FOV Width','Gaussian axial resolution'});
        hold off;
    
        clear numericalApertures;
    end
end