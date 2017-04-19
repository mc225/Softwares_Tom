% RGB=spectrumToRGB(wavelengths,wavelengthIntensities,constrainValues)
%
% Returns the red green and blue component to simulate a given spectrum.
% The number of dimensions of the returned matrix is equal to that of the
% wavelengthIntensities argument, and the dimensions are the same except
% for the last dimension which is always three.
%
% wavelengths can be a vector of wavelengths or a function handle returning the intensity for a given wavelength, in the latter case the argument wavelengthIntensities is ignored.
% If wavelengthIntensities is a matrix, its last dimension should equal the number of wavelengths.
%
% constrainValues (optional, default=true): a boolean to indicate if negative RGB values should be de-saturated to be non-negative
% 
function RGB=spectrumToRGB(wavelengths,wavelengthIntensities,constrainValues)
    if (nargin<1 || isempty(wavelengths))
        wavelengths=532e-9; %Nd-YAG, frequency-doubled
%         wavelengths=632.8e-9; %He-Ne
        wavelengths=1e-9*[300:800];
    end
    % From John Walker's http://www.fourmilab.ch/documents/specrend/specrend.c
    wavelengthsXYZTable=[380:5:780]*1e-9;
    if (isa(wavelengths,'function_handle'))
        wavelengthIntensities=wavelengths(wavelengthsXYZTable);
        wavelengths=wavelengthsXYZTable;
    else
        if (nargin>=1)
            if (nargin<2 || isempty(wavelengthIntensities))
                wavelengthIntensities=ones(size(wavelengths));
            end
        else
            wavelengthIntensities=permute(eye(length(wavelengths)),[3 1 2]);
        end
    end
    if (nargin<3 || isempty(constrainValues))
        constrainValues=true;
    end
    
    wavelengths=wavelengths(:);
    
    inputSize=size(wavelengthIntensities);
    nbWavelengths=numel(wavelengths);
    if (nbWavelengths==1 && inputSize(end)~=1)
        inputSize(end+1)=1;
    end
    wavelengthIntensities=reshape(wavelengthIntensities,[],nbWavelengths).';
    
    if (nbWavelengths~=inputSize(end))
        logMessage('spectrumToRGB Error: the number of wavelengths should be equal to the last dimension of the second argument.');
        return;
    end

    % Define Color Systems
    %  

    %White point chromaticities
    IlluminantC = [0.3101, 0.3162];	    	% For NTSC television 
    IlluminantD65 = [0.3127, 0.3291];    	% For EBU and SMPTE 
    IlluminantE = [0.33333333, 0.33333333]; % CIE equal-energy illuminant 
    % Gamma of nonlinear correction.   See Charles Poynton's ColorFAQ Item 45 and GammaFAQ Item 6 at:   http://www.poynton.com/ColorFAQ.html and http://www.poynton.com/GammaFAQ.html
    GAMMA_REC709=0;		% Rec. 709

    NTSC=struct('RGB',[0.67,0.33; 0.21,0.71; 0.14,0.08].','WhitePoint',IlluminantC,'Gamma',GAMMA_REC709);
    EBU=struct('RGB',[0.64,0.33;0.29,0.60;0.15,0.06].','WhitePoint',IlluminantD65,'Gamma',GAMMA_REC709); % PAL/SECAM
    SMPTE=struct('RGB',[0.630,0.340;0.310,0.595;0.155,0.070].','WhitePoint',IlluminantD65,'Gamma',GAMMA_REC709);
    HDTV=struct('RGB',[0.670,0.330;0.210,0.710;0.150,0.060].','WhitePoint',IlluminantD65,'Gamma',GAMMA_REC709);
    CIE=struct('RGB',[0.7355,0.2645;0.2658,0.7243;0.1669,0.0085].','WhitePoint',IlluminantE,'Gamma',GAMMA_REC709);
    CIERec709=struct('RGB',[0.64,0.33;0.30,0.60;0.15,0.06].','WhitePoint',IlluminantD65,'Gamma',GAMMA_REC709);
    
    %Pick this color system:
    colorSystem=SMPTE;

    % From John Walker's http://www.fourmilab.ch/documents/specrend/specrend.c
	XYZTable=[0.0014,0.0000,0.0065;  0.0022,0.0001,0.0105;  0.0042,0.0001,0.0201; ... 
                0.0076,0.0002,0.0362;  0.0143,0.0004,0.0679;  0.0232,0.0006,0.1102; ... 
                0.0435,0.0012,0.2074;  0.0776,0.0022,0.3713;  0.1344,0.0040,0.6456; ... 
                0.2148,0.0073,1.0391;  0.2839,0.0116,1.3856;  0.3285,0.0168,1.6230; ... 
                0.3483,0.0230,1.7471;  0.3481,0.0298,1.7826;  0.3362,0.0380,1.7721; ... 
                0.3187,0.0480,1.7441;  0.2908,0.0600,1.6692;  0.2511,0.0739,1.5281; ... 
                0.1954,0.0910,1.2876;  0.1421,0.1126,1.0419;  0.0956,0.1390,0.8130; ... 
                0.0580,0.1693,0.6162;  0.0320,0.2080,0.4652;  0.0147,0.2586,0.3533; ... 
                0.0049,0.3230,0.2720;  0.0024,0.4073,0.2123;  0.0093,0.5030,0.1582; ... 
                0.0291,0.6082,0.1117;  0.0633,0.7100,0.0782;  0.1096,0.7932,0.0573; ... 
                0.1655,0.8620,0.0422;  0.2257,0.9149,0.0298;  0.2904,0.9540,0.0203; ... 
                0.3597,0.9803,0.0134;  0.4334,0.9950,0.0087;  0.5121,1.0000,0.0057; ... 
                0.5945,0.9950,0.0039;  0.6784,0.9786,0.0027;  0.7621,0.9520,0.0021; ... 
                0.8425,0.9154,0.0018;  0.9163,0.8700,0.0017;  0.9786,0.8163,0.0014; ... 
                1.0263,0.7570,0.0011;  1.0567,0.6949,0.0010;  1.0622,0.6310,0.0008; ... 
                1.0456,0.5668,0.0006;  1.0026,0.5030,0.0003;  0.9384,0.4412,0.0002; ... 
                0.8544,0.3810,0.0002;  0.7514,0.3210,0.0001;  0.6424,0.2650,0.0000; ... 
                0.5419,0.2170,0.0000;  0.4479,0.1750,0.0000;  0.3608,0.1382,0.0000; ... 
                0.2835,0.1070,0.0000;  0.2187,0.0816,0.0000;  0.1649,0.0610,0.0000; ... 
                0.1212,0.0446,0.0000;  0.0874,0.0320,0.0000;  0.0636,0.0232,0.0000; ... 
                0.0468,0.0170,0.0000;  0.0329,0.0119,0.0000;  0.0227,0.0082,0.0000; ... 
                0.0158,0.0057,0.0000;  0.0114,0.0041,0.0000;  0.0081,0.0029,0.0000; ... 
                0.0058,0.0021,0.0000;  0.0041,0.0015,0.0000;  0.0029,0.0010,0.0000; ... 
                0.0020,0.0007,0.0000;  0.0014,0.0005,0.0000;  0.0010,0.0004,0.0000; ... 
                0.0007,0.0002,0.0000;  0.0005,0.0002,0.0000;  0.0003,0.0001,0.0000; ... 
                0.0002,0.0001,0.0000;  0.0002,0.0001,0.0000;  0.0001,0.0000,0.0000; ... 
                0.0001,0.0000,0.0000;  0.0001,0.0000,0.0000;  0.0000,0.0000,0.0000];
    
    XYZPerWavelength=interp1(wavelengthsXYZTable,XYZTable,wavelengths,'cubic',0).'; %XYZ in rows
    XYZ=XYZPerWavelength*wavelengthIntensities; % XYZ in rows, data points in columns
    XYZ=XYZ./repmat(max(eps(1),sum(XYZ,1)),[3 1]);
    
    % Transform XYZ coordinate system to RGB coordinate system
    XYZRGBTransformMatrix=cat(2,colorSystem.RGB.',1-sum(colorSystem.RGB.',2)); %XYZ in columns
    whitePoint=[colorSystem.WhitePoint 1-sum(colorSystem.WhitePoint)]; %XYZ in columns
    RGBXYZTransformMatrix=cross(XYZRGBTransformMatrix(:,[2 3 1]),XYZRGBTransformMatrix(:,[3 1 2])); %RGB in columns
    whitePointRGB=whitePoint*RGBXYZTransformMatrix/whitePoint(2);
    RGBXYZTransformMatrix=RGBXYZTransformMatrix./repmat(whitePointRGB,[3 1]);
    RGB=RGBXYZTransformMatrix*XYZ; % RGB in rows, data points in columns
    
    if (constrainValues)
        %Constrain RGB to non-negative values
        approximated=min(RGB,[],1)<0;
        RGB = RGB-repmat(min(0,min(RGB,[],1)),[3 1]);
    else
        approximated=false;
    end
        
    %Adjust to intensity of input (TODO, probably needs to be log.-scaled)
    RGB=RGB.*repmat(sum(wavelengthIntensities,1),[3 1]);
    
    RGB=reshape(RGB.',[inputSize(1:end-1) 3]);
    
    %Display
    if (nargout==0 && ndims(RGB)>=3)
        RGB=RGB./max(eps(1),max(RGB(:)));
    
        %Make sure that it can be displayed as a color.
        RGB=max(0,RGB);
        
        figure;
%         axis off;
        image(wavelengths*1e9,1,RGB);
        xlabel('\lambda [nm]','Interpreter','Tex');
        
        if (any(approximated))
            logMessage('Approximated %u%% of the colors!',round(100*mean(approximated)));
        end
        
        clear RGB;
    end
end
    