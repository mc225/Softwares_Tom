% Program for tracking the center of mass of an optically-trapped spherical particle within a flow which can
% move the particle away from the trap and possibly including other particles. The basic idea is to track the
% position of the trapped particle, the image of which is sharp, using both a sharp edge detection and the dark
% area (zone of low transmission inside the bead).
%
% External function circfit.m required.
%
% MRDSV, University of St Andrews, 30-05-2012 (v.2)
% -------------------------------------------------------------------------

clear all; close all; clc; warning off all
tic         % Starts timer     

% -------------------------------------------------------------------------
% Variables
% -------------------------------------------------------------------------

% ----- Variables related to the files to be analyzed
cd 'C:\Users\Matthieu\Documents\Matlab particle tracking\Videos\2012-05-08\serie2\s2_32'
                        % Directory where the image sequence can be found
addpath('C:\Users\Matthieu\Documents\Matlab particle tracking');
                        % Parent directory in which the external functions to be called can be found
NameRoot = 's2_32_';
FirstFrame = 1;
LastFrame = 2115;

% ----- Physical variables
dt = 1/209;             % Time interval between 2 consecutive frames (= 1/frame rate, known at +/- 1 fps)
Diameter = 3.01;        % Estimated bead diameter (µm)
Calibration = 20/465;   % Pixel size (µm)
Xexpect = 357;          % Approx expected position of the bead center (pixels)
Yexpect = 125;

% ----- Edge/dark level detection variables
Dark = 60;              % Threshold for dark detection (gray level value)
GrHi = 0.2;             % Higher threshold for the edge detection (Canny method)
GrLo = 0.1;             % Lower threshold for the edge detection (Canny method)
Filter = 3;             % Radius of the Gaussian filter for edge detection (Canny method)
Tolerance = 1;          % Maximal amplitude of the bead motion to be measured (units of bead radius)

% ----- Operating preferences
OneOverTen = 1;         % To display one frame every 10 (binary)
Display = 100;          % Alternatively, last displayed frame
DispEdge = 1;           % To display edge points (binary)
DispDark = 1;           % To display dark points (binary)
DispPoints = 0;         % To display bead points (binary)
Saving = 1;             % To save the measured data (overwrites the previous saving) (binary)

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------

Frames = (FirstFrame:LastFrame)';
EstRad = Diameter/(2*Calibration);      % Estimated radius of the bead (pixels)
Area = pi*EstRad^2;                     % Estimated surface area of the bead (pixels)
X0 = zeros(size(Frames)); Y0 = zeros(size(Frames));

for n = FirstFrame:LastFrame
    n
    
    % ---------------------------------------------------------------------
    % Step 0. Reconstructing frame name, then loading frame
    % ---------------------------------------------------------------------

    if n < 10 
        Image = [NameRoot '000' num2str(n) '.bmp'];
    else
        if n < 100
            Image = [NameRoot '00' num2str(n) '.bmp'];
        else
            if n < 1000
                Image = [NameRoot '0' num2str(n) '.bmp'];    
            else
                if n > 999
                    Image = [NameRoot '' num2str(n) '.bmp'];
                end
            end
        end
    end
    eval(['imread ' Image ';'])                                 % Loads frame 
    Photo = ans;                                                % Variable Photo = what the 'Image' contents
    Color = length(size(Photo));                                % If RGB frame, a conversion is required
    if Color > 2
        Photo = Photo(:,:);                                     % Lines up the 3 R, G and B frames
        Photo = Photo(:,1:size(Photo,2)/Color);                 % Selects the first one only (R = G = B = gray)
    end
    Xmin = max(1,floor(Xexpect-1.5*EstRad));                    % Boundaries for static cropping: when several
    Xmax = min(size(Photo,2),floor(Xexpect+1.5*EstRad));        % beads are in the field (and in focus), this
    Ymin = max(1,floor(Yexpect-1.5*EstRad));                    % cropping selects the trapped bead only.
    Ymax = min(size(Photo,1),floor(Yexpect+1.5*EstRad));
    Crop = Photo(Ymin:Ymax,Xmin:Xmax);

    % ---------------------------------------------------------------------
    % Step 1. Edge detection (rough procedure)
    % ---------------------------------------------------------------------
    
    % The first step is to detect a sharp bead and restrict the analysis to its vicinity, therefore ignoring
    % e.g. travelling beads which would appear blurry on the image (the 'if' loop swithches from dynamic (wide)
    % to static (straighter) cropping if several beads are in focus). However, there is no guarantee that the
    % detected points represent the true bead edge (edge will also be detected inside the dark ring and outside
    % the bright halo). This first estimate only allows for ruling out noisy points.
    
    BW = edge(Photo,'canny',[GrLo,GrHi],Filter);                % Detects edges in the full frame
    BWedge = find(BW == 1);                                     % Looking for the corresponding points...
    Xbw = floor(BWedge/size(BW,1)) + 1;                         % ... and their coordinates
    Ybw = ( BWedge/size(BW,1) - floor(BWedge/size(BW,1)) )*size(BW,1);
    [BWXcenter,BWYcenter,R] = circfit(Xbw,Ybw);                 % Fits the edge with a circle
    if R < 1.1*EstRad
        Xmin = max(1,floor(BWXcenter-3*EstRad));                % Boundaries for dynamic cropping: if the edge
        Xmax = min(size(Photo,2),floor(BWXcenter+3*EstRad));    % points can be circled with R < bead radius
        Ymin = max(1,floor(BWYcenter-3*EstRad));                % (+ 10 %) there is no need to restrict the edge
        Ymax = min(size(Photo,1),floor(BWYcenter+3*EstRad));    % tracking to the expected trapping position.
        Crop = Photo(Ymin:Ymax,Xmin:Xmax);
        CXcenter = BWXcenter - Xmin; CYcenter = BWYcenter - Ymin;   % Center coordinates (cropped frame)
        XbwC = Xbw - Xmin; YbwC = Ybw - Ymin;                   % Edge points coordinates (cropped frame)
    else
        BW = edge(Crop,'canny',[GrLo,GrHi],Filter);             % Another edge detection in the cropped image
        BWedge = find(BW == 1);
        XbwC = floor(BWedge/size(BW,1)) + 1;
        YbwC = ( BWedge/size(BW,1) - floor(BWedge/size(BW,1)) )*size(BW,1);
        [CXcenter,CYcenter,R] = circfit(XbwC,YbwC);
    end
    
    % ---------------------------------------------------------------------
    % Step 2. Detection of the dark points located close to the edge
    % ---------------------------------------------------------------------
    
    % The second step consists in detecting the dark points and keeping only those which are inside the detected
    % edge, assumed to be a circle bigger than, or comparable to, the bead.
    
    Black = find(Crop <= Dark);                                 % Looking for the points darker than Dark
    X1 = floor(Black/size(Crop,1)) + 1;                         % Coordinates X and Y in the cropped frame
    Y1 = ( Black/size(Crop,1) - floor(Black/size(Crop,1)) )*size(Crop,1);
    Dist_black = sqrt((X1-CXcenter).^2 + (Y1-CYcenter).^2);     % Distance to the center defined from rough edge
    Black = Black( Dist_black < min(1.1*R,1.3*EstRad) );        % Removes points out of the fitting circle + 10 %
    X1 = floor(Black/size(Crop,1)) + 1;                         % Coordinates X and Y in the cropped frame
    Y1 = ( Black/size(Crop,1) - floor(Black/size(Crop,1)) )*size(Crop,1);
    CXdark = mean(X1); CYdark = mean(Y1);                       % COM of the dark points in the cropped frame
    
    % ---------------------------------------------------------------------
    % Step 3. Edge refining
    % ---------------------------------------------------------------------
    
    % The edge points which correspond to the true bead edge are just outside the dark ring. In this step, the
    % edge points which do not effectively correspond to the physical boundaries of the bead are removed.
    % An 'if' loop is involved to account for the case where no bead is present. If the detected edge
    % effectively corresponds to a true bead, then there must be at least one edge point between 1 and 1.3 times
    % the estimated bead radius from the dark COM, where the true edge is expected to be. Another necessary
    % condition (just in case) is to detect enough dark points (say, more than Area/4). If anyone of these
    % conditions is not satisfied, there is nothing to be detected.
    
    Dist_edge = sqrt((XbwC-CXdark).^2 + (YbwC-CYdark).^2);      % Distance from edge points to the dark COM
    
    if (max(Dist_edge) >= EstRad && min(Dist_edge) <= 1.3*EstRad) && length(Black) > Area/4
        TrueEdge = find( and(Dist_edge >= EstRad, Dist_edge <= 1.3*EstRad) );   % Keeps only the true edge
        if ~isempty(TrueEdge)
            XbwC = XbwC(TrueEdge); YbwC = YbwC(TrueEdge);       % Resets edge points coordinates...
            Dist_edge = sqrt((XbwC-CXdark).^2 + (YbwC-CYdark).^2);  % ... and their distance to the dark COM
            Points = find(Crop > Dark);                         % Looking for all points brighter than Dark
            X2 = floor(Points/size(Crop,1)) + 1;                % Their coordinates
            Y2 = ( Points/size(Crop,1) - floor(Points/size(Crop,1)) )*size(Crop,1);
            Dist = sqrt((X2-CXdark).^2 + (Y2-CYdark).^2);       % Their distance to the dark COM
            Points = Points( Dist < min(Dist_edge) );           % Keeps only those inside the bead edge
            Bead = union(Black,Points);                         % Bead = dark points + other points inside edge
            Xc = floor(Bead/size(Crop,1)) + 1;                  % Coordinates of the corresponding points
            Yc = ( Bead/size(Crop,1) - floor(Bead/size(Crop,1)) )*size(Crop,1);
            X = Xc + Xmin; Y = Yc + Ymin;                       % Coordinates in the full frame
            Xcenter = mean(X); Ycenter = mean(Y);               % Bead center of mass in the full frame
        else
            XbwC = NaN; YbwC = NaN; CXcenter = NaN; CYcenter = NaN;
            X1 = NaN; Y1 = NaN; CXdark = NaN; CYdark = NaN;
            X = NaN; Y = NaN; Xcenter = NaN; Ycenter = NaN;
        end
        
        % To be meaningful, the measured COM should not change too much compared to the averaged COM position.
        % This removes the measured positions found further than 'Tolerance' radii from the average COM position.
        
        if max(X0) > 0
            Xplus = mean(X0(X0>0))+Tolerance*EstRad; Xminus = mean(X0(X0>0))-Tolerance*EstRad;
            Yplus = mean(Y0(Y0>0))+Tolerance*EstRad; Yminus = mean(Y0(Y0>0))-Tolerance*EstRad;
            if (Xcenter > Xplus || Xcenter < Xminus) || (Ycenter > Yplus || Ycenter < Yminus)
                XbwC = NaN; YbwC = NaN; CXcenter = NaN; CYcenter = NaN;
                X1 = NaN; Y1 = NaN; CXdark = NaN; CYdark = NaN;
                X = NaN; Y = NaN; Xcenter = NaN; Ycenter = NaN;
            end
        end
    else
        XbwC = NaN; YbwC = NaN; CXcenter = NaN; CYcenter = NaN;
        X1 = NaN; Y1 = NaN; CXdark = NaN; CYdark = NaN;
        X = NaN; Y = NaN; Xcenter = NaN; Ycenter = NaN;
    end

    % ---------------------------------------------------------------------
    % Optional: Displaying the frame, with the cropped area and the detected points
    % ---------------------------------------------------------------------
    
    if (OneOverTen == 1) && (n/10 == floor(n/10)) || (OneOverTen == 0) && (n < FirstFrame + Display)
        imshow(Photo,'InitialMagnification','fit');
        title(sprintf('Frame number %0.0f',n))
        hold on
        if X1 >= 0          % If something has been detected
            rectangle('position',[Xmin,Ymin,Xmax-Xmin,Ymax-Ymin])
            if DispPoints == 1
                plot(X,Y,'.b')
            end
            if DispDark == 1
                plot(X1+Xmin,Y1+Ymin,'.k', CXdark+Xmin,CYdark+Ymin,'ok')
            end
            if DispEdge == 1
                plot(XbwC+Xmin,YbwC+Ymin,'.w')
            end
            plot(Xcenter,Ycenter,'+r')
        end
        getframe;
    end
    
    % ---------------------------------------------------------------------
    % Coordinates of the bead center (perhaps the most important!)
    % ---------------------------------------------------------------------

    if FirstFrame == 0
        X0(n+1) = Xcenter; Y0(n+1) = Ycenter;
    else
        X0(n) = Xcenter; Y0(n) = Ycenter;
    end
end

% -------------------------------------------------------------------------
% Display results
% -------------------------------------------------------------------------

% ----- Removing meaningless data points
RelevantFrames = find(X0 > 0) - 1;              % Removes the frames where nothing was detected
X0 = X0(X0 > 0); Y0 = Y0(Y0 > 0);               % Removes the 'NaN' points

% ----- Resetting in physical units
Xcom = X0*Calibration; Ycom = Y0*Calibration;
Xshift = Xcom - mean(Xcom); Yshift = Ycom - mean(Ycom);
Center = [Xcom Ycom];
CenterShift = [Xshift Yshift];
Time = RelevantFrames*dt;

% ----- Plotting graphs
figure
subplot(2,1,1)
plot(Time,Xcom)
xlim([min(Time) max(Time)])
ylim([min(Xcom) max(Xcom)])
xlabel('Time (s)')
ylabel('X (um)')
subplot(2,1,2)
plot(Time,Ycom)
xlim([min(Time) max(Time)])
ylim([min(Ycom) max(Ycom)])
xlabel('Time (s)')
ylabel('Y (um)')

figure
plot(Time,Xshift, Time,Yshift)
xlim([min(Time) max(Time)])
ylim([min(min(Xshift,Yshift)) max(max(Xshift,Yshift))])
xlabel('Time (s)')
ylabel('lateral shifts (um)')
legend('X - <X>','Y - <Y>')

figure
plot(Xcom,Ycom)
xlim([min(Xcom) max(Xcom)])
ylim([min(Ycom) max(Ycom)])
xlabel('X (um)')
ylabel('Y (um)')

ToCopy = [mean(Xcom) std(Xcom) mean(Ycom) std(Ycom)];

% ----- Summary: the relevant info to be saved.
Parameters = [Dark;GrHi;GrLo;Filter;Tolerance];
if Saving == 1
    save(NameRoot,'Time','Center','Parameters');
end

toc         % Stops timer