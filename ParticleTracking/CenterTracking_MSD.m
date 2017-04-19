% Program for tracking the center of mass of an elliptical particle and calculating its 2D mean square 
% displacement, by detecting the particle edge with a (modified) gradient method.
%
% MRDSV, University of St Andrews, 17-01-2013 (v.1)
% -------------------------------------------------------------------------

clear all; close all; clc; warning off all
tic         % Starts timer     

% -------------------------------------------------------------------------
% Variables
% -------------------------------------------------------------------------

% ----- Variables related to the files to be analyzed
cd 'G:\Matthieu\Videos\diffusion'           % Directory where the image sequence can be found
Date = '2012-12-13';        % Sequence name = 'Date_Time_pXX' with XX = laser intensity (1/100 A)
Time = '18-05-43';
Power = 43;
NameRoot = [eval('Date') '_' eval('Time') '_p' num2str(Power)];
eval(['cd ' Date])
eval(['addpath ' NameRoot]);
FirstFrame = 0;
LastFrame = 500;

% ----- Physical variables
dt = 1/25;              % Time interval between 2 consecutive frames (= 1/frame rate, known at ± 1 fps)
LongAxis = 120;         % Particle long axis (pixels)
ShortAxis = 100;        % Particle short axis (pixels)
LongTol = 1.1;          % Tolerances for edge detection (relative to long & short axis)
ShortTol = 0.9;
Calibration = 20/465;   % Pixel size (µm)

% ----- Edge detection variables
GrHi = 0.20;            % Higher threshold for the edge detection (Canny method)
GrLo = 0.05;            % Lower threshold for the edge detection (Canny method)
Filter = 02;            % Radius of the Gaussian filter for edge detection (Canny method)
GrHi2 = 0.20;           % Higher threshold for the edge gap filling (Canny method)
GrLo2 = 0.05;           % Lower threshold for the edge gap filling (Canny method)
Filter2 = 02;           % Radius of the Gaussian filter for edge gap filling (Canny method)

% ----- Operating preferences
OneOverTen = 0;         % To display one frame every 10 (binary)
Display = 010;          % Alternatively, last displayed frame
DispEdge = 1;           % To display edge points (binary)
Saving = 0;             % To save the measured data (overwrites the previous saving) (binary)

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------

% ----- Definition of recurring variables

Frames = (FirstFrame:LastFrame)';
SizeEdge = zeros(size(Frames)); SizeX = zeros(size(Frames)); correct = zeros(size(Frames));
Rmax = zeros(size(Frames)); Rmin = zeros(size(Frames)); 
Xcenter = []; Ycenter = []; Step = []; t = []; MSD = 0;

for n = FirstFrame:LastFrame
    n
    if FirstFrame == 0
        m = n + 1;
    else
        m = n - FirstFrame + 1;
    end
    
    % ---------------------------------------------------------------------
    % Step 0. Reconstructing frame name, then loading frame
    % ---------------------------------------------------------------------

    if n < 10 
        Image = [NameRoot '_000' num2str(n) '.bmp'];
    else
        if n < 100
            Image = [NameRoot '_00' num2str(n) '.bmp'];
        else
            if n < 1000
                Image = [NameRoot '_0' num2str(n) '.bmp'];    
            else
                if n > 999
                    Image = [NameRoot '_' num2str(n) '.bmp'];
                end
            end
        end
    end
    eval(['imread ' Image ';'])                                 % Load frame 
    Photo = ans;                                                % Variable Photo = what the 'Image' contents
    Color = length(size(Photo));                                % If RGB frame, a conversion is required
    if Color > 2
        Photo = Photo(:,:);                                     % Lines up the 3 R, G and B frames
        Photo = Photo(:,1:size(Photo,2)/Color);                 % Select the first one only (R = G = B = gray)
    end

    % ---------------------------------------------------------------------
    % Step 1. Edge detection (rough procedure)
    % ---------------------------------------------------------------------
    
    % The first step is to detect all the regions of high contrast, which include the particle boundary, and also
    % artifacts inside (curved surface, roughness) and outside (halo) the particle. This gives a first estimate 
    % of the particle centre and radius.
    
    BW = edge(Photo,'canny',[GrLo,GrHi],Filter);                % Detect edges in the full frame
    BWedge = find(BW == 1);                                     % Looking for the corresponding points…
    Xbw = floor(BWedge/size(BW,1)) + 1;                         % … and their coordinates
    Ybw = round((BWedge/size(BW,1) - floor(BWedge/size(BW,1)))*size(BW,1));
    BWXcenter = mean(Xbw); BWYcenter = mean(Ybw);
    R = sqrt((Xbw-BWXcenter).^2 + (Ybw-BWYcenter).^2);
    
    % ---------------------------------------------------------------------
    % Step 2. Detection of the particle boundary (removes external and internal fake edge)
    % ---------------------------------------------------------------------
    
    % The second step consists in removing the dark points which are outside a ring defined by the particle major
    % and minor semiaxes (with tolerances).
    
    ExternalBound = find(R <= LongTol*LongAxis/2); InternalBound = find(R >= ShortTol*ShortAxis/2);
    Edge = intersect(BWedge(ExternalBound), BWedge(InternalBound));
    Xbw = floor(Edge/size(BW,1)) + 1;
    Ybw = round((Edge/size(BW,1) - floor(Edge/size(BW,1)))*size(BW,1));
    SizeEdge(m) = length(Edge);                                 % Number of points kept (to check consistency)
    
    % ---------------------------------------------------------------------
    % Step 3. Finding and filling the gaps in the edge
    % ---------------------------------------------------------------------
    
    % The edge is continuous at a given point A iff exactly 3 points are detected in the 9 adjacents elements
    % (including A itself). If 2 points are detected, A is therefore an edge end. The idea is to connect together
    % 2 consecutive ends, fill the gap by finding the edge in that specific area, and repeat for every gap.
    % The case of 1 isolated point (A detected, without neighbours) is removed as it will likely be overwritten
    % when filling the gaps.
    
    % !!! Care must be taken when applying this step as it may be highly time-consuming (while loop), and it 
    % sometimes creates as many errors as it solves.
    
    % ----- Marking edge ends
    
    BadNeighboring = [];
    for k = 1:length(Edge)
        Neighboring = [Edge(k)-size(BW,1)-1:Edge(k)-size(BW,1)+1, ...
            Edge(k)-1:Edge(k)+1, Edge(k)+size(BW,1)-1:Edge(k)+size(BW,1)+1]';
        Inter = intersect(Edge,Neighboring);
        if length(Inter) == 1
            Edge(k) = NaN; Xbw(k) = NaN; Ybw(k) = NaN;  % Remove potential single points leftover
        end
        if length(Inter) == 2
            BadNeighboring = [BadNeighboring,k];        % Mark edge ends
        end
    end
    
    % ----- Filling gaps
    
    if ~isempty(BadNeighboring)
        EdgeEnds = Edge(BadNeighboring);
        Xends = Xbw(BadNeighboring); Yends = Ybw(BadNeighboring);
        for q = 1:length(EdgeEnds)/2
            XLeft = min(Xends(2*q-1), Xends(2*q)); YTop = min(Yends(2*q-1), Yends(2*q));
            XRight = max(Xends(2*q-1), Xends(2*q)); YBottom = max(Yends(2*q-1), Yends(2*q));
            Width = XRight - XLeft; Height = YBottom - YTop;
            dx = round((14 - Width)/2); dy = round((14 - Height)/2);
            CropFrame = Photo(max(YTop-10,1):min(YBottom+10,size(Photo,1)),...
                max(XLeft-10,1):min(XRight+10,size(Photo,2)));          % A wider area is used (better accuracy)
            NewEdge = edge(CropFrame,'canny',[GrLo2,GrHi2],Filter2);    % Finding edge in this cropped area…
            BWNewEdge = find(NewEdge == 1);
            XbwN = floor(BWNewEdge/size(NewEdge,1)) + 1 + XLeft-10;     % … and the points coordinates
            YbwN = round((BWNewEdge/size(NewEdge,1) - floor(BWNewEdge/size(NewEdge,1)))*size(NewEdge,1)) + YTop-10;
            WhereMatch = zeros(9,1); x = zeros(9,1); y = zeros(9,1);
            for j = 1:9
                x(j) = floor((j-1)/3) - 1; y(j) = (j/3 - floor((j-1)/3))*3 - 2;
                Match = intersect([Xbw Ybw],[XbwN+x(j) YbwN+y(j)],'rows');
                WhereMatch(j) = length(Match);
            end
            J = find(WhereMatch == max(WhereMatch),1);
            XbwN = XbwN + x(J); YbwN = YbwN + y(J);
            Gap = find(XbwN >= XLeft - max(1,dx) & XbwN <= XRight + max(1,dx)...    % Restrict the area where the
                & YbwN >= YTop - max(1,dy) & YbwN <= YBottom + max(1,dy));          % points must be kept
            XbwN = XbwN(Gap); YbwN = YbwN(Gap);
            XY = union([Xbw Ybw],[XbwN YbwN],'rows');       % Put together all the points, without repetition
            Xbw = XY(:,1); Ybw = XY(:,2);
        end

        % ----- Removing the wrong points
        
        LIB = 1;
        if n == FirstFrame
            Test = length(Xbw);
            SizeX(m) = length(Xbw);
        else
            Test = SizeX(m - 1);
            SizeX(m) = SizeX(m - 1);
        end
        while min(LIB) < 5 && Test >= SizeX(m)
            LIB = 5*ones(length(Xbw),1);
            for i = 1:length(Xbw)
                BigSquare = ones(49,2);
                SmallSquare = ones(9,2);
                for ii = 1:49
                    x(ii) = floor((ii-1)/7) - 3; y(ii) = (ii/7 - floor((ii-1)/7))*7 - 4;
                    BigSquare(ii,:) = XY(i,:) + [x(ii) y(ii)];
                end
                SmallSquare = BigSquare([17:19 24:26 31:33],:);
                InterBig = intersect(BigSquare,XY,'rows');
                InterSmall = intersect(SmallSquare,XY,'rows');
                if ~isempty(InterBig)
                    LIB(i) = length(InterBig);
                end
                if length(InterBig) < 5 || length(InterSmall) > 3
                    XY(i) = NaN; Xbw(i) = NaN; Ybw(i) = NaN;
                end
            end
            Test = length(find(Xbw > 0));
        end
    end
    
    Xbw = Xbw(Xbw > 0); Ybw = Ybw(Ybw > 0);             % Remove the NaNs        
    R = sqrt((Xbw-mean(Xbw)).^2 + (Ybw-mean(Ybw)).^2);
    ExternalBound2 = find(R <= LongTol*LongAxis/2); InternalBound2 = find(R >= ShortTol*ShortAxis/2);
    Xbw = Xbw(intersect(ExternalBound2, InternalBound2));
    Ybw = Ybw(intersect(ExternalBound2, InternalBound2));
    SizeX(m) = length(Xbw);                             % Number of points kept (to check consistency)
    
    % After these three steps you'd expect having something you can rely on. Nevertheless, there might still be 
    % errors that cannot be corrected within the present scheme, so the corresponding data points must be ignored 
    % in the MSD (as errors in position are cumulative). I assume the point must be removed if the size of the 
    % boundary differs from the average value by > 5 %.

    if xor(n == FirstFrame, SizeX(m) > 0.95*mean(SizeX(1:m-1)) && SizeX(m) < 1.05*mean(SizeX(1:m-1)))
        correct(m) = 1;
        SizeX(m) = length(Xbw);
        t = [t;n*dt]; Xcenter = [Xcenter;mean(Xbw)]; Ycenter = [Ycenter;mean(Ybw)];
        R = sqrt((Xbw-mean(Xbw)).^2 + (Ybw-mean(Ybw)).^2);
    else
        correct(m) = 0;
        SizeX(m) = SizeX(m - 1);
    end
    Rmax(m) = max(R); Rmin(m) = min(R);
    
    % ---------------------------------------------------------------------
    % Optional: Displaying the frame, with the detected points
    % ---------------------------------------------------------------------
    
    if (OneOverTen == 1) && (n/10 == floor(n/10)) || (OneOverTen == 0) && (n < FirstFrame + Display)
        imshow(Photo,'InitialMagnification','fit');
        title(sprintf('Frame number %0.0f',n))
        hold on
        if DispEdge == 1
            plot(Xbw,Ybw,'.w')
        end
        plot(mean(Xbw),mean(Ybw),'+r')
        getframe;
    end
    
    % ---------------------------------------------------------------------
    % Determination of the MSD
    % ---------------------------------------------------------------------
    
    if correct(m) == 1 && n > FirstFrame
        DeltaX = diff(Xcenter); DeltaY = diff(Ycenter);
        Step = [Step;DeltaX(length(DeltaX))^2 + DeltaY(length(DeltaY))^2];
        MSD = [MSD;sum(Step)];
    else
        Step = Step; MSD = MSD;
    end
    
end

% -------------------------------------------------------------------------
% Display and save results
% -------------------------------------------------------------------------

% ----- Resetting in physical units
X = Xcenter*Calibration; Y = Ycenter*Calibration;
MSD = MSD*Calibration^2;

% ----- Trajectory; MSD = f(lag time)
figure
subplot(1,2,1)
plot(X,Y)
xlabel('X (µm)'); ylabel('Y (µm)')
subplot(1,2,2)
plot(t,MSD)
xlabel('Time (s)'); ylabel('MSD (µm^2)')

% ----- Checking consistency: calculated Rmin and Rmax compared to the initial estimates, and evolution of the 
% edge size (removed points displayed as stars)

figure
subplot(2,1,1)
hold on
title('Routine Robustness Checking')
plot(Frames,Rmax, Frames,Rmin,'r')
plot([FirstFrame,LastFrame],[LongAxis/2,LongAxis/2],'--')
plot([FirstFrame,LastFrame],[ShortAxis/2,ShortAxis/2],':r')
legend('Meas. Rmax','Meas. Rmin','Est. Rmax','Est. Rmin')
xlabel('Frame #')
ylabel('Particle radius (px)')
subplot(2,1,2)
hold on
plot(Frames,SizeEdge, Frames,SizeX)
plot(Frames(correct==0),correct(correct==0)+mean(SizeEdge),'*r')
xlabel('Frame #')
ylabel('Edge size (points)')
legend('before gap filling','after gap filling','removed','location','East')

% ----- Summary: the relevant info to be saved.
Results = [t X Y MSD];
Parameters = [LongAxis;ShortAxis;GrHi;GrLo;Filter;GrHi2;GrLo2;Filter2;LongTol;ShortTol];
if Saving == 1
    save(NameRoot,'Results','Parameters');
end

toc         % Stops timer