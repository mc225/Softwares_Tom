classdef Cam < handle
    % Camera super class
    %
    properties (Abstract)
        integrationTime;
        gain;
        regionOfInterest;
        defaultNumberOfFramesToAverage;
    end
    properties (Abstract, SetAccess = private)
        maxSize;
        bitsPerPixel;
    end
    properties
        background=0;
        acquisitionFunctor=[];
    end
    
    methods
        function cam=Cam()
        end
        function cam=acquireBackground(cam,nbFrames)
            if (nargin<2)
                nbFrames=cam.defaultNumberOfFramesToAverage;
            end
            cam.background = cam.acquireDirect(nbFrames);
        end
        function img=acquire(cam,nbFrames)
            if (nargin<2)
                nbFrames=cam.defaultNumberOfFramesToAverage;
            end
            img=cam.acquireDirect(nbFrames)-cam.background;
            if (~isempty(cam.acquisitionFunctor))
                if (ischar(cam.acquisitionFunctor))
                    cam.acquisitionFunctor=str2func(cam.acquisitionFunctor);
                end
                switch (nargin(cam.acquisitionFunctor))
                    case 0
                        cam.acquisitionFunctor();
                    otherwise
                        cam.acquisitionFunctor(img);
                end
            end
        end
        
        function img=acquireMultiple(cam,nbFrames) 
            %new multiple frames into memory, Mingzhou
            if (nargin<2)
                nbFrames=1;
            end
            img=cam.acquireMultipleFrames(nbFrames);
            if (~isempty(cam.acquisitionFunctor))
                if (ischar(cam.acquisitionFunctor))
                    cam.acquisitionFunctor=str2func(cam.acquisitionFunctor);
                end
                switch (nargin(cam.acquisitionFunctor))
                    case 0
                        cam.acquisitionFunctor();
                    otherwise
                        cam.acquisitionFunctor(img);
                end
            end
        end
    end
    methods (Abstract=true, Access = protected)
        img=acquireDirect(cam,nbFrames);        
        img=acquireMultipleFrames(cam,nbFrames);
    end 
end