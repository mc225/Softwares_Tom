%
% A replacement for surf that is more friendly to use.
% Data is converted to double and undersampled if required so the system stays responsive.
%
% Use help surf for more information
%
function res=ssurf(varargin)
    if (min(size(varargin{1}))<=1),
        if (length(varargin{1})>1),
            if nargin>1,
                res=plot(varargin{1},varargin{2});
            else
                res=plot(varargin{1});
            end
        else
            logMessage('A scalar value instead of a matrix was specified: %f',varargin{1});
            res=varargin{1};
        end
        if (nargout==0),
            clear('res');
        end
        return;
    end
    maxPoints=128;
    surfSize=size(varargin{1});
    X=double(varargin{1});
    if (nargin>2),
        Y=double(varargin{2});
        Z=double(varargin{3});
    else        
        Z=double(varargin{1});
        [X,Y]=meshgrid([1:size(Z,2)],[1:size(Z,1)]);
    end
    if (any(size(X)~=size(Z)) || any(size(Y)~=size(Z))),
        logMessage('The X, Y and Z gridsizes should be identical, ignoring X and Y.');
        [X,Y]=meshgrid([1:size(Z,2)],[1:size(Z,1)]);
    end
    %Load the color index or function value in C.
    C=double(varargin{nargin});
    
    if(any(surfSize>maxPoints)),
        %Only reduce the number of points, never increase it.
        if (surfSize(2)>maxPoints),
            xStep=(X(end)-X(1))/(maxPoints-1);
        else
            xStep=X(1,2)-X(1,1);
        end
        if (surfSize(1)>maxPoints),
            yStep=(Y(end)-Y(1))/(maxPoints-1);
        else
            yStep=Y(2,1)-Y(1,1);
        end        
        [sX,sY]=meshgrid([X(1):xStep:X(end)],[Y(1):yStep:Y(end)]);
        Z = interp2(X,Y,Z,sX,sY,'*linear');
        C = interp2(X,Y,C,sX,sY,'*linear');
        X=sX;
        Y=sY;
    end
    
    res=surf(X,Y,single(Z),C);%Convert to single precission to avoid problems with the Axis settings
    if (nargout==0),
        clear('res');
    end
end