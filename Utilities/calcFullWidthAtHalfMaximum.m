% fullWidthAtHalfMaximum=calcFullWidthAtHalfMaximum(X,Y,noiseBackground)
%
% Fits a Gaussian to the data and calculates the full width at half the
% maximum.
% 
% X must be a vector with the sample coordinates or a matrix with the
%     coordinates for a sample set in each column. If X is a vector than it
%     must have the same number as elements as the first dimension in Y.
%     The same coordinates are assumed for all sample sets. If it is a
%     matrix it must have the same size as Y.
% Y must be a vector with the sample values or a matrix with the values of
%     a sample set in each column.
% noiseBackground: (default true), if true, a biased Gaussian is fitted to the data,
% otherwise an unbiased Gaussian is fitted
function fullWidthAtHalfMaximum=calcFullWidthAtHalfMaximum(X,Y,noiseBackground)
    if (nargin<2 || isempty(Y))
        Y=X;
        X=[1:size(Y,1)];
    end
    if (nargin<3)
        noiseBackground=true;
    end

    inputSize=size(Y);
    if (any(size(X)<inputSize))
        X=repmat(X(:),[1 inputSize(2:end)]);
    end
    if (prod(inputSize)>max(inputSize))
        fullWidthAtHalfMaximum=zeros([inputSize(2:end) 1]);
        for curveIdx=1:prod(inputSize(2:end))
            fullWidthAtHalfMaximum(curveIdx)=calcSingleFullWidthAtHalfMaximum(X(:,curveIdx),Y(:,curveIdx),noiseBackground);
        end
    else
        fullWidthAtHalfMaximum=calcSingleFullWidthAtHalfMaximum(X,Y,noiseBackground);
    end
end

function fullWidthAtHalfMaximum=calcSingleFullWidthAtHalfMaximum(X,Y,noiseBackground)
    if (noiseBackground)
        offset=min(Y);
    else
        offset=0;
    end
    magnitude=max(Y)-offset;
    center=mean(X.*(Y-offset))/mean(Y);
    sigma=sqrt(mean(((X-center).^2).*(Y-offset)));
   
    if (noiseBackground)
        x0=double([offset magnitude center sigma]);
        [x]=fminsearch(@(x) norm(gaussian(x(1),x(2),x(3),x(4),X)-Y),x0,optimset('Display','none','TolX',1e-9,'TolFun',1e6*max(Y(:))));
    else
        x0=double([magnitude center sigma]);
        [x]=fminsearch(@(x) norm(gaussian(0,x(1),x(2),x(3),X)-Y),x0,optimset('Display','none','TolX',1e-9,'TolFun',1e6*max(Y(:))));
    end
    x=num2cell(x);
    if (noiseBackground)
        [offset magnitude center sigma]=deal(x{:});
    else
        [magnitude center sigma]=deal(x{:});
    end
    sigma=abs(sigma); % only the absolute value is important really
    
    fullWidthAtHalfMaximum=2*sigma*sqrt(-2*log(.5));
    
%     figure; plot(X,Y,'-',X,gaussian(offset,magnitude,center,sigma,X),':');
end

function Y=gaussian(offset,magnitude,center,sigma,X)
    Y=offset+magnitude*exp(-(X-center).^2/(2*sigma^2));
end