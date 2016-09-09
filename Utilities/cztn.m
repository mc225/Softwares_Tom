% f=cztn(x,M,W,A,originCentered)
%
% Calculate the partial spectrum of x using the n-dimensional chirp z transform.
% x: input matrix, centered on central pixel. Implicit zero padding will
% occur at the right and will be corrected for.
% M,W, and A are vectors containing the scalars corresponding to each czt
% per dimension.
% originCentered: boolean, default false. Indicates if the input and output
% matrices have the origin in the central pixel (as defined by fftshift).
% Specify, true, 'centered' or 'originCentered' to set to true. Note that A
% remains the same.
%
% f: output matrix of size [M, size(x,length(M)+1) size(x,length(M)+2) ... size(x,ndims(x))]
%
function f=cztn(x,M,W,A,originCentered)
    if (nargin<2 || isempty(M))
        M = size(x);
    end
    if (nargin<3 || isempty(W))
        W = exp(-2i*pi./M);
    end
    if (nargin<4 || isempty(A))
        A = 1;
    end
    if (nargin<5 || isempty(originCentered))
        originCentered=false;
    end
    if (ischar(originCentered))
        switch(lower(originCentered))
            case 'centered'
            case 'centred'
            case 'origincentered'
            case 'origincentred'
                originCentered=true;
            otherwise
                originCentered=false;
        end     
    end

    nbDims=length(M);
    
    % Work back to deltaRng
    maxFieldSpFreqInInputUnits=floor(size(x)/2);
    halfDeltaRng=maxFieldSpFreqInInputUnits(1:nbDims).*log(W)./(-2i*pi);
    
    % Preshift
    if (originCentered)
        A=A.*W.^(floor(M/2));
    end
    
    f=x;
    for (dimIdx=1:nbDims)
        if (size(f,dimIdx)>1)
            f=permute(f,[dimIdx, 1:dimIdx-1, dimIdx+1:ndims(x)]);
            previousSize=size(f);
            f = czt(f(:,:),M(dimIdx),W(dimIdx),A(dimIdx));
            if (originCentered)
                f=f.*repmat(exp(2i*pi*[0:size(f,1)-1]*halfDeltaRng(dimIdx)).',[1 size(f,2)]); %Correct the pre-shift induced phase error
            end
            f=reshape(f,[size(f,1) previousSize(2:end)]);
            f=ipermute(f,[dimIdx, 1:dimIdx-1, dimIdx+1:ndims(x)]);
        else
            f=repmat(f,[ones(1,dimIdx-1) M(dimIdx) ones(1,ndims(f)-dimIdx)]);
        end
    end
end
