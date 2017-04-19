% Z=zernikeComposition(X,Y,coefficients)
%
% X and Y should be matrices of real numbers.
% coefficients: the vector of standard zernike coefficients, the first of
%               which are: piston,
%                          tip(x),tilt(y),
%                          defocus, astigmatism-diag,astigmatism-x,
%                          coma-y,coma-x,  trefoil-y,trefoil-x,
%                          spherical aberration
%               where the postscripts indicate the position of the extreme
%               value on the pupil edge.
%
% See also: zernikeDecomposition.m, and zernike.m
%
function Z=zernikeComposition(X,Y,coefficients)
    Z=zeros(size(X));
    R=sqrt(X.^2+Y.^2);
    T=atan2(Y,X);
    for j=1:length(coefficients),
        currZernike=zernike(j,R,T);
        if (j==1),
%            normalization=sum(currZernike(:).^2);          
%            coefficients=coefficients./normalization;% already normalized analytically for a rectangular grid fitted around the pupil
        end
        Z=Z+coefficients(j)*currZernike;
    end
end