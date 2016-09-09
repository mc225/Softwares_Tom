% Z=zernikeComposition(X,Y,coefficients)
%
% X and Y should be matrices of real numbers.
% coefficients is the vector of standard zernike coefficients
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