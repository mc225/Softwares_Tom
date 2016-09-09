% coefficients=zernikeDecomposition(X,Y,Z,maxTerms)
%
% X, Y and Z should be matrices of real numbers.
%
function coefficients=zernikeDecomposition(X,Y,Z,maxTerms)
    R=sqrt(X.^2+Y.^2);
    T=atan2(Y,X);
    invalidPoints=any(isnan(Z(R<=1)));
    for j=1:maxTerms,
        currZernike=zernike(j,R,T);
        if (invalidPoints || j==1),
            normalization=sum(currZernike(~isnan(Z)).^2);
        end
        coefficients(j)=sum(Z(~isnan(Z)).*currZernike(~isnan(Z)))./normalization;
    end
end