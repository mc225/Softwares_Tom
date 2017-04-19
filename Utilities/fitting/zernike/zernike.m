%
% WARNING: returns two complementary zernike polynomials at the same time as a complex function!
%
% result=zernike(m,n,rho,theta)
%    For m>=0, returns the even Zernike polynomial(cos) value as the real part,
%    and the odd polynomial(sin) value as the imaginary part. For m<0, the odd
%    Zernike value is returned as the real part, and the even is returned
%    as the imaginary part.
%
% result=zernike(j,rho,theta)
%    Returns the Zernike polynomial with standard coefficient j
%       "Zernike polynomials and atmospheric turbulence", Robert J. Noll,
%       JOSA, Vol. 66, Issue 3, pp. 207-211, doi:10.1364/JOSA.66.000207
%     The first of which are: piston,
%                             tip(x),tilt(y),
%                             defocus, astigmatism-diag,astigmatism-x,
%                             coma-y,coma-x,  trefoil-y,trefoil-x,
%                             spherical aberration
%                  where the postscripts indicate the position of the extreme
%                  value on the pupil edge.
%
%This function can handle rho and theta vectors but not m and n vectors.
%
%
% Example:
%     grid=[64 64];
%     uRange=-1:(2/(grid(2)-1)):1;
%     vRange=-1:(2/(grid(1)-1)):1;
%     [X,Y]=meshgrid(uRange,vRange);
%     R=sqrt(X.^2+Y.^2);
%     T=atan2(Y,X);
%     ssurf(zernike(4,R,T));
function result=zernike(m,n,rho,theta)
    if (nargin<4),
        theta=rho;
        rho=n;
        j=m;%Standard Zernike coefficient number: 
        n=ceil((sqrt(1+8*j)-1)/2)-1;
        m=j-n*(n+1)/2-1;
        m=m+mod(n+m,2);
    end
    normalization=sqrt((2*(n+1)/(1+(m==0))));%Make orthonormal basis on unit disk (for 2-unit square, multiply with 4/pi)
    result=normalization*zernikeR(abs(m),n,rho).*exp(1i*abs(m)*theta);
    if (nargin<4),
        if (m==0 || mod(j,2)==0),%The first column's coefficients don't have sine or cosine
            result=real(result);
        else
            result=imag(result);
        end
    end
    if (m<0),
        result=conj(result)*1i;
    end
end

function result=zernikeR(m,n,rho)
    result=0;
    if (mod(n-m,2)==0),
        rhoPow=(rho<=1).*rho.^m;
        rhoSqd=(rho<=1).*rho.*rho;
        for l=((n-m)/2):-1:0,
            subResult=(((-1)^l)*faculty(n-l))/(faculty(l)*faculty((n+m)/2-l)*faculty((n-m)/2-l));
            %For speedup: rhoPow=rho.^(n-2*l);
            subResult=subResult .* rhoPow;
            result=result+subResult;
            rhoPow=rhoPow.*rhoSqd;
        end
    end
end

function result=faculty(n)
    result=gamma(n+1);
end