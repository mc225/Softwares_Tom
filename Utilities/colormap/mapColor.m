% imgRGB=mapColor(img,colorMap)
%
% Converts a gray scale image with values between zero and 1 to an RGB
% image using a color map.
%
function imgRGB=mapColor(img,colorMap)
    if (isa(colorMap,'function_handle'))
        colorMap=colorMap(4096);
    end
    if (isinteger(img))
        switch class(img)
            case 'uint16'
                img=double(img)./(2^16-1);
            otherwise
                %Assume 8 bit
                img=double(img)./(2^8-1);
        end
    end
    
    nbColors=size(colorMap,1);
    img=1+floor(nbColors*img);
    img=min(max(1,img),nbColors);
    imgRGB=reshape(cat(3,colorMap(img,1),colorMap(img,2),colorMap(img,3)),[size(img) 3]);
end