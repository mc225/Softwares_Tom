%%%
%type=1 bessel
% Row vectors for each mode containing:
%   param(1)=omega  %optical frequency :wavelength=2 pi c/(omega)
%   param(2)=gamma  %cone angle
%   param(3)=pol    %A polaraisation 1:vert; 2:horiz; 3:longi;
%   param(4)=L      %L number
%   param(5)=n      %index of refraction
function prmodes=probemodes(type,num,theta)
    c=299792458;
    switch type
        case 1
            param(1)=1.884e15;  %optical frequency :wavelength=2 pi c/(omega)
            param(2)=.6;  %cone anlge
            param(3)=1;    %A polaraisation 1:vert; 2:horiz; 3:longi;
            param(4)=4;      %L number
            param(5)=1;     %index

            prmodes=[];
            for j0=1:2
                for j2=-0:0
                    for j1=1:num
                        param(2)=(j1-1)/(num)*theta;  
                        param(3)=j0;
                        param(4)=j2;  
                        prmodes=[prmodes;param];
                    end
                end
            end
        case 2
            param(1)=1.884e15;  %optical frequency :wavelength=2 pi c/(omega)
            param(2)=.6;  %cone anlge
            param(3)=2;    %A polaraisation 1:vert; 2:horiz; 
            param(4)=4;      %L number
            param(5)=1;     %index

            prmodes=[];
            for j0=1:1  % pol
                for j2=-0:0 % L
                    for j1=1:num
                    %    param(2)=(j1-1)/(num)*theta;  
                    %    param(2)=asin((j1-1)/(num-1)*sin(theta));  
                        param(2)=asin((j1-.5)/(num-1)*sin(theta));  
                        param(3)=j0;
                        param(4)=j2;  
                        prmodes=[prmodes;param];
                    end
                end
            end
        case 4
            param(1)=2*pi*3e8/0.8e-6;  %optical frequency :wavelength=2 pi c/(omega)
            param(2)=.6;  %cone anlge
            param(3)=2;    %A polaraisation 1:vert; 2:horiz; 
            param(4)=0;      %L number
            param(5)=1;     %index

            prmodes=[];
            for j0=1:1  % pol
                for j2=-0:0 % L
                    for j1=1:num
                    %    param(2)=(j1-1)/(num)*theta;  
                    %    param(2)=asin((j1-1)/(num-1)*sin(theta));  
                        param(2)=asin((j1-.5)/(num-1)*sin(theta));  
                        param(3)=j0;
                        param(4)=j2;  
                        prmodes=[prmodes;param];
                    end
                end
            end

        case 5
            param(1)=2*pi*c/0.8e-6;  %optical frequency :wavelength=2 pi c/(omega)
            param(2)=.6;  %cone angle
            param(3)=2;    %A polaraisation 1:vert; 2:horiz; 
            param(4)=0;      %L number
            param(5)=1;     %index

            prmodes=[];
            for j0=1:2  % pol
                for j2=-2:2 % L
                    for j1=1:num
                    %    param(2)=(j1-1)/(num)*theta;  
                    %    param(2)=asin((j1-1)/(num-1)*sin(theta));  
                        param(2)=asin((j1-.5)/(num-1)*sin(theta));  
                        param(3)=j0;
                        param(4)=j2;  
                        prmodes=[prmodes;param];
                    end
                end
            end
    end
end
