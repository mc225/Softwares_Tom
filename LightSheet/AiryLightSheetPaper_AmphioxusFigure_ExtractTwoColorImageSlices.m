%%AiryLightSheetPaper_AmphioxusFigure_ExtractTwoColorImageSlices
% This code extracts various two color image slices from acquired
% datacubes and saves them as .png files in a subfolder in the directory of
% the datacubes.
% These images are used in the Amphioxus figure in the Airy Light Sheet
% paper. The data used in the paper was Run_4, taken on 13th September
% (CHECK DATE!).
% This code makes use of the subfunction:
% "extractSubVolumeThenProject_2Colors"

all_x=[1:1024];
all_y=[1:2048];
all_z=[1:250];
matfile_data=matfile('recording_lambda488nm_alpha7_beta100');

% geneates images for the "Camera" and standard "Axial" views (x-y (with
% z=0 and x-z (with y=0) respectively) 
[airy_axialView,axial_norm]=extractSubVolumeThenProject_2Colors('deconv',7,100,510:516,all_y,all_z,'Y_transp',NaN);
[airy_cameraView,camera_norm]=extractSubVolumeThenProject_2Colors('deconv',7,100,all_x,all_y,125:127,'N_transp',NaN);
gauss_axialView=extractSubVolumeThenProject_2Colors('record',0,100,510:516,all_y,all_z,'Y_transp',axial_norm);
gauss_cameraView=extractSubVolumeThenProject_2Colors('record',0,100,all_x,all_y,25:27,'N_transp',camera_norm);
bessel_axialView=extractSubVolumeThenProject_2Colors('deconv',0,5,510:516,all_y,all_z,'Y_transp',axial_norm);
bessel_cameraView=extractSubVolumeThenProject_2Colors('deconv',0,5,all_x,all_y,125:127,'N_transp',camera_norm);

% genrates images for the "Other" axial view (y-z) at selected planes
% at x=-40um
[airy_otherView_40,norm_40]=extractSubVolumeThenProject_2Colors('deconv',7,100,all_x,528:532,all_z,'Y_transp',NaN);
gauss_otherView_40=extractSubVolumeThenProject_2Colors('record',0,100,all_x,528:532,all_z,'Y_transp',norm_40);
bessel_otherView_40=extractSubVolumeThenProject_2Colors('deconv',0,5,all_x,528:532,all_z,'Y_transp',norm_40);

% at x=-20um
[airy_otherView_20,norm_20]=extractSubVolumeThenProject_2Colors('deconv',7,100,all_x,768:772,all_z,'Y_transp',NaN);
gauss_otherView_20=extractSubVolumeThenProject_2Colors('record',0,100,all_x,768:772,all_z,'Y_transp',norm_20);
bessel_otherView_20=extractSubVolumeThenProject_2Colors('deconv',0,5,all_x,768:772,all_z,'Y_transp',norm_20);

% at x=0um
[airy_otherView0,norm0]=extractSubVolumeThenProject_2Colors('deconv',7,100,all_x,1008:1012,all_z,'Y_transp',NaN);
gauss_otherView0=extractSubVolumeThenProject_2Colors('record',0,100,all_x,1008:1012,all_z,'Y_transp',norm0);
bessel_otherView0=extractSubVolumeThenProject_2Colors('deconv',0,5,all_x,1008:1012,all_z,'Y_transp',norm0);

% at x=+40um
[airy_otherView40,norm40]=extractSubVolumeThenProject_2Colors('deconv',7,100,all_x,1516:1520,all_z,'Y_transp',NaN);
gauss_otherView40=extractSubVolumeThenProject_2Colors('record',0,100,all_x,1516:1520,all_z,'Y_transp',norm40);
bessel_otherView40=extractSubVolumeThenProject_2Colors('deconv',0,5,all_x,1516:1520,all_z,'Y_transp',norm40);

% check for or make subfolder
if exist('AmphioxusFig')~=7
    mkdir('AmphioxusFig')
end
cd('AmphioxusFig')

% save .pngs
imwrite(airy_cameraView,'airy_cameraView_z=0um.png')
imwrite(airy_axialView,'airy_cameraView_y=0um.png')
imwrite(gauss_cameraView,'airy_cameraView_z=0um.png')
imwrite(gauss_axialView,'airy_cameraView_y=0um.png')
imwrite(bessel_cameraView,'bessel_cameraView_z=0um.png')
imwrite(bessel_axialView,'bessel_cameraView_y=0um.png')

imwrite(airy_otherView_40,'airy_otherView_x=-40um.png')
imwrite(airy_otherView_20,'airy_otherView_x=-20um.png')
imwrite(airy_otherView0,'airy_otherView_x=0um.png')
imwrite(airy_otherView40,'airy_otherView_x=40um.png')
imwrite(gauss_otherView_40,'gauss_otherView_x=-40um.png')
imwrite(gauss_otherView_20,'gauss_otherView_x=-20um.png')
imwrite(gauss_otherView0,'gauss_otherView_x=0um.png')
imwrite(gauss_otherView40,'gauss_otherView_x=40um.png')
imwrite(bessel_otherView_40,'bessel_otherView_x=-40um.png')
imwrite(bessel_otherView_20,'bessel_otherView_x=-20um.png')
imwrite(bessel_otherView0,'bessel_otherView_x=0um.png')
imwrite(bessel_otherView40,'bessel_otherView_x=40um.png')

cd ../