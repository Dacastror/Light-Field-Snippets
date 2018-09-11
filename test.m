rutaBase = '';
rutaCampo = '';
rutaCalibracion = '';
addpath(rutaBase);

ima = imread(fullfile(rutaBase,rutaCampo));
calibracion = imread(fullfile(rutaBase,rutaCalibracion));

% calibracion = rgb2gray(calibracion);
% calibracion = imrotate(calibracion,1,'bilinear');
% ima = ima(42:end-41,:);
% calibracion = calibracion(42:end-41,:);

% calibra_cor = factores_correctores(calibracion,'bicubic',0.5);
% calibracion = imbinarize(calibracion,'adaptive','Sensitivity',0.4);
% calibracion = uint8(calibracion*255);
% calibracion = calibracion+20000;

[ang,vecs,periodo,nlx,nly] = factores_correctores(calibracion,0.5);

newima = corregir_campo_de_luz(ima,ang,vecs,periodo,nlx,nly,'bilinear');
calibracion = corregir_campo_de_luz(calibracion,ang,vecs,periodo,nlx,nly,'bilinear');



% newima = corregir_campo_de_luz(calibracion,ang,vecs,periodo,nlx,nly,method);
% newima = imresize(newima,2099/2100,'bilinear');
% % temp = newima(1:end-1,1:end-1);

tensor = LightFieldToTensor(newima,periodo);
% perspectivaInteractiva(tensor)
% calibracion = imresize(calibracion,2099/2100,'bilinear');
% calibracion = uint16(calibracion);
% bw = imbinarize(calibracion,'adaptive','Sensitivity',0.5);
verParcelacion(ima,periodo)
alphas = -1:0.05:3;
stack = stackReenfoque(tensor,alphas,'linear');

visualizador_stack(stack,'stack')

guardarStackTif(stack,rutaBase,namestack)
% [sizex,sizey] = size(calibracion);
% method = 'bilinear';
% calibra_rot = imrotate(calibracion,ang,method,'crop');
% newcalibra = zeros(nlx*periodo,nly*periodo);
% newcalibra(1:sizex,1:sizey) = calibra_rot;
% cor = imwarp(newcalibra,campo ,'linear');
% vmax = max(cor(:));
% cor = medfilt2(uint16(cor*65535/vmax));
% cor = imbinarize(cor,'adaptive','Sensitivity',0.5);
% % tensor = LightFieldToTensor(ima,periodo)
% figure(32)
% imshow(cor)

% imshow(cor(:,periodo*20:periodo*21))

% [ima,calibracion,periodo] = rotacion_automatica(ima,calibracion,'bicubic',0.5,21);
% tensor = LightFieldToTensor(ima,periodo);
% % focusIma = reenfoque(tensor,1.5,'linear');
% 
% stack = stackReenfoque(tensor,0:0.2:2,'linear');