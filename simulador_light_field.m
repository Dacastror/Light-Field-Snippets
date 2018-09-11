% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com> 
% version: 0.3
% License: MIT, see https://opensource.org/licenses/MIT
warning ('off','all');

%%% el stack se ajusta al campo de vision del microscopio simulado %%%

folderInput = 'Stacks-simulados';
folderOutLF = 'campos-de-luz';
folderOutFS = 'reconstrucciones';
nameStack = 'letras_ABC.tif';
nameOut   = 'LF_simulado';
nameCalib = 'calibracion\stack_calibracion_proy_fit.mat';

paso       = 0.3;   % distancia entre planos en el stack de entrada(um)
lentesx    = 101;   % numero de microlentes en direccion vertical
lentesy    = 101;   % numero de microlentes en direccion horizontal
pix_sensor = 21;    % lado en pixeles del cuadrado detras de cada microlete
size_lente = 130;   % diametro de cada microlente del arreglo (um)
magnificacion     = 50;    % magnificacion del sistema optico
apertura_numerica = 1.1;   % apertura numerica del objetivo de microscopio
indice_refraccion = 1.33;  % indice de refraccion del medio de inmersion
% muestreo utilizado para calcular el campo de luz, 1 significa que cada
% imagen de entrada se escala para encajarse en una matriz de tamano 
% = lentesx, lentesy. Esto no altera el tamano final del campo de luz.
muestreo = 3;
centerZ = 0; % desplazar el volumen observado en direccion Z (um)
method = 'linear';
methodScal = 'bilinear';

aplicarMascara = true;  % simular apertura circular
margen = 2;          % margen entre el cuadrado y el circulo inscrito
kernel = [3 3];      % kernel del filtro gaussiano para suavizar bordes
factor_gauss = 0.92; % intensidad del filtro gaussiano aplicado

alphas = -0.5:0.05:2.5;
usarCalibracion = false;
% center = [0,0,0];
ttot = tic;
% informacion del stack
fname = fullfile(folderInput,nameStack);
info  = imfinfo(fname);
sizey = info.Width;
sizex = info.Height;
sizez = numel(info);
% type  = info.BitDepth;

% lenMax: numero de microlentes en la direccion con mas microlentes
lenMax = lentesy;
if lentesx>lentesy; lenMax = lentesx; end
% angulo de apertura del objetivo
angulo = asin(apertura_numerica/indice_refraccion);
distMax = tan(angulo);
% distancia virtual entre puntos de "observacion"
deltaXY = 2*distMax/(pix_sensor-1);
% tamano de imagen en con el que se realizaran las proyecciones
new_size = lenMax*muestreo;
% factor de conversion: micras a pixeles
% um2pix = pix_sensor*size_lente/new_size;
size_new_pix = lenMax*size_lente/(magnificacion*new_size);
% comparacion entre el ratio xy del stack y el de las microlentes
compRatios = sizey/sizex > lentesy/lentesx;
% factor de escala y conversion de micras a pixeles
factor_scal = new_size*magnificacion/(size_lente*lenMax);
% paso entre imagenes en pixeles
paso_scal = paso*factor_scal;
% crear matriz de margen para igualar el ratio del numero de microlentes
% con el ratio del tamano de la imagen
if compRatios
    new_sizex = sizey*lentesx/lentesy;
    margenx = floor((new_sizex-sizex)*0.5);
    margeny = 0;
    matrizMx = zeros(margenx,sizey);
else 
    new_sizey = sizex*lentesy/lentesx;
    margeny = floor((new_sizey-sizey)*0.5);
    margenx = 0;
    matrizMy = zeros(sizex,margeny);
end
% tan(angulos) de proyeccion del volumen respecto al eje optico
tanAng = -distMax:deltaXY:distMax;
% dist del plano focal a la primera imagen del stack de entrada (en pix)
z_ini = -paso_scal*(sizez-1)*0.5+factor_scal*centerZ;
svistaX = lentesx*muestreo;
svistaY = lentesy*muestreo;
persph  = zeros(svistaX,svistaY*pix_sensor);
mosaico = zeros(svistaX*pix_sensor,svistaY*pix_sensor);
mosaic_total = mosaico;

barra = waitbar(0,'...','Name','calculando proyecciones del volumen...');
tiempos = zeros(2,1);
for k = 1:sizez
    t = tic;
    Ima = imread(fname,k);
    % aplicar margenes para igualar ratios si es necesario
    if margenx>0; Ima = [matrizMx; Ima; matrizMx]; end
    if margeny>0; Ima = [matrizMy, Ima, matrizMy]; end
    Ima = imresize(Ima, [svistaX,svistaY],methodScal);
    kz = z_ini+paso_scal*(k-1);
    trs = kz*tanAng;
    for i = 1:pix_sensor
        vista_x = imtranslate(Ima,[trs(i), 0],method);
        persph(:,(i-1)*svistaY+1:i*svistaY) = vista_x;
    end
    for i = 1:pix_sensor
        vistas_y = imtranslate(persph,[0, trs(i)],method);
        mosaico((i-1)*svistaX+1:i*svistaX,:) = vistas_y;
    end
    mosaic_total = mosaic_total + mosaico;
    % para barra de progreso
    tiempos(k) = toc(t);
    part = k/sizez;
    t = round(mean(tiempos)*(sizez-k));
    formato = 'Tiempo restante %is. (procesado %i de %i)';
    waitbar(part,barra,sprintf(formato,t,k,sizez));
end
final_size = [lentesx*pix_sensor,lentesy*pix_sensor];
mosaic_total = imresize(mosaic_total,final_size,methodScal);
delete(barra)
% 
% fig = crear_figura('Mosaico');
% imshow(mosaic_total,[])
% 
lx = lentesx;
ly = lentesy; 
tensor = zeros(lentesx,lentesy,pix_sensor,pix_sensor);
for j=1:pix_sensor
    for i=1:pix_sensor
        tensor(:,:,i,j) = mosaic_total(1+(i-1)*lx:i*lx,1+(j-1)*ly:j*ly);
    end
end

campo = zeros(lentesx*pix_sensor,lentesy*pix_sensor);
L = pix_sensor;
for j = 1:lentesy
    for i = 1:lentesx
        campo(1+(i-1)*L:i*L,1+(j-1)*L:j*L) = tensor(i,j,:,:);
    end
end

% aplicar mascara para simular la apertura del objetivo
if aplicarMascara
    t = linspace(0,2*pi,360);  
    r = (L-margen)/2;              
    circle = poly2mask(r*cos(t)+L/2-0.5, r*sin(t)+L/2-0.5, L, L);
    filtro = fspecial('gaussian', kernel, factor_gauss);
    circle_b = imfilter(im2double(circle), filtro);
    mascara_c = repmat(circle_b,lentesx ,lentesy);
    campo = campo.*mascara_c;
end

toc(ttot)

% Guardar campo de luz
mini = min(campo(:));
maxi = max(campo(:));
campo = uint16(((campo-mini)/(maxi-mini))*2^16-1);
formato = '%s_(paso %1.3fum centro %1.2fum)%s';
nameLF = sprintf(formato,nameOut,paso,centerZ,nameStack);
imwrite(campo,fullfile(folderOutLF,nameLF));
figure(1)
imshow(campo,[])

if usarCalibracion
    load(nameCalib)
    desp = -(sizez-1)*paso/2:paso:(sizez-1)*paso/2;
    alphas = fitAlpha(desp);
    if isempty(alphas); error('archivo de calibracion incorrecto'); end
    alphas = alphas';
end

perspectivaInteractiva(tensor)
stack = stackReenfoque(tensor,alphas,method);
visualizador_stack(stack,'stack')
[~,name,~] = fileparts(nameStack);
formato = '%s_(paso %1.3fum centro %1.2fum)%s';
nameFS = sprintf(formato,'Focus_stack',paso,centerZ,name);
guardarStackTif(stack,folderOutFS,nameFS)






