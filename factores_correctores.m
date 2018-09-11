% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function [ang,campo,periodo,nlx,nly] = factores_correctores(calibra,umbral)
    % FACTORES_CORRECTORES devuelve los factores necesarios para corregir
    % la rotacion y la distorsion en las imagenes de campo de luz.
    % 
    % [ang,campo,periodo,nlx,nly] = factores_correctores(calibra,umbral);
    % 
    % Inputs:
    %   calibra: Imagen de entrada 2d (escala de grises). 'calibra' es una 
    %   imagen de calibracion del campo de luz, esta imagen debe tomarse 
    %   con una iluminacion homogenea de tal manera que todos los
    %   microlentes reciban la misma iluminacion. Debe configurarse la
    %   apertura de la luz entrante tal forma que la imagen generada por
    %   cada microlente no se toque con las imagenes vecinas.
    %
    %   umbral: es un escalar con rango [0, 1]. Es un factor de 
    %   sensibilidad para umbrales adaptativos 
    %
    % Outputs:
    %   ang: angulo de inclinación de las imagenes generadas por los
    %   microlentes respecto a la horizontal
    %
    %   campo: valores de desplazamiento para cada pixel de la imagen que
    %   compensan la aberracion de distorsion, las distancias no enteras
    %   en el registro de las imagenes de los microlentes y la fase 
    %   respecto al borde de la imagen.
    %
    %   periodo: separacion entre cada imagen que generada por los
    %   microlentes
    %
    %   nlx: numero de filas de microlentes detectados
    %
    %   nly: numero de columnas de microlestes detectados
    
    calib = imadjust(calibra);
    [sizex,sizey] = size(calib);
    bw = imbinarize(calib,'adaptive','Sensitivity',umbral);
    periodo_ini = periodo_con_fourier(bw);
    centers = centros_de_microlentes(calib,periodo_ini,umbral);
    size_micro = periodo_promedio(centers,periodo_ini);
    periodo = ceil(size_micro);
    
    ang = angulo_rot_aprox(centers,sizex,sizey,size_micro);
    centers_rot = rotar_centers(centers,ang,[sizex,sizey]*0.5);
    
    nlx = floor(sizex/size_micro);
    nly = floor(sizey/size_micro);
    indis_rect = sortRect(centers_rot,sizex,sizey,size_micro,0.5,0.5);
    
%     p1 = centers_rot(2,:)
%     [size_micro*0.5 size_micro*(0.5+1)]
    
    [gx,gy] = ajuste_a_distorsion(centers_rot,indis_rect,nlx,nly,size_micro);
    
%     Gx = gx(size_micro*0.5, size_micro*(0.5+1));
%     Gy = gy(size_micro*0.5, size_micro*(0.5+1));
%     [Gy+p1(1),Gx+p1(2)]
    
    campo = campo_vectorial_corrector(gx,gy,size_micro,nlx,nly);
%     [campo(1000,1000,1) campo(1000,1000,2)]
   campo = reescalar_campo(campo,size_micro);
%     [campo(1000,1000,1) campo(1000,1000,2)]
%     assignin('base','campo',campo);
end

function campo = reescalar_campo(campo,size_micro)
    [scx,scy,~] = size(campo);
    peri = ceil(size_micro);
    fac = peri/size_micro;
    idealx = (1:scx).';
    idealy = (1:scy);
    scalx = (fac*idealx)-idealx;
    scaly = (fac*idealy)-idealy;
    campo(:,:,1) = campo(:,:,1)-repmat(scaly,scx,1);
    campo(:,:,2) = campo(:,:,2)-repmat(scalx,1,scy);
end

function imag = double2uint16(ima)
    vmax = max(ima(:));
    imag = uint16(ima*65535/vmax);
end

function periodo = periodo_con_fourier(bw)
    porc = 0.01; % porcion de la imagen que abarca el máximo de la tr de Fourier
    trf = abs(fft2(bw));
    [sx,sy] = size(trf);
    trf(1:round(sx*porc),1:round(sy*porc)) = 0;
    trf = trf(1:round(sx/2),1:round(sy/2));
%         imshow(log(trf),[])
    % encontrar las coordenadas del maximo en el espacio de Fourier 2d
    [~,I] = max(trf(:));
    [I_row, I_col] = ind2sub(size(trf),I);
    % periodo de separacion entre microlentes
    periodo = 1/sqrt(((I_col-1)/sy)^2+((I_row-1)/sx)^2);
end

function centers = centros_de_microlentes(ima,periodo_ini,umbral)
    bw = imbinarize(ima,'adaptive','Sensitivity',umbral);
    stats = regionprops('table',bw,'Centroid','MajorAxisLength','MinorAxisLength');
    centers = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2;
    [centers,radii] = discriminar_puntos(centers,radii,periodo_ini);
%     [centers,radii] = refinar_centros(ima,centers,radii,periodo_ini);
% dibujar puntos y radios detectados  
    figure(22); imshow(bw)
    hold on
    plot(centers(:,1),centers(:,2),'+','LineWidth',2,'MarkerEdgeColor','r')
    viscircles(centers,radii);
    hold off
end

function [centros,radios] = discriminar_puntos(centers,radii,periodo)
    % discrimina puntos 3 sigmas por debajo de 
    % la media y mayores a la mitad del periodo
%     figure(12)
%     histogram(radii,500)
    rmax = periodo/2;
    sigma = std(radii(radii<rmax));
    rmin = mean(radii)-3*sigma;
    if rmin<0.75; rmin=0.75; end
    mask = radii>rmin & radii<rmax;
    centros = [centers(mask,1), centers(mask,2)];
    radios = radii(mask); 
end

function periodo = periodo_promedio(centers,periodo_ini)
    [~,D] = rangesearch(centers,centers,2*periodo_ini);
    periodos = zeros(1,length(D));
    for i = 1:length(D)
        dist_vecinos = cell2mat(D(i));
        if length(dist_vecinos)> 2
            periodos(i) = dist_vecinos(2);
        end
    end
    periodos = nonzeros(periodos);
    periodo = mean(periodos);
end

function angulo = angulo_rot_aprox(centers,sx,sy,size_micro)
    centro_ima = [sx,sy]*0.5;
    centers2 = vertcat(centers,centro_ima);
    dist = size_micro*(sqrt(2)+1)*0.5;
    tolerancia = size_micro*(sqrt(2)-1)*0.5;
    [idx,D] = rangesearch(centers2,centers2,dist);
    local_indis = cell2mat(idx(end));
    local_dists = cell2mat(D(end));
    [~,ind] = min(local_dists(2:end));
    prox_indi = local_indis(ind+1);
    pts_recta = zeros(1,2);
    pts_recta(1,:) = centers2(prox_indi,:);
    margen_int = size_micro*1.5;
    num_lens_y = round(0.5*(sy-margen_int*2)/size_micro);
    for i = 1:num_lens_y-1
        local_indis = cell2mat(idx(prox_indi));
        vecinos = centers2(local_indis,:);
        vecs = vecinos-repmat(centers2(prox_indi,:),length(vecinos),1);
        [val,I] = min(abs(size_micro-vecs(:,1)));
        if (val > tolerancia); break; end
        prox_indi = local_indis(I);
        pts_recta(i,:) = centers2(prox_indi,:);
    end
    p = polyfit(pts_recta(:,1),pts_recta(:,2),1);
    angulo = atan(p(1))*180/pi;
end

function centers_rot = rotar_centers(centers,ang,origen)
    ang = -ang*pi/180;
    trT = [1 0 0; 0 1 0; origen(2) origen(1) 1];
    trR = [cos(ang) sin(ang) 0; -sin(ang) cos(ang) 0; 0 0 1];
    tr = (inv(trT))*trR*trT; 
    tr(:,3) = [0;0;1];
    tr = affine2d(tr);
    centers_rot = transformPointsForward(tr,centers);
end

function indis_cen = sortRect(centers,sx,sy,size_micro,cx,cy)
    % cx y cy deben estar en el rango (0,1), definen el punto inicio de
    % busqueda de esta funcion de ordenamiento de los centros de los  
    % microlentes donde cx,cy=0.5 significa iniciar en el centro de la 
    % imagen (0.5 val recomendado)
    dist = size_micro*(sqrt(2)+1)*0.5;
    nx = floor(sx/size_micro);
    ny = floor(sy/size_micro);
    pt_inicio = [sy*cx, sx*cy]; 
    centers2 = vertcat(centers,pt_inicio);
    [idx,D] = rangesearch(centers2,centers2,dist);
    local_indis = cell2mat(idx(end));
    local_dists = cell2mat(D(end));
    [~,ind] = min(local_dists(2:end));
    prox_indi_l = local_indis(ind+1);
    indis_cen = zeros(nx,ny);
    indis_cen = subRect(indis_cen,centers2,idx,prox_indi_l,nx,ny,[1,1],size_micro);
    indis_cen = subRect(indis_cen,centers2,idx,prox_indi_l,nx,ny,[-1,-1],size_micro);
    indis_cen = subRect(indis_cen,centers2,idx,prox_indi_l,nx,ny,[1,-1],size_micro);
    indis_cen = subRect(indis_cen,centers2,idx,prox_indi_l,nx,ny,[-1,1],size_micro);
end

function indis_cen = subRect(indis_cen,centers2,idx,prox_indi_l,nx,ny,dir,size_micro)
    tolerancia = size_micro*(sqrt(2)-1)*0.5;
    fac = (dir+1)*0.5;
    for j = round(nx*0.5):dir(2):1+(nx-1)*fac(2)
        prox_indi = prox_indi_l;
        for i = round(ny*0.5):dir(1):1+(ny-1)*fac(1)
            indis_cen(j,i) = prox_indi;
            local_indis = cell2mat(idx(prox_indi));
            vecinos = centers2(local_indis,:);
            vecs = vecinos-repmat(centers2(prox_indi,:),length(vecinos),1);
            [val,I] = min(abs(size_micro-vecs(:,1)*dir(1)));
            if (val > tolerancia); break; end
            prox_indi = local_indis(I);
        end
        local_indis = cell2mat(idx(prox_indi_l));
        vecinos = centers2(local_indis,:);
        vecs = vecinos-repmat(centers2(prox_indi_l,:),length(vecinos),1);
        [val,I] = min(abs(size_micro-vecs(:,2)*dir(2)));
        if (val > tolerancia); break; end
        prox_indi_l = local_indis(I);
    end
end

function [gx,gy] = ajuste_a_distorsion(centers_rot,indis_rect,nlx,nly,size_micro)
    mascara = indis_rect>0;
    [rx,ry] = resta_con_cuadricula_ideal(centers_rot,indis_rect,nlx,nly,size_micro);
    Cx = rx(mascara);
    Cy = ry(mascara);
    [l, k] = ind2sub( [nlx,nly], find(mascara));
    l = l*size_micro;
    k = k*size_micro;
    gx = fit( [l, k], Cx, 'poly32','Normalize','on'); %32
    gy = fit( [l, k], Cy, 'poly23','Normalize','on');  %23
    % mostrar ajustes de distorsion
    figure(24); plot( gx, [l, k], Cx);
    figure(25); plot( gy, [l, k], Cy);
end

 function [rx,ry] = resta_con_cuadricula_ideal(centers,indis_rect,nlx,nly,size_micro)
    % cuadricula ideal
%     peri = ceil(size_micro);
    peri = size_micro;
    lx = (peri/2:peri:peri*(2*nlx-1)/2).';
    ly = (peri/2:peri:peri*(2*nly-1)/2);
    lx = repmat(lx,1,nly);
    ly = repmat(ly,nlx,1);
    % resta de cuadricula ideal y posiciones reales
    rx = zeros(nlx,nly);
    ry = zeros(nlx,nly);
    for jy=1:nly
       for ix=1:nlx
          if indis_rect(ix,jy)>0
             rx(ix,jy) = lx(ix,jy)-centers(indis_rect(ix,jy),2);
             ry(ix,jy) = ly(ix,jy)-centers(indis_rect(ix,jy),1);
          end
       end
    end
 end

function campo = campo_vectorial_corrector(gx,gy,size_micro,nlx,nly)
    peri = ceil(size_micro);
    fact = size_micro/peri;
    fact = 1;
   [L,K] = meshgrid((1:nlx*peri)*fact,(1:nly*peri)*fact);
    Gx = gx(L, K);
    Gy = gy(L, K);
    campo = zeros(nlx*peri,nly*peri,2);
    campo(:,:,2) = -Gx.';
    campo(:,:,1) = -Gy.';
end

function E = Gaussiana2D(x,xdata)
    E = x(1)*exp(-((xdata(:,:,1)-x(2)).^2+(xdata(:,:,2)-x(3)).^2)/(2*x(4)^2));
end

function [centros,radios] = refinar_centros(ima,centers,radii,size_micro)
    h_s = round(size_micro/2);
    cR = round(centers);
    centros = zeros(1,2);
    radios = zeros(1,1);
    opts = optimset('Display','off');
    for i = 1:length(cR)
        [sx,sy] = size(ima);
        minX = max([1,cR(i,2)-h_s]);
        maxX = min([sx,cR(i,2)+h_s]);
        minY = max([1,cR(i,1)-h_s]);
        maxY = min([sy,cR(i,1)+h_s]);
        mini = min([cR(i,2)-minX, maxX-cR(i,2), cR(i,1)-minY, maxY-cR(i,1)]);
        [X,Y] = meshgrid(-floor(mini):floor(mini));
        xdata = zeros(size(X,1),size(Y,2),2);
        xdata(:,:,1) = X;
        xdata(:,:,2) = Y;
        Z = double(ima(cR(i,2)-mini:cR(i,2)+mini, cR(i,1)-mini:cR(i,1)+mini));
        x0 = [double(ima(cR(i,2),cR(i,1))),0,0,radii(i)*2];
        x = lsqcurvefit(@Gaussiana2D,x0,xdata,Z,[],[],opts);
        centros(i,:)=[centers(i,1)+x(3),centers(i,2)+x(2)];
        radios(i) = x(4);
    end
end
