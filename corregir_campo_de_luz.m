% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function newima = corregir_campo_de_luz(ima,ang,vecs,periodo,nlx,nly,method)
    [sizex,sizey] = size(ima);
    ima_rot = imrotate(ima,ang,method,'crop');
    newima = zeros(nlx*periodo,nly*periodo);
    newima(1:sizex,1:sizey) = ima_rot;
    newima = imwarp(newima,vecs ,'linear');
end