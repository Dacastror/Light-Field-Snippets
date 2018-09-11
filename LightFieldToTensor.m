% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function tensor = LightFieldToTensor2(ima,periodo)

    % LIGHTFIELDTOTENSOR transforma un arreglo 2d en un arreglo 4d 
    % o en otras palabras, convierte una imagen de campo de luz en un
    % tensor de campo de luz
    %
    % tensor = LightFieldToTensor(ima,periodo);
    %
    % Inputs:
    %   ima: arreglo 2d que contiene al campo de luz
    %
    %   periodo: distancia de separación entre las imagenes
    %   generedas por el arreglo de microlentes en unidades de píxeles,
    %   debe ser un valor entero.
    %
    % Outputs:
    %   tensor: arreglo 4d que representa el tensor de campo de luz
    
    p = periodo; 
    [sx,sy] = size(ima);
    nlx = floor(sx/periodo);
    nly = floor(sy/periodo);
    tensor = zeros(nlx,nly,p,p);
    
    for j=1:nly
        for i=1:nlx
            tensor(i,j,:,:) = ima(1+(i-1)*p:i*p,1+(j-1)*p:j*p);
        end
    end
end

