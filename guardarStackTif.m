% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function guardarStackTif(stack,folder,namestack)
    % GUARDARSTACKTIF permite guardar un arreglo 3d en un unico archivo tif
    % multipagina
    %
    % guardarStackTif(stack,folder,namestack);
    %
    % Inputs:
    %   stack: arreglo 3d qu se desea guardar
    %
    %   folder: ruta relativa o absoluta donde se guardara el archivo
    %
    %   namestack; nombre del archivo que se desea crear sin incluir su
    %   extension.
    
    size_s = size(stack,3);
    v_max = double(max(stack(:)));
    stack16 = uint16((double(stack)/v_max)*65535);
    path_tif = fullfile(folder,strcat(namestack,'.tif'));
    imag = squeeze(stack16(:,:,1));
    imwrite(imag,path_tif)
    for num_ima = 2:size_s
        imag = squeeze(stack16(:,:,num_ima));
        imwrite(imag,path_tif,'WriteMode','append')
    end
end


 

    
    
    
    
    