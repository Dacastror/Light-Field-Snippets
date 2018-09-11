% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function stack = stackReenfoque(tensor,alphas,method)
    % STACKREENFOQUE crea un stack aplicando reenfoque sintético
    % sobre el tensor de campo de luz para los valores de alpha
    % deseados.
    %
    % stack = stackReenfoque(tensor,alphas,method)
    % 
    % Inputs:
    %   tensor: arreglo 4d (x,y,u,v) donde las dos primeras coordenadas
    %   estan relacionadas con las posiciones de las imágenes generadas por los
    %   microlentes y las dos segundas con la posición de los píxeles detras de 
    %   la imagen de un microlente.
    %   
    %   alphas: arreglo 1d con los valores asociados a planos de reenfoque sintéticos
    %   que se desean calcular.
    %
    %   method: string que define el metodo de interpolación a utilizar, 
    %   las opciones son:
    %   'linear'(recomendada)| 'nearest' | 'pchip' | 'cubic' | 'spline' | 'makima'
    %   para más información consultar:
    %   la.mathworks.com/help/matlab/ref/interpn.html#bt2rb08-1-method
    %
    % Outputs:
    %   stack: arreglo 3d que contiene las imagenes reenfocadas sintéticamente,
    %   las primeras dos dimensiones de este son iguales a las dos primeras del
    %   tensor de campo de luz, y la tercera es igual a la longitud del arreglo 
    %   de alphas.

    [sx,sy,~,~] = size(tensor);
    sz = length(alphas);
    tiempos = zeros(2,1);
    stack = zeros(sx,sy,sz);
    count = 1;
    h = waitbar(0,'...','Name','Procesando stack focal...');
    for alpha = alphas
        tic;
        imag = reenfoque(tensor,alpha,method);
        stack(:,:,count) = imag;
        mostrar_reenfoque_actual(imag,alpha);
        tiempos = progreso_del_proceso(h,tiempos,count,sz);
        count = count + 1;
    end
    delete(h);
end

function mostrar_reenfoque_actual(imag,alpha)
    [sx,sy] = size(imag);
    top = 22; % margen para el titulo
    crear_figura('procesando stack',[2,2,2,top],[400,400,sy*2,sx*2+top]);
    imshow(imag,[]); title(sprintf('Reenfoque (alfa = %2.2f)',alpha))
end

function tiempos = progreso_del_proceso(h,tiempos,count,sz)
    tiempos(count) = toc;
    part = count/sz;
    t = round(mean(tiempos)*(sz-count));
    format = 'Tiempo restante %is. (Imagen %i de %i)';
    waitbar(part,h,sprintf(format,t,count,sz));
end

