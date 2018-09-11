% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function fig = crear_figura(varargin)
    % CREAR_FIGURA hace una ventana con nombre en lugar de numero
    % reemplazando una figura con el mismo nombre previamente creada y
    % eliminando las margenes por defecto en la imagen mostrada.
    %
    % fig = crear_figura();
    % fig = crear_figura(namewindow);
    % fig = crear_figura(namewindow,margenes);
    % fig = crear_figura(namewindow,margenes,posiciones);
    %
    % Inputs:
    %   namewindow: (opcional) nombre de la ventana de la figura creada.
    %   
    %   margenes: (opcioneal) arreglo 1d de la forma [left down right top]
    %   con las distancias minimas entre los ejes de la figura y el borde 
    %   de la imagen.
    %
    %   posiciones: (opcional) arreglo 1d de la forma [x1 y1 x2 y2] o 
    %   [x2, y2]. x1,y1 es la posicion en pixeles de la ventana en   
    %   direccion horizontal y en direccion vertical. x2,y2 es el ancho y 
    %   el alto de la ventana.
    %
    % Outputs;
    %   fig: handle de la figura creada.
    
    x2 = 400;  % ancho de la figura por defecto
    y2 = 400;  % alto de la figura por defecto
    margenes = [0,0,0,0];
    

    if nargin > 3
        error('Too many input arguments');
    end
    switch nargin
        case 0                                         
            namewindow = 'Figura';
        case 1
            namewindow = varargin{1};
        case 2                                                                  
            namewindow = varargin{1};
            margenes = varargin{2};
        case 3
            namewindow = varargin{1};
            margenes = varargin{2};
            posiciones = varargin{3};
    end
    
    array_h = findobj( 'Type', 'Figure', 'Name', namewindow );
    if isempty(array_h)
        fig = figure('Visible','off','Name',namewindow,'NumberTitle',...
            'off','SizeChangedFcn',{@tama,margenes});
        if nargin < 3; posiciones = tamano_fig(fig,x2,y2); end
        posiciones = colocar_posiciones(posiciones,fig);
        set(fig,'position',posiciones)
        m = normalizar_margenes(margenes,posiciones);
        set(gca,'units','normalized','position',[m(1) m(2) 1-m(3)-m(1) 1-m(4)-m(2)])
        fig.Visible = 'on';
    else
        fig = array_h(1);
%         if nargin < 3; posiciones = tamano_fig(fig,x2,y2); end
%         posiciones = colocar_posiciones(posiciones,fig);
% %         set(fig,'position',posiciones)
%         m = normalizar_margenes(margenes,posiciones);
%         set(gca,'position',[m(1) m(2) 1-m(3)-m(1) 1-m(4)-m(2)])
        figure(fig);
    end
end

function posiciones = tamano_fig(fig,x2,y2)
    xy = get(fig,'position');
    posiciones = [xy(1) xy(2) x2 y2];
end

function posiciones = colocar_posiciones(posiciones,fig)
    if length(posiciones)==2
        posiciones = tamano_fig(fig,posiciones(1),posiciones(2)); 
    end
end

function m = normalizar_margenes(margen,posiciones)
    x2 = posiciones(3);
    y2 = posiciones(4);
    m = [margen(1)/x2 margen(2)/y2 margen(3)/x2 margen(4)/y2];
end

function tama(h_figu,~,margenes)
    fpos = get(h_figu,'Position');
    L = margenes(1)/fpos(3);
    D = margenes(2)/fpos(4);
    R = margenes(3)/fpos(3);
    T = margenes(4)/fpos(4);
    h_axes = h_figu.CurrentAxes;
    set(h_axes,'units','normalized','Position',[L D 1-R-L 1-T-D]);
end



