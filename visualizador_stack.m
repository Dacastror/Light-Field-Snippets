% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function visualizador_stack(stack,namewindow)
    % VISUALIZADOR_STACK permite ver cada una de las imagenes x,y de un
    % stack x,y,z por medio de una ventana con un slider en su parte
    % inferior.
    %
    % visualizador_stack(stack,namewindow);
    % 
    % Inputs:
    %   stack: arreglo 3d de una secuencia de imagenes
    % 
    %   namewindow: nambre de la ventana donde se visualiza el stack

    % valores por defecto
    size_v = 400;  % tamano de la ventana en direccion vertical
    alto_b = 19;   % tamano del slider en direccion vertical
    alto_t = 42;   % tamano del espacio del titulo en dir vertical
    
    [sx,sy,sz] = size(stack);
    vmax = max(stack(:));
    vmin = min(stack(:));
    ratio = sy/sx;
    fig = crear_figura_stack(alto_b,alto_t,namewindow);
    y = get(gcf,'position');
    set(gcf,'position',[y(1) y(2) ratio*size_v size_v+alto_t])
    set(gca,'units','normalized','position',[0 0.05 1 1])
    imshow(stack(:,:,1),[vmin vmax]);  
    title(sprintf('Frame %d/%d',1,size(stack,3)))
    h_control = uicontrol(fig,'Style', 'slider',...
        'Min',1,'Max',sz,'Value',1,'units', 'normalized');
    addlistener(h_control,'ContinuousValueChange',...
        @(o,e) selectframe(o,e,stack,vmin,vmax));
    tama(fig,0,alto_b,alto_t);
end

function fig = crear_figura_stack(alto_b,alto_t,namewindow)
    array_h = findobj( 'Type', 'Figure', 'Name', namewindow );
    if isempty(array_h)
        fig = figure('Visible','off','Name',namewindow,'NumberTitle',...
            'off','Toolbar','none','SizeChangedFcn',{@tama,alto_b,alto_t}); 
        fig.Visible = 'on';
    else
        fig = array_h(1);
        figure(fig);
        sliders = findobj(fig, 'type', 'uicontrol', 'style', 'slider');
        if not(isempty(sliders)); delete(sliders(1)); end
    end
end

function tama(h_figu,~,alto_barra,alto_titulo)
    fpos = get(h_figu,'Position');
    alto = fpos(4);
    q = alto_barra/alto;
    p = alto_titulo/alto;
    handles = findobj(h_figu, 'type', 'uicontrol', 'style', 'slider');
    h_axes = h_figu.CurrentAxes;
    if not(isempty(handles))
        set(handles(1),'Position', [0 0 1 q]);
    end
    set(h_axes,'units','normalized','Position',[0 q 1 1-p]);
end

function selectframe(h_slider,~,stack,vmin,vmax)
    num = round(get(h_slider,'value'));
    imshow(squeeze(stack(:,:,num)),[vmin,vmax]);
    h_figu = ancestor(h_slider,'figure');
    h_axes = h_figu.CurrentAxes;
    title(h_axes,sprintf('Frame %d/%d',num,size(stack,3)))
end
