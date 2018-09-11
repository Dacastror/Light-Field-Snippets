% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT
% Esta función es una implementacion del reenfoque sintetico, basada en el articulo:
% Ng, R., Levoy, M., Brédif, M., Duval, G., Horowitz, M., & Hanrahan, P. (2005). 
% Light field photography with a hand-held plenoptic camera. Computer Science 
% Technical Report CSTR, 2(11), 1-11.

function focusIma = reenfoque(tensor,alpha,method)
    % REENFOQUE construye una imagen reenfocada sintéticamente a partir de un
    %           tensor de campo de luz.
    %
    % focusIma = reenfoque(tensor,alpha,method);
    %
    % Inputs:
    %   tensor: es un areglo de 4 dimensiones
    %
    %   alpha: es un parámetro que define la profundidad del nuevo plano imagen, 
    %   si alpha es igual a 1, el plano imagen coincide con el plano original.
    %   Los valores que puede tomar alpha dependen del campo de visión que posea
    %   el instrumento óptico utilizado para registrar el campo de luz, los valores
    %   usuales de alpha son mayores a -2 y menores a 3.
    %   Para más información ver:
    %   Ng, R., Levoy, M., Brédif, M., Duval, G., Horowitz, M., & Hanrahan, P. (2005). 
    %   Light field photography with a hand-held plenoptic camera. Computer Science 
    %   Technical Report CSTR, 2(11), 1-11. 
    %   o tambien: NG, Ren. Fourierslice photography. 
    %   En ACM transactions on graphics (TOG). ACM, 2005. p. 735-744.
    %   
    %   method: se refiere al método de interpolación a utilizar, puede ser:
    %   'linear' | 'nearest' | 'pchip' | 'cubic' | 'spline' | 'makima'
    %   para mas información consultar:
    %   la.mathworks.com/help/matlab/ref/interpn.html#bt2rb08-1-method
    %
    % Outputs:
    %   focusIma: arreglo bidimensional correspondiente a la imagen reenfocada
    %   sintéticamente sin normalización

    [sx,sy,su,sv] = size(tensor);
    fx = @(xp,u) (alpha-1)*(u-0.5*(su+1)) + xp;  
    fy = @(yp,v) (alpha-1)*(v-0.5*(sv+1)) + yp;
    [xp,u] = ndgrid(1:sx,1:su);
    [yp,v] = ndgrid(1:sy,1:sv);
    x = fx(double(xp),double(u));
    y = fy(double(yp),double(v));
    xn = permute(repmat(x,[1 1 sy sv]),[1 3 2 4]);
    yn = permute(repmat(y,[1 1 sx su]),[3 1 4 2]);
    [xg,yg,u,v] = ndgrid(1:sx,1:sy,1:su,1:sv);
    inter = interpn(xg,yg,u,v,tensor,xn,yn,u,v,method);
    inter(isnan(inter)) = 0;
    focusIma = trapz(trapz(inter,4),3);
end