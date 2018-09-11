% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function verParcelacion(ima,periodo)
    figure(3); imshow(ima,[])
    hold on
    drawVerticalLines(ima,periodo);
    drawHorizontalLines(ima,periodo);
    hold off
end

function drawVerticalLines(ima,periodo)
    [sizex,sizey] = size(ima);
    for i=1:round(sizey/periodo)
        ind = periodo*i+0.5;
        line([ind,ind],[0.5,sizex+0.5],'Color','red')
    end
end

function drawHorizontalLines(ima,periodo)
    [sizex,sizey] = size(ima);
    for i=1:round(sizex/periodo)
        ind = periodo*i+0.5;
        line([0.5,sizey+0.5],[ind,ind],'Color','red')
    end
end