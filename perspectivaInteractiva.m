% Copyright (C) 2018 Diego Alberto Castro Rodríguez <dacastror@gmail.com>
% License: MIT, see https://opensource.org/licenses/MIT

function perspectivaInteractiva(tensor)
    [~,~,su,sv] = size(tensor);
    p = [round(su*0.5),round(sv*0.5)];
    selecPerspectiva(p,tensor);
    control = zeros(su,su)+100; % el numero sumado es el color
    figure(6);subplot(1,2,1);imshow(control,[0,255]);
    xlabel('v'); ylabel('u'); title('Control');
    fcn = makeConstrainToRectFcn('impoint',[1,su],[1,sv]);
    pt = impoint(gca,p);
    pt.setColor('r');
    pt.setPositionConstraintFcn(fcn);
    pt.addNewPositionCallback(@(p) observPt(p,tensor));
end

function observPt(p,tensor)
    q = round(p);
    selecPerspectiva(q,tensor);
end

function selecPerspectiva(q,tensor)
    u = q(2); v = q(1); % coordenadas asignadas de forma invertida
    nIma = tensor(:,:,u,v);
    t = strcat('(u, v) = (',num2str(q(1)),{', '},num2str(q(2)),')');
    figure(6); subplot(1,2,2);imshow(nIma,[]); title(t);
    xlabel('x'); ylabel('y');
end
