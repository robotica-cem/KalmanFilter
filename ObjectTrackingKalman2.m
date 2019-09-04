function ObjectTrackingKalman2

close all;
T=1;
Ax=[1 T 0.5*T^2;0 1 T;0 0 0]; Bx=[T^3/6;T^2/2;T]; Cx=[1 0 0]; Fx=Bx;
Ay=[1 T 0.5*T^2;0 1 T;0 0 0]; By=[T^3/6;T^2/2;T]; Cy=[1 0 0]; Fy=By;

x=1;y=1; errorMed=0; errorEst=0;
Xest=[x;0;0]; % valor incial del estado estimado
Yest=[y;0;0];
desvW = 0.4; % stdv of acceleration: pixeles/sec^2
desvV = 0.5; % stdv of readings: pixeles
R=desvV^2;
Qx=desvW^2;     Qy=desvW^2;
Px=Fx*Qx*Fx';   Py=Fy*Qy*Fy';  % Inicializando Px y Py;


base_dir = '.\hexbug_frames_compressed\';
cd(base_dir);
f_list =  dir('*png'); % obtiene la lista de las imágenes png en la carpeta.

S_frame = 10; %Obtendremos una imagen promedio de estos frames
[ren,col,cap]=size(imread(f_list(1).name));
imagen=uint16(zeros(ren,col,cap));
for k=1:S_frame
    laimagen=uint16(imread(f_list(k).name));
    imagen = imadd(imagen,laimagen,'uint16');
end
backImagen=uint8(imagen/S_frame);
%{
disp('Primero veamos el video original')
for k=S_frame+1:length(f_list)
    % Esta sección lee la posición de insecto y produce su posición [x,y]
    imagen=imread(f_list(k).name);
    imshow(imagen);
end
disp('Ahora veamos el video con filtro de Kalman')
pause
%}
for k=S_frame+1:length(f_list)
    % Esta sección lee la posición de insecto y produce su posición [x,y]
    imagen=imread(f_list(k).name);
    diffImagen=imabsdiff(backImagen,imagen);
    BW = im2bw(diffImagen,0.15);
    s=regionprops(BW,'centroid','area');
    [~,donde]=max(cat(1,s.Area));
    imshow(imagen)
    hold on
    if ~isempty(donde)
        x=s(donde).Centroid(1);
        y=s(donde).Centroid(2);
    end
    if (k>240)&&(k<250)
        y1=NaN; y2=NaN;
        imshow(backImagen)
    else
        y1=x+desvV*randn; y2=y+desvV*randn;
        plot(y1,y2,'r*'); 
    end
    L=covarellipse(diag([Px(1,1)+eps,Py(1,1)+eps]));
    plot(Xest(1),Yest(1),'g*',L(1,:)+Xest(1),L(2,:)+Yest(1),'g-');
    hold off
    drawnow
    errorEst=errorEst+norm([x;y]-[Xest(1);Yest(1)]);
    errorMed=errorMed+norm([x;y]-[y1;y2]);
    [Xest,Px]=filtroKalman(Ax,Bx,Fx,Cx,Qx,R,Xest,Px,0,y1);
    [Yest,Py]=filtroKalman(Ay,By,Fy,Cy,Qy,R,Yest,Py,0,y2);
end
cd ..
disp(['Error total entre la medición y el verdadero valor ',num2str(errorMed)]);
disp(['Error total entre la estimación y el verdadero valor ',num2str(errorEst)]);
end


function [xk1,Pk1]=filtroKalman(A,B,F,C,Q,R,xk,Pk,uk,yk)
    %Corrección
    if ~isnan(yk)
        Gk=Pk*C'/(C*Pk*C'+R);
        xk=xk+Gk*(yk-C*xk);
        Pk=(eye(size(A))-Gk*C)*Pk;
    end

    %Predicción
    xk1=A*xk+B*uk;
    Pk1=A*Pk*A'+F*Q*F';
end