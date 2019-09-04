function ObjectTrackingKalman

close all;
T=1;
A=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];
B=[0.5*T^2 0;0 0.5*T^2;T 0;0 T];  F=B;
C=[1 0 0 0;0 1 0 0];

x=1;y=1;    errorMed=0; errorEst=0;
Xest= [x; y; 0; 0]; % valor incial del estado estimado
desvW = 1.5; % stdv of acceleration: pixeles/sec^2
desvV=  0.1; % stdv of readings: pixeles
R=diag([desvV^2,desvV^2]);
Q=diag([desvW^2,desvW^2]);

P=F*Q*F';  % Inicializando P;


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
        Y=[NaN;NaN];
        imshow(backImagen)
    else
        Y=[x+desvV*randn; y+desvV*randn];
        plot(Y(1),Y(2),'r*'); 
    end
    L=covarellipse(P(1:2,1:2));
    plot(Xest(1),Xest(2),'g*',L(1,:)+Xest(1),L(2,:)+Xest(2),'g-');
    hold off
    drawnow
    errorEst=errorEst+norm([x;y]-Xest(1:2));
    errorMed=errorMed+norm([x;y]-Y);
    [Xest,P]=filtroKalman(A,B,F,C,Q,R,Xest,P,[0;0],Y);
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