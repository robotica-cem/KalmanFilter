function ObjectTrackingKalman3

close all;
T=1;
Ax=[1 T 0.5*T^2;0 1 T;0 0 0]; Bx=[T^3/6;T^2/2;T]; Cx=[1 0 0]; Fx=eye(3);
Ay=[1 T 0.5*T^2;0 1 T;0 0 0]; By=[T^3/6;T^2/2;T]; Cy=[1 0 0]; Fy=Fx;

errorMed=0; errorEst=0;
Inicio=1;
    
desvW=0.5;  % stdv of acceleration: pixeles/sec^2)
desvV=0.5;
R=desvV^2;
Qx=desvW^2;     Qy=desvW^2;
Px=Fx*Qx*Fx';   Py=Fy*Qy*Fy';  % Inicializando Px y Py;


videoFReader=vision.VideoFileReader('singleball.mp4');
disp('Primero veamos el video original')
while ~isDone(videoFReader)
    Frame2 = step(videoFReader);
    pause(0.1)
    imshow(Frame2);
end
disp('Ahora veamos el video con filtro de Kalman')
pause



videoFReader=vision.VideoFileReader('singleball.mp4');
Frame1 = step(videoFReader); 
while ~isDone(videoFReader)
    Frame2 = step(videoFReader);
    Resta=imabsdiff(Frame1,Frame2);
    BW=im2bw(Resta,0.2);
    s=regionprops(BW,'centroid','area');
    [~,donde]=max(cat(1,s.Area));
    imshow(Frame2)
    hold on
    if Inicio==1
        if ~isempty(donde)
            x=s(donde).Centroid(1);
            y=s(donde).Centroid(2);
            Xest=[x;25;-3];  Yest=[y;-2;0];
            Inicio=2;     
        end;
    end
    if Inicio~=1
        if ~isempty(donde)
            x=s(donde).Centroid(1);
            y=s(donde).Centroid(2);
            y1=x+desvV*randn; y2=y+desvV*randn;
            plot(y1,y2,'r*');
        else
            y1=NaN; y2=NaN;
        end
        L=covarellipse(diag([Px(1,1),Py(1,1)]));
        plot(Xest(1),Yest(1),'g*',L(1,:)+Xest(1),L(2,:)+Yest(1),'g-');
        hold off
        drawnow
        errorEst=errorEst+norm([x;y]-[Xest(1);Yest(1)]);
        errorMed=errorMed+norm([x;y]-[y1;y2]);
        [Xest,Px]=filtroKalman(Ax,Bx,Fx,Cx,Qx,R,Xest,Px,0,y1);
        [Yest,Py]=filtroKalman(Ay,By,Fy,Cy,Qy,R,Yest,Py,0,y2);
    end

    pause(0.1)
end
release(videoFReader);

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