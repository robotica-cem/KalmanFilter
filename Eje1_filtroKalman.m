function Eje1_filtroKalman
% Definión del modelo
A=1; B=1; F=1; C=1;
% Definición del ruido
desvW=0.5;      Q=desvW^2*eye(size(A));  desvV=0.2;      R=desvV^2;
% Definición de los valores iniciales
P=100*eye(size(A));
Xreal_0=30;    Xest_0=20;

% Inicio del ciclo de simulación
kf=50; Xreal=zeros(1,kf); Xest=Xreal; 
Xreal(1)=Xreal_0; Xest(1)=Xest_0;

z=zeros(kf);   z(45)=100;
for k=1:kf-1
    w=desvW*randn(1);    v=desvV*randn(1);    u=1;
    Xreal(k+1) = A*Xreal(k) + B*u + w;
    y = C*Xreal(k) + v  + z(k);
    [Xest(k+1),P]=filtroKalman2(A,B,C,F,Q,R,Xest(k),P,u,y);
end

% Graficación de señales
plot(0:kf-1,Xreal,'-k.',0:kf-1,Xest,'-g.');
xlabel('k');
ylabel('Respuesta');
legend('Xreal','Xest');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xreal-Xest);
disp(['El error cometido entre Xreal y Xest es: ',...
    num2str(err)])
end



function [xk1,Pk1]=filtroKalman1(A,B,F,C,Q,R,xk,Pk,uk,yk)
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


function [xk1,Pk1]=filtroKalman2(A,B,F,C,Q,R,xk,Pk,uk,yk)
    %Corrección
    if ~isnan(yk)       
        Gk=Pk*C'/(C*Pk*C'+R);

        ek=yk-C*xk;
        ICK =eye(size(A))-C*Gk;
        S= ICK'/R*ICK+Gk'/Pk*Gk;
        Sinv=diag(1.0/diag(S));
        se=sign(ek);
        lambda=1;
        z= ek-lambda*Sinv*se;
        z(sign(z) ~= se) = 0;
        xk=xk+Gk*(ek-z);
        
        %xk=xk+Gk*(yk-C*xk);
        Pk=(eye(size(A))-Gk*C)*Pk;
    end

    %Predicción
    xk1=A*xk+B*uk;
    Pk1=A*Pk*A'+F*Q*F';
end