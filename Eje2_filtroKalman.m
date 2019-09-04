function Eje2_filtroKalman
% Definión del modelo
R=30; L=1; C=1e-3; T=0.01;
[A,B,C,~]=ssdata(c2d(ss([0 1;-1/(L*C) -R/L],[0;1/(L*C)],[1 0],0),T));
F=B;
% Definición del ruido
desvW=0.05;  Q=desvW^2;    desvV=0.5; R=desvV^2;
% Definición de los valores iniciales
P=F*Q*F';
Xmed_0=[0;0];    Xest_0=[0;0];   Xnom_0=[0;0];

% Inicio del ciclo de simulación
kf=80; Xmed=zeros(2,kf); Xest=Xmed;  Xnom=Xmed; y=zeros(1,kf);
Xmed(:,1)=Xmed_0; Xest(:,1)=Xest_0;  Xnom(:,1)=Xnom_0; y(1)=0;
for k=1:kf-1
    w=desvW*randn(1,1);    v=desvV*randn(1,1);    
    if  k<40, u=1;  else u=2;  end
    Xnom(:,k+1) = A*Xnom(:,k) + B*u;
    Xmed(:,k+1) = A*Xmed(:,k) + B*u + F*w;
    y(k) = C*Xmed(:,k) + v;
    [Xest(:,k+1),P]=filtroKalman(A,B,F,C,Q,R,Xest(:,k),P,u,y(k));
end

% Graficación de señales
figure(1);
plot(0:kf-1,Xnom(1,:),'-k.',0:kf-1,y,'-r.',0:kf-1,Xest(1,:),'g-');
xlabel('k');
ylabel('Respuesta');
legend('Enominal','Emedido','Xestimado');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xnom(1,:)-y);
disp(['El error cometido entre Enominal y Emedido es: ',...
    num2str(err)])
err = norm(Xnom(1,:)-Xest(1,:));
disp(['El error cometido entre Enominal y Eestimado es: ',...
    num2str(err)])

figure(2);
plot(0:kf-1,Xnom(2,:),'-k.',0:kf-1,Xest(2,:),'-g.');
xlabel('k');
ylabel('Respuesta');
legend('dEnominal','dXestimado');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xnom(2,:)-Xest(2,:));
disp(['El error cometido entre dEnominal y dEestimado es: ',...
    num2str(err)])
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