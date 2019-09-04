function Eje4_filtroKalman
% Definición del modelo
T=0.01;
A=[1 T 0.5*T^2;0 1 T;0 0 0]; B=[0;0;1]; C=[1 0 0;0 1 0];
F=B;
% Definición del ruido
desvW=5; Q=desvW^2;
desvV1=5; desvV2=5; R=diag([desvV1^2, desvV2^2]);
% Definición de los valores iniciales
P=F*Q*F';
Xreal_0=[0;0;0];    Xest_0=[0;0;0];   y_0=[0;0];

% Inicio del ciclo de simulación
kf=3000; Xreal=zeros(3,kf); Xest=Xreal; u=zeros(1,kf); y=zeros(2,kf);
Xreal(:,1)=Xreal_0; Xest(:,1)=Xest_0; y(:,1)=y_0;
for k=1:kf-1
    w=desvW*randn(1,1);    v=sqrt(R)*randn(2,1);    
    if (k>500)&&(k<1500), u(k)=10; else u(k)=0; end;
    Xreal(:,k+1) = A*Xreal(:,k) + B*u(k) + F*w;
    y(:,k) = C*Xreal(:,k) + v;
    [Xest(:,k+1),P]=filtroKalman(A,B,F,C,Q,R,Xest(:,k),P,u(k),y(:,k));
end

% Graficación de señales
figure(1);
plot(0:kf-1,Xreal(1,:),'k-',0:kf-1,Xest(1,:),'-g.',0:kf-1,y(1,:),'r:');
xlabel('k');
ylabel('Posicion(m)');
legend('Xreal','Xestimado','Xmedido');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xreal(1,:)-Xest(1,:));
disp(['El error cometido entre Xreal y Xestimado es: ',...
    num2str(err)])

figure(2);
plot(0:kf-1,Xreal(2,:),'k-',0:kf-1,Xest(2,:),'-g.',0:kf-1,y(2,:),'r:');
xlabel('k');
ylabel('Velocidad (m/2)');
legend('dXreal','dXestimado','dXmedido');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xreal(2,:)-Xest(2,:));
disp(['El error cometido entre dXreal y dXestimado es: ',...
    num2str(err)])

figure(3);
plot(0:kf-1,Xreal(3,:),'k-',0:kf-1,Xest(3,:),'g-');
xlabel('k');
ylabel('Aceleración (m/s^2)');
legend('ddXreal','ddXestimado');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xreal(3,:)-Xest(3,:));
disp(['El error cometido entre ddXreal y ddXestimado es: ',...
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