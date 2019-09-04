function Eje3_filtroKalman
% Definión del modelo
A=[0.81464 0.08174;-0.0817 0.81464]; T=0.05;
B=[0.00219;0.04525]; F=[0.04306;-0.0474]; C=[1 0];
 
% Definición del ruido
desvW=0.25; Q=desvW^2;    desvV=0.09; R=desvV^2;
% Definición de los valores iniciales
P=F*Q*F';
Xreal_0=[0;0];    Xest_0=[0;0];
 
% Inicio del ciclo de simulación
kf=200; Xreal=zeros(2,kf); Xest=Xreal;
Xreal(:,1)=Xreal_0; Xest(:,1)=Xest_0;
for k=1:kf-1
    w=desvW*randn(1,1);    v=desvV*randn(1,1);    
    u=sin(k*T);
    Xreal(:,k+1) = A*Xreal(:,k) + B*u + F*w;
    y = C*Xreal(:,k) + v;
    [Xest(:,k+1),P]=filtroKalman(A,B,F,C,Q,R,Xest(:,k),P,u,y);
end
 
% Graficación de señales
figure(1);
plot(0:kf-1,Xreal(1,:),'k-',0:kf-1,Xest(1,:),'g-');
xlabel('k');
ylabel('Respuesta');
legend('X1real','X1estimado');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xreal(1,:)-Xest(1,:));
disp(['El error cometido entre X1real y X1estimado es: ',...
    num2str(err)])
 
figure(2);
plot(0:kf-1,Xreal(2,:),'k-',0:kf-1,Xest(2,:),'g-');
xlabel('k');
ylabel('Respuesta');
legend('X2real','X2estimado');
title('Simulación de un sistema dinámico con Filtro de Kalman');
err = norm(Xreal(2,:)-Xest(2,:));
disp(['El error cometido entre X2real y X2estimado es: ',...
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