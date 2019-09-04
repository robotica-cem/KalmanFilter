function Robot_1D_FiltroKalman(tf,T)

% Esta función simula el comportamiento de un robot moviéndose en una sola
% dimensión utilizando un sensor de posición ruidoso y una acción de
% control también ruidosa. Se incorpora el uso de un filtro de Kalman para
% estimar de la mejor manera posible tanto la posición como la velocidad
% del robot. 
%
% La función es de la forma:        Robot_1D_FiltroKalman(tf,T)
% donde:
% tf:   tiempo final para la simulación
% T:    periodo de muestreo
close all

vect = 0 : T : tf;              % Este es el vector de tiempo
accel = 1.5;       % Aceleración fija
u = accel*ones(length(vect)); % Vector de aceleraciones en el tiempo

A = [1 T;0 1];          % Matriz de transición de estado
B = [0.5*T^2;T];        % Matriz de control
C = [1 0];              % Matriz de medición
F=B;

deltaV = 10;  % Desv. estándar del ruido de medición de distancia
deltaA = 0.05;  % Desv. estándar del ruido de aceleración
v = deltaV*randn(1,length(vect));   % Vector de ruido en las mediciones
w = deltaA*randn(1,length(vect));   % Vector de ruido en la aceleración
% Nótese que la función randn produce números aleatorios de distribución
% estandar normal con media cero y desviación estandar unitaria.

Q = deltaA^2;      % Covarianza de ruido en el estado
R = deltaV^2  *eye(size(C,1));  % Covarianza de error de medición
P = F*Q*F';   % Covarianza inicial de estimación del estado

x    = zeros(2,length(vect));    % Vector real
xhat = zeros(2,length(vect));    % Vector estimado
y    = zeros(1,length(vect));    % Vector de medidas

for k = 1:length(vect)-1
    % Simulación del proceso que entrega la lectura del sensor
    x(:,k+1) = A*x(:,k) + B*(u(k)+w(k));
    y(k)     = C*x(:,k) + v(k);
    %Uso del filtro de Kalman
    [xhat(:,k+1),P]=filtroKalman(A,B,B,C,Q,R,xhat(:,k),P,u(k),y(k));
    
    figure(1)
    Dx = -10:0.05:90;
    DnormMedi = normpdf(Dx,     y(k),deltaV); 
    DnormEsti = normpdf(Dx,xhat(1,k),sqrt(P(1,1))); 
    DnormReal = normpdf(Dx,   x(1,k),0.05);
    DnormReal = max(DnormEsti)*DnormReal/norm(DnormReal);
    plot(Dx,DnormReal,'k-',Dx,DnormMedi,'b-',Dx,DnormEsti,'g-');
    axis tight
    legend('PposReal','PposMedi','PposEsti');
    xlabel('y');
    ylabel('Probabilidad');
    title(['Gráfica del densidad de probabilidad (',num2str(k),'/',num2str(length(vect)),')']);
    if k < 5
        pause
    else
        pause(0.05)
    end
end

% Graficar los resultados
figure(2)
plot(vect,x(1,:),'-k.',vect,y,'-g.',vect,xhat(1,:),'-r.');
legend('Posición real','Posición medida','Posición estimada');

grid;
xlabel('Tiempo (s)');
ylabel('Posición (m)');
title('Gráfica del filtro de Kalman');
err = norm(x(1,:)-xhat(1,:));
disp(['El error cometido entre la posicion real y la estimada es: ',...
    num2str(err)])
err = norm(x(1,:)-y);
disp(['El error cometido entre la posicion real y la medida es: ',...
    num2str(err)])

figure(3)
plot(vect,x(2,:),'-b.',vect,xhat(2,:),'-m.');
legend('Velocidad real','Velocidad estimada');
grid;
xlabel('Tiempo (s)');
ylabel('Velocidad (m/s)');
title('Gráfica del filtro de Kalman');
err = norm(x(2,:)-xhat(2,:));
disp(['El error cometido entre la velocidad real y la estimado es: ',...
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