function Experimento1
T=0.1; A = [1 0;0 1]; Q = eye(2);
noise=measurenoise(10,T,2);
P=cov(noise');  L=covarellipse(P);
plot(L(1,:),L(2,:),'k-',noise(1,:),noise(2,:),'r*')
hold on
for i=1:6
    P=A*P*A'+Q;     L=covarellipse(P);
    plot(L(1,:),L(2,:),'g-');     pause(0.5)
end
hold off
end

function MeasNoise = measurenoise(tf,T,n)
t = 0 : T : tf;
MeasNoise = randn(n,size(t,2));
end