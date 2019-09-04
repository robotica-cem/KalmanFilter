function [x,y]=dummy_gps_reading(n)
rng('shuffle');
x=11+0.5*randn(n,1);
y=22+2.0*randn(n,1);