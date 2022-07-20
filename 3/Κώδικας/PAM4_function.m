function [tx,Xt]= PAM4_function (X,Ts,over,phi,t_phi) 

X_delta= 1/Ts * upsample(X,over); 
t =  0:Ts:length(X) - Ts ; 

Xt = conv(phi,X_delta) * Ts ; 
tx = linspace(t(1)+ t_phi(1),t(end)+ t_phi(end),length(Xt)) ;
end