function [tx,Xt] = PAM2_function (N,Ts,over,phi,t_phi) 

b = (sign(randn(N,1))+1)/2; 
X = bits_to_2PAM(b)  ;

X_delta= 1/Ts * upsample(X,over); 
t =  0:Ts:N - Ts ; 

Xt = conv(phi,X_delta)*Ts  ; 
tx = linspace(t(1)+ t_phi(1),t(end)+ t_phi(end) ,length(Xt));