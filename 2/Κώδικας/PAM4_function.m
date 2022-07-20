function [tx,Xt]= PAM4_function (N,Ts,over,phi,t_phi) 

b = zeros(N/2,2);
b(:,1) = (sign(randn(N/2,1))+1)/2 ;
b(:,2) = (sign(randn(N/2,1))+1)/2 ;
X = bits_to_4PAM(b)  ;

X_delta= 1/Ts * upsample(X,over); 
t =  0:Ts:(N/2) - Ts ; 

Xt = conv(phi,X_delta) * Ts ; 
tx = linspace(t(1)+ t_phi(1),t(end)+ t_phi(end) ,length(Xt)) ;
