clc ; 
clear all ; 
close all ; 

%C1

N = 100 ; 
b = (sign(randn(N,1))+1)/2; 

%C2

%(α)

X = bits_to_2PAM(b)  ;

%(β)

T = 1 ;
over = 10;
T_s = T/over  ;
X_delta= 1/T_s * upsample(X,over); 
t =  0:T_s:N - T_s ; 

figure(1) ; 
stem (t,X_delta) ; 
title('(t,Xδ(t))');
xlabel('time(sec)');
ylabel('Xδ(t)');

%(γ)

A = 5 ;
a = 0 ;

[phi,t_phi] = srrc_pulse(T, T_s, A, a); % Τα τρία σήματα έχουν τις ίδιες παραμέτρους (T, Ts, A), άρα τα ίδια διανύσμτα χρόνου
X_t = conv(phi,X_delta)*T_s  ; 
t_conv = linspace(t(1)+ t_phi(1),t(end)+ t_phi(end) ,length(X_t)) ;

figure(2) ; 
plot (t_conv,X_t) ;
title('conv(Xδ(t),φ(t))');
xlabel('time(sec)');
ylabel('X(t)');

%(δ)

phi_rev = phi(end:-1:1); % Ανάκλαση φ(-t)
Z_t = conv(X_t,phi_rev) * T_s ; 
t_z = linspace(t_conv(1)+ t_phi(1),t_conv(end)+ t_phi(end) ,length(Z_t)) ;
figure(3);
plot (t_z,Z_t) ;
title('conv(X(t),φ(-t))');
xlabel('time(sec)');
ylabel('Z(t)');
hold on ; 

stem([0:N-1],X);


%C3

%(α)

b = zeros(N/2,2);
b(:,1) = (sign(randn(N/2,1))+1)/2 ;
b(:,2) = (sign(randn(N/2,1))+1)/2 ;
X = bits_to_4PAM(b)  ;

%(β)

T = 1 ;
over = 10;
T_s = T/over  ;
X_delta= 1/T_s * upsample(X,over); 
t =  0:T_s:(N/2) - T_s ; 

figure(4) ; 
stem (t,X_delta) ; 
title('(t,Xδ(t))');
xlabel('time(sec)');
ylabel('Xδ(t)');

X_t = conv(phi,X_delta) * T_s ; 
t_conv = linspace(t(1)+ t_phi(1),t(end)+ t_phi(end) ,length(X_t)) ;

figure(5) ; 
plot (t_conv,X_t) ;
title('conv(Xδ(t),φ(t))');
xlabel('time(sec)');
ylabel('X(t)');

phi_rev = phi(end:-1:1); % Ανάκλαση φ(-t)
Z_t = conv(X_t,phi_rev) * T_s ; 
t_z = linspace(t_conv(1)+ t_phi(1),t_conv(end)+ t_phi(end) ,length(Z_t)) ;

figure(6);
plot (t_z,Z_t) ;
title('conv(X(t),φ(-t))');
xlabel('time(sec)');
ylabel('Z(t)');
hold on ; 

stem([0:((N/2)-1)], X);


