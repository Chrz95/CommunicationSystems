clc ; 
clear all ; 
close all ; 

%B4

T = 1 ;
over = 10 ;
Ts = T/over ; % Περίοδος δειγματοληψίας
Fs = 1/Ts ; % Συχνότητα δειγματοληψίας
A = 6;
a = 0.8 ;

first = 1+a/(2*T); % 1/(2*T)
last = (Fs/2) - first;
f0 = (last-first).*rand(1,1) + first;

theta = unifrnd(0-Ts,2*pi,1,1) ; % [0,2π)

[phi,t_phi] = srrc_pulse(T, Ts, A, a); 

Nf = 2048; % Αριθμός δειγμάτων 
f_axis = [-1/2 : 1/Nf: 1/2 - (1/Nf)]; % Ανοικτό διάστημα δεξια , κλειστό αριστερά
F_axis = f_axis*Fs;

PHIF = fftshift(fft(phi,Nf)*Ts);

% (α)

N = 100 ; 
K = 1000 ;

for i=1:K  
    [tx ,X_t] = PAM2_function(N,Ts,over,phi,t_phi) ;
    Ttotal = length(tx);
    Y_t = X_t .* cos(2*pi*f0*tx + theta) ;
    YF = fftshift(fft(Y_t,Nf)*Ts);
    PYFs(i,:) = power(abs(YF),2)/Ttotal ;    
end

SYFexperiment = (sum(PYFs,1)./ K) ; % Είναι η μέση τιμή κάθε στήλης του πίνακα δηλαδή μιας τιμής απο κάθε περιοδόγραμμα

%(β)

Sxftheory = (var(X_t)/T) * power(abs(PHIF),2)*Ts  ; 

F_axis1 = F_axis - f0 ;
F_axis2 = F_axis + f0 ; 

f_axis = [-1/2 : 1/Nf: 1/2 - (1/Nf)]; % Ανοικτό διάστημα δεξια , κλειστό αριστερά
F_axis = f_axis*Fs;

SXFTH_LEFT = 0.25 * Sxftheory ;
SXFTH_RIGHT = 0.25 * Sxftheory ;

% Αφαίρεση των τιμών , όπου τα δύο SxF επικαλύπτονται

k1 = find(F_axis1 >= 0) ;
k2 = find(F_axis2 <= 0) ;

for i=k1 
    SXFTH_LEFT(i) = 0 ; 
end

for i=k2
    SXFTH_RIGHT(i) = 0 ; 
end

% Αφαίρεση των τιμών που τα δύο SxF βγαίνουν εκτός του φάσματος του
% πειραματικου SyF

k3 = find(F_axis1 < min(F_axis)) ;
k4 = find(F_axis2 > max(F_axis)) ;

for i=k3
    SXFTH_LEFT(i) = 0 ; 
end

for i=k4
    SXFTH_RIGHT(i) = 0 ; 
end

figure(1)
semilogy(F_axis,SYFexperiment);
title(['Φασματική Πυκνότητα Ισχύος (Sy(F)) - f0 = ',num2str(f0) , ' Hz']);
xlabel('frequency (Hz)');
ylabel('Sy(F)');
hold on
semilogy(F_axis1,SXFTH_LEFT,'g');
semilogy(F_axis2,SXFTH_RIGHT,'g');
legend ('Πειραματική','Θεωρητική');
hold off ;


