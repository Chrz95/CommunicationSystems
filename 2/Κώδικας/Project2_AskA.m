clc ;
clear all ;
close all ;

%A1

Nf = 2048; % Αριθμός δειγμάτων 

T = 0.01 ;
over = 10 ;
Ts = T/over ; % Περίοδος δειγματοληψίας
Fs = 1/Ts ; % Συχνότητα δειγματοληψίας
A = 4;
a = 0.4 ;

f_axis = [-1/2 : 1/Nf : 1/2 - (1/Nf)]; % Ανοικτό διάστημα δεξια , κλειστό αριστερά
F_axis = f_axis*Fs;

[phi,t_phi] = srrc_pulse(T, Ts, A, a); 
PHIF = fftshift(fft(phi,Nf)*Ts);

figure(1)
semilogy(F_axis,power(abs(PHIF),2))
title('Energy Spectrum of SRRC(t) with parameters (T = 0.01,Ts = 0.001,A = 4,a = 0.4)');
xlabel('frequency (Hz)');
ylabel('|Φ(F)|^2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 1 ;
over = 10 ;
Ts = T/over ; % Περίοδος δειγματοληψίας
Fs = 1/Ts ; % Συχνότητα δειγματοληψίας
A = 6;
a = 0.8 ;

f_axis = [-1/2 : 1/Nf : 1/2 - (1/Nf)]; % Ανοικτό διάστημα δεξια , κλειστό αριστερά
F_axis = f_axis*Fs;

[phi,t_phi] = srrc_pulse(T, Ts, A, a); 
PHIF = fftshift(fft(phi,Nf)*Ts);

figure(2)
semilogy(F_axis,power(abs(PHIF),2))
title('Energy Spectrum of SRRC(t) with parameters (T = 1,Ts = 0.1,A = 6,a = 0.8)');
xlabel('frequency (Hz)');
ylabel('|Φ(F)|^2');

%A2

N = 100 ;

%%%%%%%%%%%% Δημιουργία 2-PAM συνάρτησης %%%%%%%%%%%%%%% 
[tx,X_t] = PAM2_function(N,Ts,over,phi,t_phi) ;

figure(3)
plot(tx,X_t)
title('2-PAM function');
xlabel('time (sec)');
ylabel('X(t) = Sum(Χn*Φ(t-kΤ))');

%A3

%(α)

NUMOFIMPL = 4 ;

for i=1:NUMOFIMPL
    [tx ,X_t] = PAM2_function(N,Ts,over,phi,t_phi) ;
    Ttotal = length(tx);
    XF = fftshift(fft(X_t,Nf)*Ts);

    PXF = power(abs(XF),2)/Ttotal ;

    figure(4)
    subplot(2,2,i)
    plot(F_axis,PXF)
    title (['Περιοδόγραμμα της Χ(t)'])
    legend([num2str(i) ,'η Υλοποιήση'])
    xlabel('frequency (Hz)');
    ylabel('Px(F)');

    figure(5)
    subplot(2,2,i)
    semilogy(F_axis,PXF)
    title ('Περιοδόγραμμα της Χ(t) - (Ημιλογαριθμική Κλίμακα)')
    legend([num2str(i) ,'η Υλοποιήση'])
    xlabel('frequency (Hz)');
    ylabel('Px(F)');
end

%(β), %(γ)

N = 100 ; 
K = 1000 ;

for i=1:K  
   
    [tx ,X_t] = PAM2_function(N,Ts,over,phi,t_phi) ;
    Ttotal = length(tx);
    XF = fftshift(fft(X_t,Nf)*Ts);
    PXFs(i,:) = power(abs(XF),2)/Ttotal ;    
end

SxFapprox = (sum(PXFs,1)./ K) ; % Είναι η μέση τιμή κάθε στήλης του πίνακα δηλαδή μιας τιμής απο κάθε περιοδόγραμμα
Sxftheory = (var(X_t)/T ) * power(abs(PHIF),2)*Ts  ; 

figure(6);
semilogy(F_axis,SxFapprox);
title('Φασματική Πυκνότητα Ισχύος (Sx(F))');
xlabel('frequency (Hz)');
ylabel('Sx(F)');
hold on ;
semilogy(F_axis,Sxftheory);
legend ('Πειραματική','Θεωρητική');
c2 = T/power(10,5)* ones(length(F_axis)) ;
plot (F_axis,c2)
hold off ;

% figure(12);
% plot(F_axis,Sxftheory);
% title('Φασματική Πυκνότητα Ισχύος (Sx(F)) (2-PAM)');
% xlabel('frequency (Hz)');
% ylabel('Sx(F)');

%A4

%(α)

N = 100 ;

%%%%%%%%%%%% Δημιουργία 4-PAM συνάρτησης %%%%%%%%%%%%%%% 
[tx ,X_t] = PAM4_function (N,Ts,over,phi,t_phi) ;

figure(7)
plot(tx,X_t)
title('4-PAM function');
xlabel('time (sec)');
ylabel('X(t) = Sum(Χn*Φ(t-kΤ))');

%(β)

PXF = power(abs(XF),2)/Ttotal ;

N = 100 ; 
K = 1000 ;

for i=1:K  
    % PXFs = zeros (100,length(F_axis)) ;
    [tx ,X_t] = PAM4_function(N,Ts,over,phi,t_phi) ;
    Ttotal = length(tx);
    XF = fftshift(fft(X_t,Nf)*Ts);
    PXFs(i,:) = power(abs(XF),2)/Ttotal ;    
end

SxFapprox = (sum(PXFs,1)./ K) ; 
Sxftheory = (var(X_t)/T ) * power(abs(PHIF),2)*Ts ;

figure(8);
semilogy(F_axis,SxFapprox);
title('Φασματική Πυκνότητα Ισχύος (Sx(F))');
xlabel('frequency (Hz)');
ylabel('Sx(F)');
hold on ;
semilogy(F_axis,Sxftheory);
legend ('Πειραματική','Θεωρητική');
hold off ;

% figure(13);
% plot(F_axis,Sxftheory);
% title('Φασματική Πυκνότητα Ισχύος (Sx(F)) (4-PAM)');
% xlabel('frequency (Hz)');
% ylabel('Sx(F)');

%A5

N = 100 ; 

Tnew = 2*T ;

[phi,t_phi] = srrc_pulse(Tnew, Ts, A, a); 
PHIF = fftshift(fft(phi,Nf)*Ts);

%(α)

for i=1:4
    [tx ,X_t] = PAM2_function(N,Ts,over,phi,t_phi) ;
    Ttotal = length(tx);
    XF = fftshift(fft(X_t,Nf)*Ts);

    PXF = power(abs(XF),2)/Ttotal ;

    figure(9)
    subplot(2,2,i)
    plot(F_axis,PXF)
    title (['Περιοδόγραμμα της Χ(t)'])
    legend([num2str(i) ,'η Υλοποιήση'])
    xlabel('frequency (Hz)');
    ylabel('Px(F)');

    figure(10)
    subplot(2,2,i)
    semilogy(F_axis,PXF)
    title ('Περιοδόγραμμα της Χ(t) - (Ημιλογαριθμική Κλίμακα)')
    legend([num2str(i) ,'η Υλοποιήση'])
    xlabel('frequency (Hz)');
    ylabel('Px(F)');
end

K = 1000 ;

for i=1:K      
    [tx ,X_t] = PAM2_function(N,Ts,over,phi,t_phi) ;
    Ttotal = length(tx);
    XF = fftshift(fft(X_t,Nf)*Ts);
    PXFs(i,:) = power(abs(XF),2)/Ttotal ;    
end

SxFapprox = (sum(PXFs,1)./ K) ; 
Sxftheory = (var(X_t)/Tnew ) * power(abs(PHIF),2)*Ts ;

figure(11);
semilogy(F_axis,SxFapprox);
title('Φασματική Πυκνότητα Ισχύος (Sx(F))');
xlabel('frequency (Hz)');
ylabel('Sx(F)');
hold on ;
semilogy(F_axis,Sxftheory);
legend ('Πειραματική','Θεωρητική');
c2 = T/power(10,5)* ones(length(F_axis)) ;
plot (F_axis,c2)
hold off ;
