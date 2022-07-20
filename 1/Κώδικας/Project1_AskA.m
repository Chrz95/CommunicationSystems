clc ;
clear all ;
close all ;

%A1

T = 1 ;
over = 10 ;
Ts = T/over ; % Περίοδος δειγματοληψίας
Fs = 1/Ts ; % Συχνότητα δειγματοληψίας
A = 4 ;

[phi1,t_phi] = srrc_pulse(T, Ts, A, 0); % Τα τρία σήματα έχουν τις ίδιες παραμέτρους (T, Ts, A), άρα τα ίδια διανύσμτα χρόνου
[phi2,t_phi] = srrc_pulse(T, Ts, A, 0.5);
[phi3,t_phi] = srrc_pulse(T, Ts, A, 1);

figure(1);
plot(t_phi,phi1);
title('SRRC(t) pulses with different roll-off factors (a) ');
xlabel('time (sec)');
ylabel('Φ(t)');
hold on ;
plot(t_phi,phi2);
plot(t_phi,phi3);
legend ('a = 0','a = 0.5','a = 1');
hold off ;

%A2

Nf = 2048; % Αριθμός δειγμάτων 
f_axis = [-1/2 : 1/Nf: 1/2 - (1/Nf)]; % Ανοικτό διάστημα δεξια , κλειστό αριστερά
F_axis = f_axis*Fs;

PH1F = fftshift(fft(phi1,Nf)*Ts);
PH2F = fftshift(fft(phi2,Nf)*Ts);
PH3F = fftshift(fft(phi3,Nf)*Ts);

figure(2)
plot(F_axis,power(abs(PH1F),2))
title('|Φ(F)|^2:Energy Spectrum of SRRC / Raised Cosine (RC(F)) pulses');
xlabel('frequency (Hz)');
ylabel('|Φ(F)|^2');
hold on
plot(F_axis,power(abs(PH2F),2),'.r')
plot(F_axis,power(abs(PH3F),2),'.g')
legend ('a = 0','a = 0.5','a = 1');
hold off

figure(3)
semilogy(F_axis,power(abs(PH1F),2))
title('Energy Spectrum of SRRC / Raised Cosine (RC(F)) pulses (logarithmic scale)');
xlabel('frequency (Hz)');
ylabel('|Φ(F)|^2');
hold on
semilogy(F_axis,power(abs(PH2F),2),'.r')
semilogy(F_axis,power(abs(PH3F),2),'.g')
legend ('a = 0','a = 0.5','a = 1');
hold off

%A3

figure(4)
semilogy(F_axis,power(abs(PH1F),2))
title('Energy Spectrum of SRRC / Raised Cosine (RC(F)) pulses (logarithmic scale)');
xlabel('frequency (Hz)');
ylabel('|Φ(F)|^2');
hold on
semilogy(F_axis,power(abs(PH2F),2),'.r')
semilogy(F_axis,power(abs(PH3F),2),'.g')

c1 = T/power(10,3)* ones(length(F_axis))  ;
c2 = T/power(10,5)* ones(length(F_axis)) ;
plot (F_axis,c1)
plot (F_axis,c2)

legend ('a = 0','a = 0.5','a = 1','c=T/10^3','c=T/10^5');
hold off

