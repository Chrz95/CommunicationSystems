clear all;
close all;
clc;

% A1

N = 500 ; 
bit_seq = (sign(randn(4*N,1))+1)/2;

% A2

A = 1;
X = bits_to_4_PAM(bit_seq, A);

% A3

XI = X(1:N);
XQ = X(N+1:end);

% A4

T = 1 ;
over = 30 ; 
Ts = T/over ; % Περίοδος δειγματοληψίας
Fs = 1/Ts ; % Συχνότητα δειγματοληψίας

A_phi = 6 ; 
roll_off = 1 ; 

[phi,t_phi] = srrc_pulse(T,Ts,A_phi,roll_off);

[tx,XIt] = PAM4_function (XI,Ts,over,phi,t_phi) ;
[tx,XQt] = PAM4_function (XQ,Ts,over,phi,t_phi) ;

figure(1);
plot (tx,XIt) ; 
title('4-PAM function - XI(t)');
xlabel('time (sec)');
ylabel('XI(t) = Sum(ΧI,n*Φ(t-kΤ))');

figure(2)
plot (tx,XQt) ; 
title('4-PAM function - XQ(t)');
xlabel('time (sec)');
ylabel('XQ(t) = Sum(ΧQ,n*Φ(t-kΤ))');

Nf = 2048 ; % Αριθμός δειγμάτων στο πεδίο της συχνότητας 
f_axis = [-1/2 : 1/Nf : 1/2 - (1/Nf)]; % Ανοικτό διάστημα δεξια , κλειστό αριστερά
F_axis = f_axis*Fs;

Ttotal = length(tx);
XIF = fftshift(fft(XIt,Nf)*Ts);
PXIF = power(abs(XIF),2)/Ttotal ;
figure(3)
semilogy(F_axis,PXIF)
title (['Περιοδόγραμμα της ΧI(t)'])
xlabel('frequency (Hz)');
ylabel('PXI(F)');

XQF = fftshift(fft(XQt,Nf)*Ts);
PXQF = power(abs(XQF),2)/Ttotal ;
figure(4)
semilogy(F_axis,PXQF)
title (['Περιοδόγραμμα της ΧQ(t)'])
xlabel('frequency (Hz)');
ylabel('PXQ(F)');

% A5

f0 = 2;

XImod = XIt.*(2*cos(2*pi*f0*tx)); % Διαμόρφωση
XQmod = XQt.*((-2)*sin(2*pi*f0*tx));

figure(5);
plot (tx,XImod) ; 
title('XImod = 2*XIt.*cos(2*pi*f0*tx)');
xlabel('time (sec)');
ylabel('XImod');

figure(6)
plot (tx,XQmod) ; 
title('XQmod = -2*XQt.*sin(2*pi*f0*tx)');
xlabel('time (sec)');
ylabel('XQmod');

XImodF = fftshift(fft(XImod,Nf)*Ts);
PXImodF = power(abs(XImodF),2)/Ttotal ;
figure(7)
semilogy(F_axis,PXImodF)
title (['Περιοδόγραμμα της XImod(t)'])
xlabel('frequency (Hz)');
ylabel('PXImodF(F)');

XQmodF  = fftshift(fft(XQmod,Nf)*Ts);
PXQmodF = power(abs(XQmodF),2)/Ttotal ;
figure(8)
semilogy(F_axis,PXQmodF)
title (['Περιοδόγραμμα της XQmod(t)'])
xlabel('frequency (Hz)');
ylabel('PXQmodF(F)');

% A6

Xmod = XImod + XQmod;

figure(9)
plot (tx,Xmod) ; 
title('Xmod = XImod + XQmod');
xlabel('time (sec)');
ylabel('Xmod');

XmodF  = fftshift(fft(Xmod,Nf)*Ts);
PXmodF = power(abs(XmodF),2)/Ttotal ;
figure(10)
semilogy(F_axis,PXmodF)
title (['Περιοδόγραμμα της Xmod(t)'])
xlabel('frequency (Hz)');
ylabel('PXQmodF(F)');

% A7 , A8

mu = 0 ;
SNRdb = 20 ; 
% Δημιουργία του θορύβου
varianceW = (10*A)/(2*Ts*power(10,SNRdb/10)); 
sigmaW = sqrt(varianceW) ; 
Noise = normrnd(mu,sigmaW,1,length(tx)) ;

XmodNoise = Xmod + Noise ; 

% A9

XIunmod = XmodNoise .* cos(2*pi*f0*tx) ;
XQunmod = -XmodNoise .* sin(2*pi*f0*tx) ; 

figure(11);
plot (tx,XIunmod) ; 
title('XImod * cos(2*pi*f0*tx)');
xlabel('time (sec)');
ylabel('XIunmod');

figure(12)
plot (tx,XQunmod) ; 
title('XQmod * -sin(2*pi*f0*tx)');
xlabel('time (sec)');
ylabel('XQunmod');

XIunmodF = fftshift(fft(XIunmod,Nf)*Ts);
PXIunmodF = power(abs(XIunmodF),2)/Ttotal ;
figure(13)
semilogy(F_axis,PXIunmodF)
title (['Περιοδόγραμμα της XIunmod(t)'])
xlabel('frequency (Hz)');
ylabel('PXImodF(F)');

XQunmodF  = fftshift(fft(XQunmod,Nf)*Ts);
PXQunmodF = power(abs(XQunmodF),2)/Ttotal ;
figure(14)
semilogy(F_axis,PXQunmodF)
title (['Περιοδόγραμμα της XQunmod(t)'])
xlabel('frequency (Hz)');
ylabel('PXQmodF(F)');

% A10

XIFiltered = conv(phi,XIunmod)*Ts ; 
XQFiltered = conv(phi,XQunmod)*Ts ;

tnew = linspace(t_phi(1) + tx(1),t_phi(end) + tx(end),length(XIFiltered)) ; % λόγω συνέλιξης
Ttotal = length(tnew) ;

figure(15);
plot (tnew,XIFiltered) ; 
title('conv(phi,XI_noise)');
xlabel('time (sec)');
ylabel('XIFiltered');

figure(16)
plot (tnew,XQFiltered) ; 
title('conv(phi,XQ_noise)');
xlabel('time (sec)');
ylabel('XQFiltered');

XIFilteredF = fftshift(fft(XIFiltered,Nf)*Ts);
PXIFilteredF = power(abs(XIFilteredF),2)/Ttotal ;
figure(17)
semilogy(F_axis,PXIFilteredF)
title (['Περιοδόγραμμα της XIFiltered(t)'])
xlabel('frequency (Hz)');
ylabel('PXIFiltered(F)');

XQFilteredF  = fftshift(fft(XQFiltered,Nf)*Ts);
PXQFilteredF = power(abs(XQFilteredF),2)/Ttotal ;
figure(18)
semilogy(F_axis,PXQFilteredF)
title (['Περιοδόγραμμα της XQFiltered(t)'])
xlabel('frequency (Hz)');
ylabel('PXQFiltered(F)');

% A11

Y = zeros(N,2);
i=1;

for k = 2*A_phi*(T/Ts) :(T/Ts) :(length(XIFiltered)-1) - 2*A_phi*(T/Ts) % N δειγματα ανα Fs άρα 2Ν δείγματα συνολικά
    Y(i,1) = XIFiltered(k);
    Y(i,2) = XQFiltered(k);
    i=i+1;
end

scatterplot (Y)

% A12

est_XI_symb = detect_4_PAM((Y(:,1)), A); % inphase
est_XQ_symb = detect_4_PAM((Y(:,2)), A); % quadrature
 
% A13

est_X_symb = [est_XI_symb est_XQ_symb] ;
Symbol_Errors = symerr(X,est_X_symb);

% A14

est_XI_bits = PAM_4_to_bits(est_XI_symb,A) ; 
est_XQ_bits = PAM_4_to_bits(est_XQ_symb,A) ;

% A15

est_X_bits = [est_XI_bits est_XQ_bits] ; 
Bit_Errors= biterr (bit_seq',est_X_bits) ;

% B1

K = 1000 ; 
SNRdb = 0:2:16 ; 

SentSymbols = N ; 
SentBits = 4*N ; 

Pesymb_exp = zeros (length(SNRdb),1) ;
Pebit_exp = zeros (length(SNRdb),1) ;
Pesymb_th = zeros (length(SNRdb),1) ;
Pebit_th = zeros (length(SNRdb),1) ;

j=1 ; 
SumOfSymbErr = 0 ; 
SumOfBitErr = 0 ;           

for SNR = 0:2:16    % Experiment
           
    SumOfSymbErr = 0 ;
    SumOfBitErr = 0 ; 
    
    for i = 1:K       
        [Symbol_Errors,Bit_Errors] = Error_Possibility (N,A,T,A_phi,roll_off,over,f0,SNR)   ;      
        SumOfSymbErr = Symbol_Errors + SumOfSymbErr ;
        SumOfBitErr = Bit_Errors + SumOfBitErr ;
    end  
     
    % H διαδικασία επαναλαμβάνεται Κ φορές
    Pesymb_exp(j) = SumOfSymbErr /(SentSymbols*K) ;
    Pebit_exp(j) = SumOfBitErr /(SentBits*K) ;    
    
    j = j + 1 ; 
end

% B2

i = 1 ;  % Theoritical Possibilities
for SNR = 0:2:16  
    varianceW = (10*A)/(2*Ts*power(10,SNR/10)); 
    varianceN  = varianceW * Ts * 0.5 ; 
    sigmaN = sqrt(varianceN) ;    
    Pesymb_th(i) = 3*Q(A/sigmaN) ; % Πιθανότητα σφάλματος συμβόλου 16-QAM
    i = i+1 ; 
end

figure(20)
semilogy (SNRdb,Pesymb_exp) ; 
hold on ; 
semilogy (SNRdb,Pesymb_th) ;
title (['Πιθανότητα σφάλματος συμβόλου'])
legend ('Πειραματική','Θεωρητική')
xlabel('SNR(db)');
ylabel('P(Esymbol)');

% B3

for i = 1:length(Pesymb_th)   
    Pebit_th(i) = Pesymb_th(i) / log2(16) ; % log2(16) = bps 16 - QAM
end

figure(21)
semilogy (SNRdb,Pebit_exp) ; 
hold on ; 
semilogy (SNRdb,Pebit_th) ; 
title (['Πιθανότητα σφάλματος bit'])
legend ('Πειραματική','Θεωρητική')
xlabel('SNR(db)');
ylabel('P(Ebit)');
