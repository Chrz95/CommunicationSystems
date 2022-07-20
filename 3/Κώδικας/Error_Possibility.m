function [Symbol_Errors,Bit_Errors] = Error_Possibility(N,A,T,A_phi,roll_off,over,f0,SNRdb)  

bit_seq = (sign(randn(4*N,1))+1)/2;

X = bits_to_4_PAM(bit_seq, A);

XI = X(1:N);
XQ = X(N+1:end);

Ts = T/over ; % Περίοδος δειγματοληψίας
Fs = 1/Ts ; % Συχνότητα δειγματοληψίας

[phi,t_phi] = srrc_pulse(T,Ts,A_phi,roll_off);

[tx,XIt] = PAM4_function (XI,Ts,over,phi,t_phi) ;
[tx,XQt] = PAM4_function (XQ,Ts,over,phi,t_phi) ;

XImod = XIt.*(2*cos(2*pi*f0*tx));
XQmod = XQt.*((-2)*sin(2*pi*f0*tx));

Xmod = XImod + XQmod;

mu = 0 ;
varianceW = (10*A)/(2*Ts*power(10,SNRdb/10));
sigmaW = sqrt(varianceW) ; 
Noise = normrnd(mu,sigmaW,1,length(tx)) ;

XmodNoise = Xmod + Noise ; 

XIunmod = XmodNoise .* cos(2*pi*f0*tx) ;
XQunmod = -XmodNoise .* sin(2*pi*f0*tx) ; 

XIFiltered = conv(phi,XIunmod)*Ts ; 
XQFiltered = conv(phi,XQunmod)*Ts ;

Y = zeros(N,2);
i=1;

for k = 2*A_phi*(T/Ts) :(T/Ts):(length(XIFiltered)-1) - 2*A_phi*(T/Ts) % N δειγματα ανα Fs άρα 2Ν δείγματα
    Y(i,1)=XIFiltered(k);
    Y(i,2)=XQFiltered(k);
    i=i+1;
end

est_XI_symb = detect_4_PAM((Y(:,1)), A); % inphase
est_XQ_symb = detect_4_PAM((Y(:,2)), A); % quadrature
est_X_symb = [est_XI_symb est_XQ_symb] ; 

Symbol_Errors = symerr(X,est_X_symb) ;

est_XI_bits = PAM_4_to_bits(est_XI_symb,A) ; 
est_XQ_bits = PAM_4_to_bits(est_XQ_symb,A) ;
est_X_bits = [est_XI_bits est_XQ_bits] ; 

Bit_Errors = biterr (bit_seq',est_X_bits) ;