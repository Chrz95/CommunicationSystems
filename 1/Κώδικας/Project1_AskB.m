 clc ; 
clear all ; 
close all ; 

%B1

%(α)

T = 1 ;
over = 10 ;
Ts = T/over ; % Περίοδος δειγματοληψίας
Fs = 1/Ts ; % Συχνότητα δειγματοληψίας
A = 4 ;

[phi,t_phi] = srrc_pulse(T, Ts, A, 0);

%(α)
figure(1)
k=1 ; 
phi_moved=[zeros(1,(1/Ts)*k*T) phi(1:end-(1/Ts)*k*T)];
plot (t_phi,phi_moved)
hold on ;
k=4 ; 
phi_moved=[zeros(1,(1/Ts)*k*T) phi(1:end-(1/Ts)*k*T)];
plot (t_phi, phi_moved)
title('Μετατοπισμένοι παλμοί Φ δεξιά κατά κT');
xlabel('time (sec)');
ylabel('Φ(t),Φ(t - κΤ)');
legend('k = 1','k = 4')
hold off ; 

%(β),(γ)

integrals = zeros(1,length(0:2*A));

disp('The integrals of Φ(t)*Φ(t - κΤ) are : ')

for k = 0:2*A % k = [0 1 2 3 4 5 6 7 8]
    figure(2) ; 
    title (['k = ' , num2str(k-1)])
    xlabel('time(sec)');
    ylabel('Φ(t)*Φ(t - κΤ)'); 
    subplot(3,3,k+1) 
        
    phi_moved=[zeros(1,(1/Ts)*k*T) phi(1:end-(1/Ts)*k*T)];
    phi_product = phi.*phi_moved;  
     
    integrals(k+1) = sum(phi_product)*Ts;
    disp(['k = ' , num2str(k)])
    disp(['Sum(Φ(t)*Φ(t - κΤ)) = ' ,num2str(integrals(k+1)) ])
        
    plot(t_phi,phi_product)
    
end



