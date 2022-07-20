function [est_X] = detect_4_PAM(Y,A )

est_X = zeros (1,length(Y)) ;

for i=1:length(Y)
   if (Y(i)>= 2*A) % D3
       est_X(i) = +3*A;       
   elseif ((Y(i)< 2*A) && (Y(i) >= 0)) % D2
       est_X(i) = +A  ;              
   elseif ((Y(i)> -2*A) && (Y(i) < 0)) % D1
       est_X(i) = -A  ;   
   elseif (Y(i)<= -2*A) % D0
       est_X(i) = -3*A ; 
   end
end

