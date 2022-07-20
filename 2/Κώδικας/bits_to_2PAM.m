function [PAM_2]=bits_to_2PAM(bits) 

PAM_2 = zeros(1,length(bits)) ;

for i=1:length(bits)
  if bits(i) == 1;
    PAM_2(i) = -1;
  elseif bits(i) == 0 
    PAM_2(i) = 1 ;
  else
    PAM_2(i) = 0;
  end
end