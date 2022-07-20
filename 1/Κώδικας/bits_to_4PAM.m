function [PAM_4]=bits_to_4PAM(bits) 

PAM_4 = zeros(1,length(bits(:,1))) ;

for i=1:length(bits(:,1))
  if bits(i,1) == 0 && bits(i,2) == 0 ;
    PAM_4(i) = -3;
  elseif bits(i,1) == 0 && bits(i,2) == 1 ;
    PAM_4(i) = 1 ;
  elseif bits(i,1) == 1 && bits(i,2) == 1;
    PAM_4(i) = -1 ;
  elseif bits(i,1) == 1 && bits(i,2) == 0 ;
    PAM_4(i) = 3 ;
  else
    PAM_4(i) = 0;
  end
end