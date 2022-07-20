function [X] = bits_to_4_PAM(bit_seq,A )

X = zeros(1,length(bit_seq)/2);
j=1;

for i=1:2:length(bit_seq)
    if(bit_seq(i)==0 && bit_seq(i+1)==0)
        X(j)=-3*A;
    elseif(bit_seq(i)==0 && bit_seq(i+1)==1)
        X(j)=(-1)*A;
    elseif(bit_seq(i)==1 && bit_seq(i+1)==1)
        X(j)=A;
    elseif(bit_seq(i)==1 && bit_seq(i+1)==0)
        X(j)=3*A;
    else
         X(j)=0*A;
    end
   j=j+1;
end