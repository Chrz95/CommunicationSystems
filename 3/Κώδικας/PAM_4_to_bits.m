function [est_X] = PAM_4_to_bits(X,A)

est_X = ones(1,2*length(X));

for i=1:length(X)
    if(X(i)==-3*A)
        est_X(2*i-1) = 0;
        est_X(2*i) =0;
    elseif(X(i)==-A)
        est_X(2*i-1) = 0;
        est_X(2*i)= 1;
    elseif(X(i)==A)
        est_X(2*i-1) = 1;
        est_X(2*i)= 1;
    elseif(X(i)==3*A)
        est_X(2*i-1) = 1;
        est_X(2*i)= 0;
    end
end
end
