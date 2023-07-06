function [lambda,phi]=get_lambda_phi(x,R,P)
% Hugo Esquivel, 2023.
% -

lambda=zeros(P+1,1);

for i=0:P
    lambda(i+1)=x(i+1);
end

phi=zeros(R,P+1);

for u=1:R
    for i=0:P
        phi(u,i+1)=x(u+R*i+P+1);
    end
end
end
