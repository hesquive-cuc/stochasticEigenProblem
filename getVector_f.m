function f=getVector_f(K,I,eta,x)
% Hugo Esquivel, 2023.
% -

R=size(K,1);
P=size(K,3)-1;
Q=(R+1)*(P+1);

[lambda,phi]=get_lambda_phi(x,R,P);

f=zeros(Q,1);

% f-block 1: g^u^i
for u=1:R
    for i=0:P
        for j=0:P
            for v=1:R
                f(u+R*i)=f(u+R*i)+K(u,v,i+1,j+1)*phi(v,j+1);
            end
            
            for k=0:P
                f(u+R*i)=f(u+R*i)-I(i+1,j+1,k+1)*lambda(k+1)*phi(u,j+1);
            end
        end
    end
end

% f-block 2: h^i
for i=0:P
    for j=0:P
        for k=0:P
            for u=1:R
                for v=1:R
                    f(i+R*(P+1)+1)=f(i+R*(P+1)+1)+I(i+1,j+1,k+1)*eta(u,v)*phi(u,j+1)*phi(v,k+1);
                end
            end
        end
    end

    f(i+R*(P+1)+1)=f(i+R*(P+1)+1)-I(i+1,1,1);
end
end
