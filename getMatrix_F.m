function F=getMatrix_F(K,I,eta,x)
% Hugo Esquivel, 2023.
% -

R=size(K,1);
P=size(K,3)-1;
Q=(R+1)*(P+1);

delta=eye(R);

[lambda,phi]=get_lambda_phi(x,R,P);

F=zeros(Q);

% F-block 11: F^u^i_beta
for u=1:R
    for i=0:P
        for beta=0:P
            for j=0:P
                F(u+R*i,beta+1)=F(u+R*i,beta+1)-I(i+1,beta+1,j+1)*phi(u,j+1);
            end
        end
    end
end

% F-block 12: F^u_alpha^i_beta
for u=1:R
    for alpha=1:R
        for i=0:P
            for beta=0:P
                F(u+R*i,alpha+R*beta+P+1)=K(u,alpha,i+1,beta+1);

                for j=0:P
                    F(u+R*i,alpha+R*beta+P+1)=F(u+R*i,alpha+R*beta+P+1)-I(i+1,beta+1,j+1)*lambda(j+1)*delta(u,alpha);
                end
            end
        end
    end
end

% F-block 21: F^i_beta
% for i=0:P
%     for beta=0:P
%         F(i+R*(P+1)+1,beta+1)=0;
%     end
% end

% F-block 22: F_alpha^i_beta
for alpha=1:R
    for i=0:P
        for beta=0:P
            for u=1:R
                F(i+R*(P+1)+1,alpha+R*beta+P+1)=F(i+R*(P+1)+1,alpha+R*beta+P+1)-2*eta(alpha,u)*F(u+R*i,beta+1);
            end
        end
    end
end
end
