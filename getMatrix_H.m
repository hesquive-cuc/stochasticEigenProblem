function H=getMatrix_H(I,K,eta,lambda,phi,X,f)
% Hugo Esquivel, 2023.
% -

R=size(K,1);
P=size(K,3)-1;

delta=eye(R);

H=zeros((R+1)*(P+1));

lambdabar=zeros(P+1,1);
phibar=zeros(R,P+1);

% lambdabar:
for j=0:P
    lambdabar(j+1)=lambda(j+1);

    for rho=0:P
        for gamma=1:R
            lambdabar(j+1)=lambdabar(j+1)-1/2*X(j+1,gamma+R*rho)*f(gamma+R*rho);
        end
    end

    for rho=0:P
        lambdabar(j+1)=lambdabar(j+1)-1/2*X(j+1,rho+R*(1+P)+1)*f(R*(1+P)+rho+1);
    end
end

% phibar:
for u=1:R
    for j=0:P
        phibar(u,j+1)=phi(u,j+1);

        for rho=0:P
            for gamma=1:R
                phibar(u,j+1)=phibar(u,j+1)-1/2*X(u+R*j+P+1,gamma+R*rho)*f(gamma+R*rho);
            end
        end

        for rho=0:P
            phibar(u,j+1)=phibar(u,j+1)-1/2*X(u+R*j+P+1,rho+R*(1+P)+1)*f(R*(1+P)+rho+1);
        end
    end
end

% H-block 11: H^u^i_beta
for u=1:R
    for i=0:P
        for beta=0:P
            for j=0:P
                H(u+R*i,beta+1)=H(u+R*i,beta+1)-I(i+1,beta+1,j+1)*phibar(u,j+1);
            end
        end
    end
end

% H-block 12: H^u_alpha^i_beta
for u=1:R
    for alpha=1:R
        for i=0:P
            for beta=0:P
                H(u+R*i,alpha+R*beta+P+1)=K(u,alpha,i+1,beta+1);

                for j=0:P
                    H(u+R*i,alpha+R*beta+P+1)=H(u+R*i,alpha+R*beta+P+1)-I(i+1,beta+1,j+1)*lambdabar(j+1)*delta(u,alpha);
                end
            end
        end
    end
end

% H-block 21: H^i_beta
for i=0:P
    for beta=0:P
        H(i+R*(1+P)+1,beta+1)=0;
    end
end

% H-block 22: H_alpha^i_beta
for alpha=1:R
    for i=0:P
        for beta=0:P
            for u=1:R
                H(i+R*(1+P)+1,alpha+R*beta+P+1)=H(i+R*(1+P)+1,alpha+R*beta+P+1)-2*eta(alpha,u)*H(u+R*i,beta+1);
            end
        end
    end
end
end
