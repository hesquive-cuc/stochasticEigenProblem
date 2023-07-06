function x=getVector_x(lambda,phi)
% Hugo Esquivel, 2023.
% -

R=size(phi,1);
P=size(phi,2)-1;
Q=(R+1)*(P+1);

x=zeros(Q,1);

% x-block 1: lambda^i
for i=0:P
    x(i+1)=lambda(i+1);
end

% x-block 2: phi^u^i
for u=1:R
    for i=0:P
        x(u+R*i+P+1)=phi(u,i+1);
    end
end
end
