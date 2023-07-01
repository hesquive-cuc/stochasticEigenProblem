function Psi=get_Psi(d,p,distr)
% Hugo Esquivel, 2023.
% -
% d = dimensionality of random domain = number of random variables.
% p = maximum polynomial order.

Psi1=cell(p+1,d);

for i=1:d
    switch distr{i}.name
        case 'uniform'
            for j=0:p
                xi=sym(sprintf('xi%d',i));

                Psi1{j+1,i}(xi)=expand(legendreP(j,xi));
            end

        case 'beta'
            alpha=distr{i}.alpha;
            beta=distr{i}.beta;

            for j=0:p
                xi=sym(sprintf('xi%d',i));

                Psi1{j+1,i}(xi)=expand(jacobiP(j,beta-1,alpha-1,xi));
            end

        case 'gamma'
            alpha=distr{i}.alpha;

            for j=0:p
                xi=sym(sprintf('xi%d',i));

                Psi1{j+1,i}(xi)=expand(laguerreL(j,alpha-1,xi));
            end

        case 'normal'
            for j=0:p
                xi=sym(sprintf('xi%d',i));

                Psi1{j+1,i}(xi)=expand(hermiteH(j,xi/sqrt(sym(2)))/2^(j/2));
            end
    end
end

midx=ind2subb((p+1)*ones(d,1),1:(p+1)^d); % obtained from performing a full tensor product among indices...

P=factorial(d+p)/(factorial(d)*factorial(p))-1;

Psi=cell(P+1,1); % obtained from performing a total-order tensor product for simplicity...

k=0;

for order=0:p
    for j=1:size(midx,1)
        if sum(midx(j,:))-d==order
            k=k+1;

            Psi{k}=1;

            for i=1:d
                xi=sym(sprintf('xi%d',i));

                Psi{k}=Psi{k}*Psi1{midx(j,i),i}(xi);
            end

            Psi{k}=expand(Psi{k});
        end
    end
end
end
