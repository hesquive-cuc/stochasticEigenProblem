function I=getTensor_I(d,pdf,Psi,Upsilon,distr)
% Hugo Esquivel, 2023.
% -
% d = dimensionality of random domain = number of random variables.

P=length(Psi)-1;

I=zeros(P+1,P+1,P+1);

for i=0:P
    for j=0:P
        for k=j:P
            tmp=Psi{i+1}*Psi{j+1}*Psi{k+1}*pdf;

            for ii=1:d
                tmp=int(tmp,sym(sprintf('xi%d',ii)),distr{ii}.support);
            end

            I(i+1,j+1,k+1)=double(tmp);
        end
    end

    for j=1:P
        for k=0:j-1
            I(i+1,j+1,k+1)=I(i+1,k+1,j+1);
        end
    end

    I(i+1,:,:)=I(i+1,:,:)/Upsilon(i+1,i+1);
end
end
