function Upsilon=get_Upsilon(d,pdf,Psi,distr)
% Hugo Esquivel, 2023.
% -
% d = dimensionality of random domain = number of random variables.

P=length(Psi)-1;

Upsilon=zeros(P+1);

% for i=0:P
%     for j=i:P
%         tmp=Psi{i+1}*Psi{j+1}*pdf;
% 
%         for ii=1:d
%             tmp=int(tmp,sym(sprintf('xi%d',ii)),distr{ii}.support);
%         end
% 
%         Upsilon(i+1,j+1)=double(tmp);
%     end
% end
% 
% for i=1:P
%     for j=0:i-1
%         Upsilon(i+1,j+1)=Upsilon(j+1,i+1);
%     end
% end

% Or taking advantage of orthogonality of Psi functions...
for i=0:P
    tmp=Psi{i+1}*Psi{i+1}*pdf;

    for ii=1:d
        tmp=int(tmp,sym(sprintf('xi%d',ii)),distr{ii}.support);
    end

    Upsilon(i+1,i+1)=double(tmp);
end
end
