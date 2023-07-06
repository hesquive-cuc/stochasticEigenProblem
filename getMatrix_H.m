function H=getMatrix_H(I,K,eta,x,X,f)
% Hugo Esquivel, 2023.
% -

xbar=x-1/2*X*f;

H=getMatrix_F(I,K,eta,xbar);
end
