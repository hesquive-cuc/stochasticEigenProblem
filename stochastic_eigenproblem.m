% Hugo Esquivel, 2023.
% -
% Illustrative example.

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
eigen=1; % Options: 1, 2 or 3
caseStudy='uniform'; % Options: uniform, unsymmetricBeta, symmetricBeta, gamma, normal or uniformNormal

numRandomVariables=2;
maxPolynomialOrder=2;

runMethod='Halley'; % Options: Newton or Halley

maxNumIterations=100;
tol=1e-12;
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY 1: Performs the spectral-chaos simulation...
% ----------------------------------------------------------------------------------------------------------------------
p=maxPolynomialOrder;
d=numRandomVariables;

P=factorial(d+p)/(factorial(d)*factorial(p))-1;

data=load(sprintf('info_%s_P%d',caseStudy,P));
% obtained from: gen_stochastic_eigenproblem.m

Psi=data.Psi;
Upsilon=data.Upsilon;
I=data.I;
K=data.K;
eta=data.eta;

R=size(K,1);
P=size(K,3)-1;

data=load(sprintf('trial_eigen%d_%s_P%d.mat',eigen,caseStudy,P));
% obtained from: stochastic_eigenproblem_montercarlo.m

lambda=data.lambda;
phi=data.phi;

x=getVector_x(lambda,phi);
f=getVector_f(I,K,eta,x);

df=norm(f);

fprintf('Iteration %d: %e...\n',0,df)

greenFlag=false;

tic

switch runMethod
    case 'Newton' % Old approach: Ghanem and Gosh's approach (2007 paper)...
        % "Efficient characterization of the random eigenvalue problem in a polynomial chaos decomposition."
        for ii=1:maxNumIterations
            F=getMatrix_F(I,K,eta,x);

            x=x-F\f;

            f=getVector_f(I,K,eta,x);

            df=norm(f);

            fprintf('Iteration %d: %e...\n',ii,df)

            if df<tol
                greenFlag=true;

                break
            end
        end

    case 'Halley' % New approach...
        for ii=1:maxNumIterations
            F=getMatrix_F(I,K,eta,x);
            X=getMatrix_X(F);
            H=getMatrix_H(I,K,eta,x,X,f);

            x=x-H\f;

            f=getVector_f(I,K,eta,x);

            df=norm(f);

            fprintf('Iteration %d: %e...\n',ii,df)

            if df<tol
                greenFlag=true;

                break
            end
        end
end

fprintf('\n')

toc

fprintf('\n')

if greenFlag==false
    fprintf('>>>>> Simulation DID NOT converge <<<<<\n')
    fprintf('\n')
end

[lambda,phi]=get_lambda_phi(x,R,P);
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY 2: Obtains the probability moments for lambda and phi spectrally...
% ----------------------------------------------------------------------------------------------------------------------
lambdamean=lambda(1);
lambdavariance=sum(diag(Upsilon(2:end,2:end)).*lambda(2:end).^2);

phimean=phi(:,1);

phivariance=zeros(R,1);

for j=1:R
    phivariance(j)=sum(diag(Upsilon(2:end,2:end)).*transpose(phi(j,2:end)).^2);
end

fprintf('lambda mean (:,1) and lambda variance (:,2):\n')
disp([lambdamean,lambdavariance])

fprintf('phi mean (:,1) and phi variance (:,2):\n')
disp([phimean,phivariance])
% ----------------------------------------------------------------------------------------------------------------------
