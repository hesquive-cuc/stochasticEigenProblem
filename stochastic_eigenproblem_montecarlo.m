% Hugo Esquivel, 2023.
% -
% Illustrative example.

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
eigen=1; % Options: 1, 2 or 3
caseStudy='uniform'; % Options: uniform, unsymmetricBeta, symmetricBeta, gamma, normal or uniformNormal

numSimulations=100;
maxPolynomialOrder=2;
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY 1: Performs the Monte Carlo simulation...
% ----------------------------------------------------------------------------------------------------------------------
d=2; % do not change... this is the number of random variables.
R=3; % do not change... this is the number of ODEs in the illustrative example.

[distr,rvar]=getDistribution(caseStudy);

switch caseStudy
    case 'uniform'
        z1=2*rand([numSimulations,1])-1;
        z2=2*rand([numSimulations,1])-1;

        aMin=rvar{1}.min;
        aMax=rvar{1}.max;

        bMin=rvar{2}.min;
        bMax=rvar{2}.max;

        aMC=1/2*(aMax+aMin)+1/2*(aMax-aMin)*z1;
        bMC=1/2*(bMax+bMin)+1/2*(bMax-bMin)*z2;

    case {'symmetricBeta','unsymmetricBeta'}
        z1=2*betarnd(distr{1}.alpha,distr{1}.beta,[numSimulations,1])-1;
        z2=2*betarnd(distr{2}.alpha,distr{2}.beta,[numSimulations,1])-1;

        aMin=rvar{1}.min;
        aMax=rvar{1}.max;

        bMin=rvar{2}.min;
        bMax=rvar{2}.max;

        aMC=1/2*(aMax+aMin)+1/2*(aMax-aMin)*z1;
        bMC=1/2*(bMax+bMin)+1/2*(bMax-bMin)*z2;

    case 'gamma'
        alpha=5;
        beta=1; % do not change because laguerreL function assumes beta=1...

        z1=gamrnd(distr{1}.alpha,1/distr{1}.beta,[numSimulations,1]);
        z2=gamrnd(distr{2}.alpha,1/distr{2}.beta,[numSimulations,1]);

        aMin=rvar{1}.min;
        bMin=rvar{2}.min;
        
        aMC=aMin+z1;
        bMC=bMin+z2;

    case 'normal'
        z1=normrnd(distr{1}.mu,distr{1}.sigma,[numSimulations,1]);
        z2=normrnd(distr{2}.mu,distr{2}.sigma,[numSimulations,1]);

        aMu=rvar{1}.mu;
        aSigma=rvar{1}.sigma;

        bMu=rvar{2}.mu;
        bSigma=rvar{2}.sigma;

        aMC=aMu+aSigma/sqrt(2)*(z1.^2-1);
        bMC=bMu+bSigma/sqrt(2)*(z2.^2-1);

    case 'uniformNormal'
        z1=2*rand([numSimulations,1])-1;
        z2=normrnd(distr{2}.mu,distr{2}.sigma,[numSimulations,1]);

        aMin=rvar{1}.min;
        aMax=rvar{1}.max;

        bMu=rvar{2}.mu;
        bSigma=rvar{2}.sigma;

        distr{1}.name='uniform';
        distr{1}.support=[-1,1];

        aMC=1/2*(aMax+aMin)+1/2*(aMax-aMin)*z1;
        bMC=bMu+bSigma/sqrt(2)*(z2.^2-1);
end

C=[2 -1 0; -1 1 0; 0 0 0];
D=[0 0 0; 0 1 -1; 0 -1 1];

lambdaMC=zeros(R,numSimulations);
phiMC=zeros(R,R,numSimulations);

for i=1:numSimulations
    KMC=aMC(i)*C+bMC(i)*D;

    [VMC,LMC]=eig(KMC);

    lambdaMC(:,i)=diag(LMC);
    phiMC(:,:,i)=VMC;
end
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY 2: Obtains the trial solution for use in the spectral-chaos simulation...
% ----------------------------------------------------------------------------------------------------------------------
syms xi1 xi2

p=maxPolynomialOrder;

P=factorial(d+p)/(factorial(d)*factorial(p))-1;

data=load(sprintf('info_%s_P%d',caseStudy,P)); 
% obtained from: gen_stochastic_eigenvalue.m

Psi=cell(P+1,1);

for i=0:P
    Psi{i+1}=matlabFunction(data.Psi{i+1},'Vars',[xi1,xi2]);
end

lambda=zeros(P+1,1);
phi=zeros(R,P+1);

lambdaMCstar=permute(lambdaMC,[2,1]);
phiMCstar=permute(phiMC,[3,1,2]);

zerovec=0*z1+0*z2;

for i=0:P
    lambda(i+1)=sum((Psi{i+1}(z1,z2)+zerovec).*lambdaMCstar(:,eigen))/sum((Psi{i+1}(z1,z2)+zerovec).^2);

    for j=1:R
        phi(j,i+1)=sum((Psi{i+1}(z1,z2)+zerovec).*phiMCstar(:,j,eigen))/sum((Psi{i+1}(z1,z2)+zerovec).^2);
    end
end

save(sprintf('trial_eigen%d_%s_P%d',eigen,caseStudy,P),'lambda','phi');
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY 3: Obtains the probability moments for lambda and phi statistically...
% ----------------------------------------------------------------------------------------------------------------------
% Note: The larger numSimulations is, the more accurate the probability moments are.

lambdamean1=sum(lambdaMCstar(:,eigen))/numSimulations;
lambdavariance1=sum((lambdaMCstar(:,eigen)-lambdamean1).^2)/(numSimulations-1);

phimean1=zeros(R,1);
phivariance1=zeros(R,1);

for j=1:R
    phimean1(j)=sum(phiMCstar(:,j,eigen))/numSimulations;
    phivariance1(j)=sum((phiMCstar(:,j,eigen)-phimean1(j)).^2)/(numSimulations-1);
end

fprintf('------ Probability moments computed statistically (N = %d) ------\n',numSimulations)
fprintf('\n')

fprintf('lambda mean (:,1) and lambda variance (:,2):\n')
disp([lambdamean1,lambdavariance1])

fprintf('phi mean (:,1) and phi variance (:,2):\n')
disp([phimean1,phivariance1])
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY 4: Obtains the probability moments for lambda and phi spectrally...
% ----------------------------------------------------------------------------------------------------------------------
% Note: The larger numSimulations and maxPolynomialOrder are, the more accurate the probability moments are.

Upsilon=data.Upsilon;

lambdamean2=lambda(1);
lambdavariance2=sum(diag(Upsilon(2:end,2:end)).*lambda(2:end).^2);

phimean2=phi(:,1);

phivariance2=zeros(R,1);

for j=1:R
    phivariance2(j)=sum(diag(Upsilon(2:end,2:end)).*transpose(phi(j,2:end)).^2);
end

fprintf('----- Probability moments computed spectrally (N = %d and P = %d) -----\n',numSimulations,P)
fprintf('\n')

fprintf('lambda mean (:,1) and lambda variance (:,2):\n')
disp([lambdamean2,lambdavariance2])

fprintf('phi mean (:,1) and phi variance (:,2):\n')
disp([phimean2,phivariance2])
% ----------------------------------------------------------------------------------------------------------------------
