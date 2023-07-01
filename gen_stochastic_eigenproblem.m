% Hugo Esquivel, 2023.
% -
% Illustrative example.

clearvars; close all; clc;

% ----------------------------------------------------------------------------------------------------------------------
% INPUT:
% ----------------------------------------------------------------------------------------------------------------------
caseStudy='uniform'; % Options: uniform, unsymmetricBeta, symmetricBeta, gamma, normal or uniformNormal

numRandomVariables=2;
maxPolynomialOrder=2;
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% BODY:
% ----------------------------------------------------------------------------------------------------------------------
d=numRandomVariables;
p=maxPolynomialOrder;

[distr,rvar]=getDistribution(caseStudy);

tic

pdf=get_pdf(d,distr);
Psi=get_Psi(d,p,distr);
Upsilon=get_Upsilon(d,pdf,Psi,distr);
I=getTensor_I(d,pdf,Psi,Upsilon,distr);
K=getTensor_K(d,pdf,Psi,Upsilon,caseStudy,distr,rvar);
eta=getTensor_eta();

toc

P=factorial(d+p)/(factorial(d)*factorial(p))-1;

save(sprintf('info_%s_P%d',caseStudy,P),'pdf','Psi','Upsilon','I','K','eta');
% ----------------------------------------------------------------------------------------------------------------------
