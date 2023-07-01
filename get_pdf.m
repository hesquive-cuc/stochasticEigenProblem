function pdf=get_pdf(d,distr)
% Hugo Esquivel, 2023.
% -
% d = dimensionality of random domain = number of random variables.

pdf=1;

for i=1:d
    switch distr{i}.name
        case 'uniform' % on [-1,1]
            pdf=pdf*1/2;

        case 'beta' % on [-1,1]
            alpha=distr{i}.alpha;
            beta=distr{i}.beta;

            xi=sym(sprintf('xi%d',i));

            pdf=pdf*(1+xi)^(alpha-1)*(1-xi)^(beta-1)/(2^(alpha+beta-1)*betafun(alpha,beta));

        case 'gamma' % on [0,Inf]
            alpha=distr{i}.alpha;
            beta=distr{i}.beta;

            xi=sym(sprintf('xi%d',i));

            pdf=pdf*beta^alpha/gamma(alpha)*xi^(alpha-1)*exp(-beta*xi);

        case 'normal' % on [-Inf,Inf]
            mu=distr{i}.mu;
            sigma=distr{i}.sigma;

            xi=sym(sprintf('xi%d',i));

            pdf=pdf*1/(sigma*sqrt(2*sym(pi)))*exp(-1/2*((xi-mu)/sigma)^2);
    end
end
end

function betaval=betafun(z1,z2)
betaval=beta(z1,z2);
end
