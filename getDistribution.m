function [distr,rvar]=getDistribution(caseStudy)
% Hugo Esquivel, 2023.
% -
% d = dimensionality of random domain = number of random variables.

d=2; % do not change...

distr=cell(d,1);
rvar=cell(d,1);

switch caseStudy
    case 'uniform'
        for i=1:d
            distr{i}.name='uniform';
            distr{i}.support=[-1,1];
        end

        rvar{1}.min=15;
        rvar{1}.max=25;

        rvar{2}.min=5;
        rvar{2}.max=25;

    case 'unsymmetricBeta'
        for i=1:d
            distr{i}.name='beta';
            distr{i}.alpha=2;
            distr{i}.beta=5;
            distr{i}.support=[-1,1];
        end

        rvar{1}.min=15;
        rvar{1}.max=25;

        rvar{2}.min=5;
        rvar{2}.max=25;

    case 'symmetricBeta'
        for i=1:d
            distr{i}.name='beta';
            distr{i}.alpha=2;
            distr{i}.beta=2;
            distr{i}.support=[-1,1];
        end

        rvar{1}.min=15;
        rvar{1}.max=25;

        rvar{2}.min=5;
        rvar{2}.max=25;

    case 'gamma'
        for i=1:d
            distr{i}.name='gamma';
            distr{i}.alpha=5;
            distr{i}.beta=1;
            distr{i}.support=[0,Inf];
        end

        rvar{1}.min=15;
        rvar{2}.min=5;

    case 'normal'
        for i=1:d
            distr{i}.name='normal';
            distr{i}.mu=0;
            distr{i}.sigma=1;
            distr{i}.support=[-Inf,Inf];
        end

        rvar{1}.mu=15;
        rvar{1}.sigma=3;

        rvar{2}.mu=20;
        rvar{2}.sigma=5;

    case 'uniformNormal'
        distr{1}.name='uniform';
        distr{1}.support=[-1,1];

        rvar{1}.min=15;
        rvar{1}.max=25;

        distr{2}.name='normal';
        distr{2}.mu=0;
        distr{2}.sigma=1;
        distr{2}.support=[-Inf,Inf];

        rvar{2}.mu=20;
        rvar{2}.sigma=5;
end

checkDistribution(distr)
end
