function checkDistribution(distr)
% Hugo Esquivel, 2023.
% -

d=length(distr);

for i=1:d
    switch distr{i}.name
        case 'uniform'
            if any(distr{i}.support~=[-1,1])
                error('The code requires that uniform''s support is [-1,1].')
            end

        case 'beta'
            if distr{i}.alpha<=0
                error('By definition, beta''s alpha>0.')
            end

            if distr{i}.beta<=0
                error('By definition, beta''s beta>0.')
            end

            if any(distr{i}.support~=[-1,1])
                error('The code requires that beta''s support is [-1,1].')
            end

        case 'gamma'
            if distr{i}.alpha<=0
                error('By definition, gamma''s alpha>0.')
            end

            if distr{i}.beta<=0
                error('By definition, gamma''s beta>0.')
            end

            if distr{i}.beta~=1
                error('The code requires that gamma''s beta = 1.') % this is because laguerreL function assumes beta = 1...
            end

            if any(distr{i}.support~=[0,Inf])
                error('The code requires that gamma''s support is [0,Inf].')
            end

        case 'normal'
            if distr{i}.sigma<=0
                error('By definition, normal''s sigma>0.')
            end

            if distr{i}.mu~=0
                error('The code requires that normal''s mu = 0.')
            end

            if distr{i}.sigma~=1
                error('The code requires that normal''s sigma = 1.')
            end

            if any(distr{i}.support~=[-Inf,Inf])
                error('The code requires that normal''s support is [-Inf,Inf].')
            end
    end
end
end
