function K=getTensor_K(d,pdf,Psi,Upsilon,caseStudy,distr,rvar)
% Hugo Esquivel, 2023.
% -
% d = dimensionality of random domain = number of random variables.

R=3; % do not change... this is the number of ODEs in the illustrative example.
P=length(Psi)-1;

if d~=2
    error('d must be equal to 2.')
end

K=zeros(R,R,P+1,P+1);

C=[2 -1 0; -1 1 0; 0 0 0];
D=[0 0 0; 0 1 -1; 0 -1 1];

switch caseStudy
    case 'uniform'
        aMin=rvar{1}.min;
        aMax=rvar{1}.max;

        bMin=rvar{2}.min;
        bMax=rvar{2}.max;

        a=get_rvar_uniform(1,d,aMin,aMax,pdf,Psi,Upsilon,distr);
        b=get_rvar_uniform(2,d,bMin,bMax,pdf,Psi,Upsilon,distr);

    case {'unsymmetricBeta','symmetricBeta'}
        aMin=rvar{1}.min;
        aMax=rvar{1}.max;

        bMin=rvar{2}.min;
        bMax=rvar{2}.max;

        a=get_rvar_beta(1,d,aMin,aMax,pdf,Psi,Upsilon,distr);
        b=get_rvar_beta(2,d,bMin,bMax,pdf,Psi,Upsilon,distr);

    case 'gamma'
        aMin=rvar{1}.min;
        bMin=rvar{2}.min;

        a=get_rvar_gamma(1,d,aMin,pdf,Psi,Upsilon,distr);
        b=get_rvar_gamma(2,d,bMin,pdf,Psi,Upsilon,distr);

    case 'normal'
        aMu=rvar{1}.mu;
        aSigma=rvar{1}.sigma;

        bMu=rvar{2}.mu;
        bSigma=rvar{2}.sigma;

        a=get_rvar_normal(1,d,aMu,aSigma,pdf,Psi,Upsilon,distr);
        b=get_rvar_normal(2,d,bMu,bSigma,pdf,Psi,Upsilon,distr);

    case 'uniformNormal'
        aMin=rvar{1}.min;
        aMax=rvar{1}.max;

        bMu=rvar{2}.mu;
        bSigma=rvar{2}.sigma;

        a=get_rvar_uniform(1,d,aMin,aMax,pdf,Psi,Upsilon,distr);
        b=get_rvar_normal(2,d,bMu,bSigma,pdf,Psi,Upsilon,distr);
end

for i=0:P
    for j=0:P
        K(:,:,i+1,j+1)=a(i+1,j+1)*C+b(i+1,j+1)*D;
    end
end
end

function rvar=get_rvar_uniform(idx,d,rvarmin,rvarmax,pdf,Psi,Upsilon,distr)
P=length(Psi)-1;

UpsilonMoment1=get_UpsilonMoment(1,idx,d,pdf,Psi,distr);

rvar=zeros(P+1);

for i=0:P
    for j=i:P
        rvar(i+1,j+1)=1/2*(rvarmax+rvarmin)*Upsilon(i+1,j+1)+1/2*(rvarmax-rvarmin)*UpsilonMoment1(i+1,j+1);
    end
end

for i=1:P
    for j=0:i-1
        rvar(i+1,j+1)=rvar(j+1,i+1);
    end
end

for i=0:P
    rvar(i+1,:)=rvar(i+1,:)/Upsilon(i+1,i+1);
end
end

function rvar=get_rvar_beta(idx,d,rvarmin,rvarmax,pdf,Psi,Upsilon,distr)
P=length(Psi)-1;

UpsilonMoment1=get_UpsilonMoment(1,idx,d,pdf,Psi,distr);

rvar=zeros(P+1);

for i=0:P
    for j=i:P
        rvar(i+1,j+1)=1/2*(rvarmax+rvarmin)*Upsilon(i+1,j+1)+1/2*(rvarmax-rvarmin)*UpsilonMoment1(i+1,j+1);
    end
end

for i=1:P
    for j=0:i-1
        rvar(i+1,j+1)=rvar(j+1,i+1);
    end
end

for i=0:P
    rvar(i+1,:)=rvar(i+1,:)/Upsilon(i+1,i+1);
end
end

function rvar=get_rvar_gamma(idx,d,rvarmin,pdf,Psi,Upsilon,distr)
P=length(Psi)-1;

UpsilonMoment1=get_UpsilonMoment(1,idx,d,pdf,Psi,distr);

rvar=zeros(P+1);

for i=0:P
    for j=i:P
        rvar(i+1,j+1)=rvarmin*Upsilon(i+1,j+1)+UpsilonMoment1(i+1,j+1);
    end
end

for i=1:P
    for j=0:i-1
        rvar(i+1,j+1)=rvar(j+1,i+1);
    end
end

for i=0:P
    rvar(i+1,:)=rvar(i+1,:)/Upsilon(i+1,i+1);
end
end

function rvar=get_rvar_normal(idx,d,rvarmu,rvarsigma,pdf,Psi,Upsilon,distr)
P=length(Psi)-1;

UpsilonMoment2=get_UpsilonMoment(2,idx,d,pdf,Psi,distr);

rvar=zeros(P+1);

for i=0:P
    for j=i:P
        rvar(i+1,j+1)=(rvarmu-rvarsigma/sqrt(2))*Upsilon(i+1,j+1)+rvarsigma/sqrt(2)*UpsilonMoment2(i+1,j+1);
    end
end

for i=1:P
    for j=0:i-1
        rvar(i+1,j+1)=rvar(j+1,i+1);
    end
end

for i=0:P
    rvar(i+1,:)=rvar(i+1,:)/Upsilon(i+1,i+1);
end
end

function UpsilonMoment=get_UpsilonMoment(power,idx,d,pdf,Psi,distr)
xi=sym(sprintf('xi%d',idx));

P=length(Psi)-1;

UpsilonMoment=zeros(P+1);

for i=0:P
    for j=i:P
        tmp=xi^power*Psi{i+1}*Psi{j+1}*pdf;

        for ii=1:d
            tmp=int(tmp,sym(sprintf('xi%d',ii)),distr{ii}.support);
        end

        UpsilonMoment(i+1,j+1)=double(tmp);
    end
end

for i=1:P
    for j=0:i-1
        UpsilonMoment(i+1,j+1)=UpsilonMoment(j+1,i+1);
    end
end
end
