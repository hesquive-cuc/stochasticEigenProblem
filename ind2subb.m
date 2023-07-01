function subs=ind2subb(siz,ind)
ndims=length(siz);
subs=zeros(length(ind),ndims);
ind=ind-1;

for i=1:ndims
  r=rem(ind,siz(i));
  subs(:,i)=r;
  ind=(ind-r)/siz(i);
end

subs=subs+1;
end
