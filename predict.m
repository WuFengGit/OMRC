function pre = predict(cls, fea)
fea(:,:,end+1) = 1;
[M,~] = size(cls);
N = size(fea,2);
pre = zeros(N,M);

for i=1:N
    pre(i,:) = sum(cls.*permute(fea(:,i,:),[1,3,2]),2);
end
pre = pre';
pre = exp(pre)./(exp(pre)+exp(-pre))-0.5;
end