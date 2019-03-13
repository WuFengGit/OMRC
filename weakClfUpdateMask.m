function cls = weakClfUpdateMask(cls, posx, negx, mask)
numpos = size(posx.feature_hog,2);
numneg = size(negx.feature_hog,2);
M_hog = size(posx.feature_hog,1);
dim = size(posx.feature_hog,3);
if isempty(cls)
    cls = zeros(M_hog, dim+1);
end
for i=1:M_hog 
    if mask(i)==0
        continue;
    end
    %% update with current frame
    x = [permute(posx.feature_hog(i,:,:), [2,3,1]); permute(negx.feature_hog(i,:,:), [2,3,1])];
    y = [zeros(numpos,1)+1;zeros(numneg,1)-1];
    w = [posx.w'/sum(posx.w); zeros(numneg,1)+1/numneg];
    [cls(i,1:end-1)] = BPA(cls(i,1:end-1)', x, y, w*0.5);
end
end