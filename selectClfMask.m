function [clspool_hog, clspool_hist, posx, negx, selector] = selectClfMask(clspool_hog, clspool_hist, posx, negx, mask1, mask2)
mask1(1) = 1;mask2(1) = 1;
%--------Update the weak classifiers with mask1==1
clspool_hog = weakClfUpdateMask(clspool_hog, posx, negx, mask1);
tmpposx.feature_hog = posx.feature_hist;tmpposx.w = posx.w;
tmpnegx.feature_hog = negx.feature_hist;
clspool_hist = weakClfUpdateMask(clspool_hist, tmpposx, tmpnegx, mask1);


posx.pospred = [predict(clspool_hog, posx.feature_hog);predict(clspool_hist, posx.feature_hist)];
negx.negpred = [predict(clspool_hog, negx.feature_hog);predict(clspool_hist, negx.feature_hist)];


%--------select features with mask2==1
mask2 = [mask2;mask2];

idx = find(mask2>0);
tmp_posx.w = posx.w;
tmp_posx.pospred = posx.pospred(idx,:);
tmp_negx.negpred = negx.negpred(idx,:);
if ~isempty(posx.first)
    tmp_posx.first.pospred = posx.first.pospred(idx,:);
    tmp_posx.first.w = posx.first.w;
else
    tmp_posx.first = [];
end
if sum(mask2)<=15
    selector = zeros(length(mask1)+1,1);
    selector(idx+1) = 1;
else
    tmp_selector = sel(tmp_posx, tmp_negx);

    selector = zeros(length(mask1)+1,1);
    selector(1) = tmp_selector(1);
    selector(idx+1) = tmp_selector(2:end);
end

end

function selector =  sel(posx, negx)
tmp_posx = posx;
if ~isempty(posx.first)
    tmp_posx.pospred = [posx.pospred posx.first.pospred];
    tmp_posx.w = [posx.w,posx.first.w];
end
tmp_posx.w = tmp_posx.w/sum(tmp_posx.w);

selectorI = clfMilBoostUpdate(tmp_posx,negx,15);
selector = zeros(size(posx.pospred,1)+1, 1);
selector(selectorI+1) = 1;
end


function [selector] = clfMilBoostUpdate(posx,negx,NumSel)
[~,numpos] = size(posx.pospred);
[M,numneg] = size(negx.negpred);

selector = zeros(1,0);

hmp = posx.pospred;
hmn = negx.negpred;
Hijp = zeros(1,numpos);
Hijn = zeros(1,numneg);
Lm = zeros(M,1);

selector(end+1) = 1;
selector(end+1) = M/2+1;

for k=3:NumSel
    for m=1:M
        pijm_p = 1./(1+exp(-Hijp-hmp(m,:)));
        pijm_n = 1./(1+exp(-Hijn-hmn(m,:)));
        Lm(m) = log(sum(pijm_p)/numpos) + log(sum(1-pijm_n)/numneg);
    end  
    [~, Indexs] = sort(-Lm);
    for m=1:length(Indexs)
        if isempty(find(selector==Indexs(m), 1))
            selector(end+1) = Indexs(m);
            Hijp = Hijp + hmp(Indexs(m),:);
            Hijn = Hijn + hmn(Indexs(m),:);
            break;
        end
    end

end
end


function [selector] = clfMilBoostUpdateW(posx,negx,NumSel)
[~,numpos] = size(posx.pospred);
[M,numneg] = size(negx.negpred);

selector = zeros(1,0);

hmp = posx.pospred;
hmn = negx.negpred;
Hijp = zeros(1,numpos);
Hijn = zeros(1,numneg);
Lm = zeros(M,1);

selector(end+1) = 1;
selector(end+1) = M/2+1;

for k=3:NumSel
    for m=1:M
        pijm_p = 1./(1+exp(-Hijp-hmp(m,:)));
        pijm_n = 1./(1+exp(-Hijn-hmn(m,:)));
        Lm(m) = log(pijm_p*posx.w') + log(sum(1-pijm_n)/numneg);
    end  
    [~, Indexs] = sort(-Lm);
    for m=1:length(Indexs)
        if isempty(find(selector==Indexs(m), 1))
            selector(end+1) = Indexs(m);
            Hijp = Hijp + hmp(Indexs(m),:);
            Hijn = Hijn + hmn(Indexs(m),:);
            break;
        end
    end

end
end