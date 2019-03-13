benchmark_dir = './OTB2015/';
seq = 'David';

SHOW_PROCEDURE = true;
addpath('./tool');
addpath('./auxiliary');
addpath([ benchmark_dir seq '/img/' ]);

tic
rand('state',sum(100*clock));
img_dir = dir([benchmark_dir seq '/img/*.jpg']);
gt = load([ benchmark_dir seq '/groundtruth_rect.txt']);
initstate = gt(1,:);%initial tracker
num = length(img_dir);% number of frames
%% Parameter Settings
trparams.init_negnumtrain = 50;%number of trained negative samples
trparams.init_postrainrad = 4.0;%radical scope of positive samples; boy 8
trparams.init_postrainrad2 = 20;
trparams.initstate = initstate;% object position [x y width height]
trparams.srchwinsz = 25;% size of search window

% classifier parameters
clfparams.width = trparams.initstate(3);
clfparams.height= trparams.initstate(4);
M_hog = 150;
clspool_hog = [];clspool_hist = [];

%compute feature template
[ftrPos_hog, M_hog] = hogFtr(clfparams,M_hog);

%% initilize the first frame
img0 = imread(img_dir(1).name);
if size(img0,3)>1
    img0 = double(rgb2gray(img0));
else
    img0 = double(img0);
end
gaussianFilter = fspecial('gaussian',[5 5],1);
img = imfilter(img0,gaussianFilter,'conv');

% pixel model
bgmodel = img0;
bgmodel(initstate(2):initstate(2)+initstate(4)-1,initstate(1):initstate(1)+initstate(3)-1) = -1;
tarmodel = img0(initstate(2):initstate(2)+initstate(4)-1,initstate(1):initstate(1)+initstate(3)-1);
proimg = ones(initstate(4), initstate(3));

% extract hog & hist features
posx.sampleImage = sampleImg(img,initstate,trparams.init_postrainrad,0,100000);
negx.sampleImage = sampleImg(img,initstate,2*trparams.srchwinsz,1.5*trparams.init_postrainrad,trparams.init_negnumtrain);
samplePos.sx = [posx.sampleImage.sx,negx.sampleImage.sx];
samplePos.sy = [posx.sampleImage.sy,negx.sampleImage.sy];
feature_hog = reshape(HogExtract36(img, [samplePos.sx',samplePos.sy'],ftrPos_hog),...
    size(ftrPos_hog,1), length(samplePos.sx),9*4);
feature_hist = reshape(HistExtract80(img0, [samplePos.sx',samplePos.sy'],ftrPos_hog),...
    size(ftrPos_hog,1), length(samplePos.sx),80);
posx.feature_hog = feature_hog(:,1:length(posx.sampleImage.sx),:);
negx.feature_hog = feature_hog(:,length(posx.sampleImage.sx)+1:end,:);
posx.feature_hist = feature_hist(:,1:length(posx.sampleImage.sx),:);
negx.feature_hist = feature_hist(:,length(posx.sampleImage.sx)+1:end,:);
posx.first = {};

% train & select weak classifiers
posx.w = 1./(1+((posx.sampleImage.sx-initstate(1)).^2+(posx.sampleImage.sy-initstate(2)).^2));
[clspool_hog, clspool_hist, posx, negx, selector] = selectClfMask(clspool_hog, clspool_hist, posx, negx, ones(M_hog,1), ones(M_hog,1));

%% Start tracking
trackpos = zeros(num,4);
trackpos(1,:) = initstate;
dist = zeros(num,1);

fprintf('%d ',num); 
th = 0;
for i = 2:num    
    if mod(i,20)==0
        fprintf('.');
    end
    
    img0 = imread(img_dir(i).name); 
    if size(img0,3)>1
        img0 = double(rgb2gray(img0));
    else
        img0 = double(img0);
    end
    img = imfilter(img0,gaussianFilter,'conv');
    
    % detection
    detectx.sampleImage = sampleImg(img,initstate,trparams.srchwinsz,0,10000);
    selector_hog =  selector(2:M_hog+1);selector_hist =  selector(M_hog+2:end);
    prob = zeros(1,length(detectx.sampleImage));
    if sum(selector_hog)>0
        detectx.feature_hog = reshape(HogExtract36(img, [detectx.sampleImage.sx',detectx.sampleImage.sy'],ftrPos_hog(selector_hog>0,:)),...
            sum(selector_hog>0), length(detectx.sampleImage.sx),9*4);
        prob = selector_hog(selector_hog>0)'*predict(clspool_hog(selector_hog>0,:), detectx.feature_hog);
    end
    if sum(selector_hist)>0
        detectx.feature_hist = reshape(HistExtract80(img0, [detectx.sampleImage.sx',detectx.sampleImage.sy'],ftrPos_hog(selector_hist>0,:)),...
        sum(selector_hist>0), length(detectx.sampleImage.sx),80);
        prob = prob + selector_hist(selector_hist>0)'*predict(clspool_hist(selector_hist>0,:), detectx.feature_hist);
    end    
    if SHOW_PROCEDURE
        showDetection(i, img, prob, detectx.sampleImage, size(img), ...
            find(selector_hog>0),find(selector_hist>0), ftrPos_hog);
    end
    
    [c,index] = max(prob);
    initstate = [detectx.sampleImage.sx(index) detectx.sampleImage.sy(index) detectx.sampleImage.sw(index) detectx.sampleImage.sh(index)];

    trackpos(i,:) = initstate;
    dist(i) = sqrt((gt(i,1)+gt(i,3)/2-initstate(1)-initstate(3)/2)^2 + ...
        (gt(i,2)+gt(i,4)/2-initstate(2)-initstate(4)/2)^2);

    % sampling and update classifiers 
    posx.sampleImage = sampleImg(img,initstate,trparams.init_postrainrad,0,10000);
    negx.sampleImage = sampleImg(img,initstate,1.5*trparams.srchwinsz,4+trparams.init_postrainrad,trparams.init_negnumtrain);
    % weight of the positive instance     
    posx.w = 1./(1+((posx.sampleImage.sx-initstate(1)).^2+(posx.sampleImage.sy-initstate(2)).^2));
    samplePos.sx = [posx.sampleImage.sx,negx.sampleImage.sx];
    samplePos. sy = [posx.sampleImage.sy,negx.sampleImage.sy];
    feature_hog = reshape(HogExtract36(img, [samplePos.sx',samplePos.sy'],ftrPos_hog),...
        size(ftrPos_hog,1), length(samplePos.sx),9*4);
    feature_hist = reshape(HistExtract80(img0, [samplePos.sx',samplePos.sy'],ftrPos_hog),...
        size(ftrPos_hog,1), length(samplePos.sx),80);
    posx.feature_hog = feature_hog(:,1:length(posx.sampleImage.sx),:);
    negx.feature_hog = feature_hog(:,length(posx.sampleImage.sx)+1:end,:);
    posx.feature_hist = feature_hist(:,1:length(posx.sampleImage.sx),:);
    negx.feature_hist = feature_hist(:,length(posx.sampleImage.sx)+1:end,:);

    
    % update background template and target template
    [bgmodel, tarmodel, proimg, th] = TemplateUpdate(bgmodel, tarmodel, img0, proimg, initstate(2), initstate(1),[th i/5.0]);
    
    % get mask 
    mask = GetMask(proimg, ftrPos_hog, 0.5);
    
    if SHOW_PROCEDURE
        figure(3);
        set(gcf, 'position', [0,0,500,500]);
        subplot(1,3,1);imshow(imresize(bgmodel/255, size(bgmodel)/2));title('bgmodel');
        subplot(1,3,2);imshow(tarmodel/255);title('tarmodel');
        subplot(1,3,3);imshow(proimg);title('proimg');
        pause(0.1);
    end
    
    % train & select weak classifiers
    [clspool_hog, clspool_hist, posx, negx, selector] = selectClfMask(clspool_hog, clspool_hist, posx, negx, mask, mask);

end
times = toc;

save(['result/' benchmark '.mat'], 'trackpos', 'dist', 'ftrPos_hog', 'times');