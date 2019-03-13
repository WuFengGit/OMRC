function initstate = showDetection(frameNO, img, prob, sampleImage, imgsize, selector_hog, selector_hist, ftr_hog)
drawFea = true;
[c,index] = max(prob);
%-------------------------------------
x = sampleImage.sx(index);
y = sampleImage.sy(index);
w = sampleImage.sw(index);
h = sampleImage.sh(index);
initstate = [x y w h];

figure(1)
set(gcf, 'position', [500,0,500,500]);
imshow(imresize(img/255, size(img)/2));
rectangle('Position',initstate/2,'LineWidth',2,'EdgeColor','r');
text(5, 18, strcat('#',num2str(frameNO)), 'Color','y', 'FontWeight','bold', 'FontSize',20);
set(gca,'position',[0 0 1 1]);
% -------------------------------------------------  draw selected feature location
if drawFea
    %selector_hog = 1:150;
    for j=1:length(selector_hog)
        tmpinitstate = [ftr_hog(selector_hog(j),1)+x,ftr_hog(selector_hog(j),2)+y,ftr_hog(selector_hog(j),3),ftr_hog(selector_hog(j),4)];
        rectangle('Position',tmpinitstate/2,'LineWidth',1,'EdgeColor','r');
    end
    for j=1:length(selector_hist)
        tmpinitstate = [ftr_hog(selector_hist(j),1)+x,ftr_hog(selector_hist(j),2)+y,ftr_hog(selector_hist(j),3),ftr_hog(selector_hist(j),4)];
        rectangle('Position',tmpinitstate/2,'LineWidth',1,'EdgeColor','b');
    end
end
pause(0.1)

minp = min(prob);
maxp = max(prob);
prob = (prob-minp)/(maxp-minp);
img = zeros(imgsize);
for i=1:length(sampleImage.sx)
    img(sampleImage.sy(i),sampleImage.sx(i)) = prob(i);
end

figure(2)
set(gcf, 'position', [1000,0,500,500]);
imshow(imresize(uint8(img*255), size(img)/2));
rectangle('Position',initstate/2,'LineWidth',2,'EdgeColor','r');
title('problity')
pause(0.02);
end
