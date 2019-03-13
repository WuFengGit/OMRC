function [pos, num] = hogFtr(clfparams,M)
% xywh
width = clfparams.width;
height = clfparams.height;

pos = zeros(1,4);


w = width;
h=height;
idx = 1;

pos(idx,:) = [1,1,w,h];
idx = idx +1;


lastsum = h+w;
L = 12;
while w>=L && h>=L
    w = w/2;
    h = h/2;
    if w<L
        w = L;
    end
    if h<L
        h = L;
    end
    if h+w==lastsum
        break;
    end
    
    for x = 1:w/4:width-w
        for y=1:h/4:height-h
            pos(idx,:) = [x,y,w,h];
            idx =idx+1;
            if idx>M
                num = M;
                return;
            end
        end
    end
    lastsum = h+w;
end
pos = floor(pos);
pos(pos<1) = 1;
num = idx - 1;


