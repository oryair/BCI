close all;
% clear;

load BPF.mat;
dirPath = 'D:\Or\Workarea\Technion\BCI\BCICIV_2a_gdf\';

% fileName = 'A01T.gdf';
fileName = 'A02T.gdf';

[mS, H] = sload([dirPath, fileName]);

vTrig  = H.TRIG;
vClass = H.Classlabel;

% load([dirPath, 'A01E.mat']);
% vClass = classlabel;

Fs  = 250;
Ts  = 1 / Fs;
T   = 3;
vT  = 0 : Ts : 7; vT(end) = [];

nFull        = length(vT);
startIdx = find(vT == 3);
endIdx   = find(vT == 6) - 1;

N  = endIdx - startIdx + 1;

tC = nan(22, 22, 288);
for ii = 1 : 288
    mEvent = mS(vTrig(ii) : (vTrig(ii) + nFull - 1), 1 : 22);
    mEvent = conv2(mEvent, BPF', 'same');
    mEvent = mEvent(startIdx : endIdx, :);
    
    mC         = cov(mEvent);
    tC(:,:,ii) = mC;
end

%%
tC1 = tC(:,:,vClass == 1);
tC2 = tC(:,:,vClass == 2);

mMean1 = RiemannianMean(tC1);
mMean2 = RiemannianMean(tC2);

%%
N1     = size(tC1, 3);
N2     = size(tC2, 3);
vDist1 = nan(2, N1);
vDist2 = nan(2, N2);

for ii = 1 : N1
    vDist1(1,ii) = RiemannianDist(tC1(:,:,ii), mMean1);
    vDist1(2,ii) = RiemannianDist(tC1(:,:,ii), mMean2);
end
   
for ii = 1 : N2
    vDist2(1,ii) = RiemannianDist(tC2(:,:,ii), mMean1);
    vDist2(2,ii) = RiemannianDist(tC2(:,:,ii), mMean2);
end

%%
    
figure; hold on; grid on;
scatter(vDist1(1,:), vDist1(2,:), 100, 'r+');
scatter(vDist2(1,:), vDist2(2,:), 100, 'bo');
plot([0, 10], [0, 10], ':k', 'LineWidth', 2);
axis equal;

%%
sum(vDist1(1,:) < vDist1(2,:))
sum(vDist2(1,:) > vDist2(2,:))
