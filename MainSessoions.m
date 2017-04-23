close all;
clear;

load BPF.mat;
dirPath = 'D:\Or\Workarea\Technion\BCI\BCICIV_2a_gdf\';

Fs  = 250;
Ts  = 1 / Fs;
T   = 3;
vT  = 0 : Ts : 7; vT(end) = [];

nFull    = length(vT);
startIdx = find(vT == 3);
endIdx   = find(vT == 6) - 1;
N        = endIdx - startIdx + 1;

%%
fileName = 'A01T.gdf';

[mS, H] = sload([dirPath, fileName]);

vTrig  = H.TRIG;
vClass = H.Classlabel;

tC = nan(22, 22, 288);
for ii = 1 : 288
    mEvent = mS(vTrig(ii) : (vTrig(ii) + nFull - 1), 1 : 22);
    mEvent = conv2(mEvent, BPF', 'same');
    mEvent = mEvent(startIdx : endIdx, :);
    
    mC         = cov(mEvent);
    tC(:,:,ii) = mC;
end

tC = tC(:,:,((vClass == 1) | (vClass == 2)));

mRiemannianMean = RiemannianMean(tC);
mCref           = mRiemannianMean^(-1/2);

K  = size(tC, 3);
M  = 22;
MM = M * (M + 1) / 2;
mX = zeros(MM, K);

mW = sqrt(2) * ones(M) - (sqrt(2) - 1) * eye(M);
for kk = 1 : K
    Skk      = logm(mCref * tC(:,:,kk) * mCref) .* mW;
    mX(:,kk) = Skk(triu(true(size(Skk))));
end

%%
% mMean  = mean(mX, 2);
% mX     = mX - mMean;
mCoeff = pca(mX');
mY     = mCoeff(:,1:2)' * mX;
mCov   = cov(mY');
vMu    = mCoeff(:,1:2)' * mean(mX, 2);
figure; error_ellipse(mCov, vMu); axis equal;
hold on; plot(vMu(1), vMu(2), 'b*');

%%
fileName = 'A01E.gdf';
load([dirPath, 'A01E.mat']);

[mS, H] = sload([dirPath, fileName]);

vTrig  = H.TRIG;
vClass = classlabel;

tC = nan(22, 22, 288);
for ii = 1 : 288
    mEvent = mS(vTrig(ii) : (vTrig(ii) + nFull - 1), 1 : 22);
    mEvent = conv2(mEvent, BPF', 'same');
    mEvent = mEvent(startIdx : endIdx, :);
    
    mC         = cov(mEvent);
    tC(:,:,ii) = mC;
end

tC = tC(:,:,((vClass == 1) | (vClass == 2)));

mRiemannianMean = RiemannianMean(tC);
% mCref           = mRiemannianMean^(-1/2);

K  = size(tC, 3);
mX = zeros(MM, K);
for kk = 1 : K
    Skk      = real(logm(mCref * tC(:,:,kk) * mCref)) .* mW;
    mX(:,kk) = Skk(triu(true(size(Skk))));
end

mY     = mCoeff(:,1:2)' * mX;
mCov   = cov(mY');
vMu    = mCoeff(:,1:2)' * mean(mX, 2);
hold on; error_ellipse(mCov, vMu); axis equal;
plot(vMu(1), vMu(2), 'r+');
