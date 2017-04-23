close all;
clear;

load BPF.mat;
dirPath = 'D:\Or\Workarea\Technion\BCI\BCICIV_2a_gdf\';

fileName = 'A01T.gdf';

[mS, H] = sload([dirPath, fileName]);

vTrig  = H.TRIG;
vClass = H.Classlabel;

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
mRiemannianMean = RiemannianMean(tC);
mCSR            = mRiemannianMean^(-1/2);

K  = size(tC, 3);
M  = 22;
MM = M * (M + 1) / 2;
mX = zeros(MM, K);

mW = sqrt(2) * ones(M) - (sqrt(2) - 1) * eye(M);
for kk = 1 : K
    Skk      = logm(mCSR * tC(:,:,kk) * mCSR) .* mW;
    mX(:,kk) = Skk(triu(true(size(Skk))));
end

%%
mW  = squareform( pdist(mX') );
eps = 1 * median(mW(:));
mK  = exp( -mW.^2 / eps^2 );
mA  = bsxfun(@rdivide, mK, sum(mK, 2));
[mPhi, mLam] = eig(mA);

figure; scatter3(mPhi(:,2), mPhi(:,3), mPhi(:,4), 100, vClass, 'Fill'); colorbar; axis equal;

%%
mData = [mPhi; vClass'];