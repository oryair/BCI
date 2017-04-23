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
mW1 = zeros(288, 288);
for ii = 1 : 288
    ii
    for jj = ii + 1 : 288
        mW1(ii,jj) = RiemannianDist(tC(:,:,ii), tC(:,:,jj), 1);
    end
end
mW1 = mW1 + mW1';

%%
eps = 1 * median(mW1(:));
mK  = exp( -mW1.^2 / eps^2 );
mA  = bsxfun(@rdivide, mK, sum(mK, 2));
[mPhi, mLam] = eig(mA);

figure; scatter3(mPhi(:,2), mPhi(:,3), mPhi(:,4), 100, vClass, 'Fill'); colorbar; axis equal;

%%
mData = [mPhi; vClass'];