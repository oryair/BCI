close all;
clear;

load BPF.mat;
% dirPath = 'D:\Or\Workarea\Technion\BCI\BCICIV_2a_gdf\';
dirPath = 'C:\Users\Oryair\Desktop\Workarea\BCI\BCICIV_2a_gdf\';

addpath('Connectivity\');

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
%     P          = CalcPairs(mEvent);
%     mC         = CalcSig(P);
    tC(:,:,ii) = mC;
end

% tC(:,:,((vClass == 1) | (vClass == 2))) = [];

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
% mData  = [mX; vClass'];
vClass2 = vClass;
% vClass2((vClass == 1) | (vClass == 2)) = [];
mData  = [mX; vClass2'];

% mData(:, ((vClass == 2) | (vClass == 3))) = [];


% mX(:,ii) = mC(triu(true(size(mC))));
% mData = [mX; H.Classlabel'];
% figure; plot(s(:, 1:3));

