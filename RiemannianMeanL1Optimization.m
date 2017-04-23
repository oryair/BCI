function mM = RiemannianMeanL1Optimization(tC)

    f   = @(vU0) CalcVal(vU0, tC);
    mM  = RiemannianMean(tC);
%     mM  = mean(tC, 3);
    mL0 = chol(mM, 'lower');
    vL0 = mL0(tril(true(size(mL0))));
%     vL0 = vL0 + 0.01 * randn(size(vL0));
    
%     DDD = abs(mM - mU0 * mU0');
%     max(DDD(:))
    
    options = optimset('Display','iter','PlotFcns',@optimplotfval, 'MaxIter', 100000, 'MaxFunEvals', 100000);
    vL = fminsearch(f, vL0, options);
    
    mL       = tril(ones(size(tC,1)));
    mL(~~mL) = vL;
    mM       = mL * mL';

end

function fVal = CalcVal(vL, tC)

    mL       = tril(ones(size(tC,1)));
    mL(~~mL) = vL;
    mP       = mL * mL';
    
    fVal = 0;
    for ii = 1 : size(tC, 3)
%         ii
        fVal = fVal + RiemannianDist(mP, tC(:,:,ii), 1)^2;
    end

end