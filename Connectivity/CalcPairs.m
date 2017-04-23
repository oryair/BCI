function Pairs = CalcPairs(mS)

N = size(mS, 2);

Pairs{N} = zeros(size(mS, 1));

    for ii = 1 : N
       
        vS  = mS(:,ii);
        mD  = squareform( pdist(vS) );
        
        eps = median(mD(:));
        mK  = exp( -mD.^2 / eps^2 );
%         mK  = mK ./ sum(mK, 2);
        mK = bsxfun(@rdivide, mK, sum(mK,2));
        
        Pairs{ii} = mK;
        
    end

end