function mSig = CalcSig(Pairs)

    N           = length(Pairs);
    mAnorm      = zeros(N, N);
    mEigVecs    = zeros(size(Pairs{1},1),N);
    
    for ii = 1 : N
        [mEigVecs(:,ii), ~] = eigs(Pairs{ii}',1);
    end
    
    
    for ii = 1 : N
        mPi   = Pairs{ii};
        vPsii = mEigVecs(:,ii);
        
        for jj = (ii + 1) : N
            mPj     = Pairs{jj};
            vPsij   = mEigVecs(:,jj);

%             mAnorm(ii,jj) = norm(mPj * mPi' - mPi * mPj', 'fro');
%             mAnorm(ii,jj) = norm(mPj' * vPsii - mPi' * vPsij);

%             mAnorm(ii,jj) = norm(mPj' * vPsii - vPsij)^2 + ...
%                             norm(mPi' * vPsij - vPsii)^2;
            mAnorm(ii,jj) = norm(mPj' * vPsii - vPsii)^2 + ...
                            norm(mPi' * vPsij - vPsij)^2;
            mAnorm(ii,jj) = sqrt(mAnorm(ii,jj));
            
            
            mAnorm(jj,ii) = mAnorm(ii,jj);
        end
    end
    
    eps  = median(mAnorm(:));
    mSig = exp( -mAnorm.^2 / eps^2 );
    min(eig(mSig))
end