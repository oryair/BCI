function M = RiemannianMeanL1(tC)

Np = size(tC, 3);
M  = mean(tC, 3);

for ii = 1 : 200
    A = M ^ (1/2);      %-- A = C^(1/2)
    B = A ^ (-1);       %-- B = C^(-1/2)
        
    tS = zeros(size(tC));
    for jj = 1 : Np
        C = tC(:,:,jj);
        tS(:,:,jj) = A * logm(B * C * B) * A;
    end
    S = median(tS, 3);
    S = (S + S') / 2;
%     DDD = abs(S - S');
%     max(DDD(:))
    
    M = A * expm(B * S * B) * A; 
    
    eps = norm(S, 'fro')
    if (eps < 1e-6)
        break;
    end
end

end