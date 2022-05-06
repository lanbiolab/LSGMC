function MM = retain(M)
    N = size(M,2);
    MM = zeros(N,N);
    ro= 0.032+1/(0.018*N-1.42);
    [S,Ind] = sort(abs(M),1,'descend');
    for i = 1:N
        cL1 = sum(S(:,i));
        stop = false;
        cSum = 0; t = 0;
        while (~stop)
            t = t + 1;
            cSum = cSum + S(t,i);
            if ( cSum >= ro*cL1 )
                stop = true;
                MM(Ind(1:t,i),i) = M(Ind(1:t,i),i);
            end
        end
    end
end