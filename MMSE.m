function W_MMSE = MMSE(H,P)
    [K,~] = size(H);
    alpha = K/P;% regularization parameter
    W_MMSE = H'/(H*H'+alpha*eye(K));
    W_MMSE = sqrt(P)/norm(W_MMSE,'fro')*W_MMSE;    
end

