function W_ZF = ZF(H,P)
    W_ZF = H'/(H*H');
    W_ZF = sqrt(P)/norm(W_ZF,'fro')*W_ZF;    
end

