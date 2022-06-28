%"MMSE precoding for multiuser MISO downlink transmission 
%with non-homogeneous user SNR conditions"
function W_RZF = RZF(H,P,path_loss)
    [K,~] = size(H);
    if(~exist('path_loss','var'))
        path_loss = unifrnd(-5,5,1,K);
    end
    alpha = sum(1./(10.^(path_loss/10)*P));% regularization parameter
    W_RZF = H'/(H*H'+alpha*eye(K));
    W_RZF = sqrt(P)/norm(W_RZF,'fro')*W_RZF;    
end


