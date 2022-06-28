function PHId = expCorModel(rd,dim)
	PHId = zeros(dim);
    for i = 1:dim
        for j = 1:dim
            PHId(i,j) = rd^abs(i-j);
        end
    end
end

