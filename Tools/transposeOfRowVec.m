function vec = transposeOfRowVec(vec)
    [a,b]=size(vec);
    if a<b
        vec=vec';
    end
end