function[first_aS_id] = getfirstDividingAS(T)

candidates = sort(find(strcmp(T,'A')));
for i=1:length(candidates)
    if any(strcmp(T.get(T.getchildren(candidates(i))),'A')) || any(strcmp(T.get(T.getchildren(candidates(i))),'T'))
        first_aS_id = candidates(i);
        break;
    end
end

end