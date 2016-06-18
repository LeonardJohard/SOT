function  [index, distance] = somdistance(vector1, vector2)
%vector1 = positions, vector2 = nodes
distance = [];
index = [];
for i = 1:size(vector1,2)
    relpos = vector2 - (vector1(:,i) * ones(1,length(vector2)));
    distance(i,:) = sum(relpos.^2, 1);
    [dummy, index(i)] = min(distance(i,:));
end