function D = earthMoversDistances(H1,H2)

D = zeros(size(H1,2),size(H2,2));
for i=1:size(H1,2)
    for j=1:size(H2,2)
        D(i,j) = earthMoverDist(H1(:,i),H2(:,j));
    end
end
