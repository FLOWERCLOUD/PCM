function x=earthMoverDist(hist1,hist2)

x=0;
for i=1:length(hist1)-1
    x=x+abs(hist2(i)-hist1(i));
    if hist2(i)>hist1(i)
        hist2(i+1)=hist2(i+1)+(hist2(i)-hist1(i));
        hist2(i)=hist1(i);
    else
        hist1(i+1)=hist1(i+1)+(hist1(i)-hist2(i));
        hist1(i)=hist2(i);
    end
end;