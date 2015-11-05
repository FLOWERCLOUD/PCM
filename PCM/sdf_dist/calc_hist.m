
%  Copyright (c) 2013
%      Noa Fish <noafish@post.tau.ac.il>

function h = calc_hist(data, mn, mx, nbins, farea)

sum_farea = sum(farea);

ndata = length(data);
h = zeros(nbins, 1);

if(mn == mx)
    h(1) = 1;
    return;
end

if(sum_farea == 0)
    return;
end

for i=1:ndata
    bin = ceil(nbins * (data(i) - mn) / (mx - mn));
    if(bin == 0)
        bin = 1;
    end
    if(bin == nbins+1)
        bin = nbins;
    end
    h(bin) = h(bin) + farea(i);
end

    
end
