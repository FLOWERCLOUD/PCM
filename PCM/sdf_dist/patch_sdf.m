%%% calc patch SDF histograms
%%% output:
% sdf_str: a struct containing SDF histograms for all patches as well as
%          extra data like id's of flat patches


function sdf_str = patch_sdf(points, normals, patches, patches_to_do, outliers, idx, flat_strictness, sdf_str_)


npatches = length(patches);

flats = [];
flats2 = cell(npatches,1);

raysAngle = 20 / 180 * pi; % angle of the cone

nbins = 10;
pow = 0.5; % damping factor

dists = cell(npatches,1); % stores distances of rays for the SDF histogram
dists2 = dists; % stores index of patch that was hit

sdf = zeros(length(patches),nbins);
sdf_norm = zeros(length(patches),nbins);
sdf_norm2 = zeros(length(patches),nbins);

for i=1:npatches
    
    if(isempty(find(patches_to_do == i)))
        dists{i} = [NaN];
        continue;
    end
    
    source_points = patches{i};
    
    for k=1:length(source_points)
    
        src = source_points(k);
        srcp = points(src,:);
        srcn = -1 * normals(src,:);
        
        % calc vector between src and all other points
        vecs = bsxfun(@minus, points, srcp);
        unit_vecs = sqrt(sum(vecs.^2, 2));
        unit_vecs = bsxfun(@rdivide, vecs, unit_vecs);
        
        costheta = sum(bsxfun(@times, unit_vecs, srcn), 2);
        angs = acos(costheta); 
        
        % filter out the points with a vector that isn't within the
        % designated cone
        within_cone = find(angs <= raysAngle);
        within_cone = setdiff(within_cone, outliers);
        idx_cone = idx(within_cone);
        
        if(isempty(within_cone))
            dists{i} = [dists{i}; NaN];
        else
            
            tmp = vecs(within_cone, :);
            tmp = sqrt(sum(tmp.^2, 2));
            tmp = tmp .* cos(angs(within_cone));
            
            [srdist, sri] = sort(tmp);
            
            % points that are close enough to the closest point
            valid_p = find(srdist <= srdist(1) + 1e-2);
            
            line_plane_denom = sum(unit_vecs(within_cone,:) .* normals(within_cone,:), 2);
            line_plane_denom2 = line_plane_denom(sri);
            
            % points that are facing the wrong way
            neg = find(line_plane_denom2 <= 0);
            
            if(~isempty(neg))
                ind = max(length(valid_p), neg(1)-1);
                idx_cone2 = idx_cone(sri);
                idx_cone3 = idx_cone2(1:ind);
                % take median as the distance
                dst = median(srdist(idx_cone3 == i).^pow);
                dists{i} = [dists{i}; dst];
                
                if(isnan(dst))
                    ind = min(length(valid_p), neg(1)-1);
                    idx_cone2 = idx_cone(sri);
                    idx_cone3 = idx_cone2(1:ind);
                    hst = hist([idx_cone3 1:npatches], npatches);
                    [mxh, mxi] = max(hst);
                    if(mxh == 1)
                        dists2{i} = [dists2{i}; NaN];
                    else
                        dists2{i} = [dists2{i}; mxi];
                    end
                else
                    dists2{i} = [dists2{i}; i];
                end
            else
                % take median as the distance
                dst = median(tmp(idx_cone == i).^pow);
                dists{i} = [dists{i}; dst];
                if(isnan(dst))
                    hst = hist([idx_cone 1:npatches], npatches);
                    [mxh, mxi] = max(hst);
                    dists2{i} = [dists2{i}; mxi];
                else
                    dists2{i} = [dists2{i}; i];
                end
            end          
        end  
    end    
end


% update min and max dist values
mxdist = cellfun(@max, dists);
mxdist = max(mxdist);
mndist = cellfun(@min, dists);
mndist = min(mndist);

if exist('sdf_str_','var')
    mndist = min(sdf_str_.minmax(1), mndist);
    mxdist = max(sdf_str_.minmax(2), mxdist);
end

minmax2 = [mndist mxdist];

for i=1:npatches
    
    % if this patch is not on the list, copy its previous data
    if(isempty(find(patches_to_do == i)))
        sdf(i,:) = sdf_str_.sdf(i,:);
        sdf_norm(i,:) = sdf_str_.sdf_norm(i,:);
        sdf_norm2(i,:) = sdf_str_.sdf_norm2(i,:);
        continue;
    end
    
    d = dists{i};
    non_nan = find(~isnan(d));
    nd = length(non_nan);
    
    hst = calc_hist(d(non_nan), mndist, mxdist, nbins, ones(size(non_nan)));
    
    % filter out bins that are smaller or equal to 5
    hst(hst <= 5) = 0;
    %hst = clean_hist(hst);
    
    good_bins1 = find(hst > 0);
    good_bins2 = find(hst >= 0.1*length(patches{i}));
    
    if(flat_strictness == 1)
        if(isempty(good_bins2))
            flats = [flats i];
        end
    else
        if(isempty(good_bins1))
            flats = [flats i];
        end
    end
    
    % filter out the patches that were hit by the current patch (i) at
    % least 30% of the time (used later on with merges)
    nndists2 = dists2{i};
    nndists2 = nndists2(~isnan(nndists2));
    
    tmp = [nndists2' 1:npatches];
    tmp2 = hist(tmp, npatches);
    tmp2 = tmp2 - 1;
    tmp2(i) = 0;
    [srt, srti] = sort(tmp2, 'descend');
    valid = find(srt >= 0.3*length(patches{i}));
    segs = srti(valid);
    perc = srt(valid) / length(patches{i});
    flats2{i} = [segs; perc];

    sdf(i,:) = hst;
    sdf_norm(i,:) = hst / sum(hst);
    sdf_norm2(i,:) = hst / trapz(hst);
end

flats2 = {flats2};

sdf_str = struct('sdf', sdf, 'sdf_norm', sdf_norm, 'sdf_norm2', sdf_norm2, 'flats', flats, 'flats2', flats2, 'minmax', minmax2);

end



function hst = clean_hist(hist)

hst = hist;
[mh, imh] = max(hist);

clean = 0;
for i=imh+1:length(hist)
    if(clean)
        hst(i) = 0;
    end
    if(hist(i) == 0)
        clean = 1;
    end
end

clean = 0;
for i=imh-1:-1:1
    if(clean)
        hst(i) = 0;
    end
    if(hist(i) == 0)
        clean = 1;
    end
end

end


