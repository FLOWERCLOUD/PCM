%%% input:
% points: 3D point set
% normals: 3D normals corresponding to points
% patches: cell array containing the candidate points of all patches
% patch_to_do: id's of patches that need SDF histogram calculation
% outliers: outlier points
% idx: index of patch containing each point
% flat_strictness: 0: a patch is considered flat if it has an empty SDF
%                     histogram
%                  1: a patch is considered flat if it doesn't have a bin
%                     that is at least 10% of its size
%%%output:
% dist_matrix: distance matrix between patches
function dist_matrix = sdf_dist_matrix( points, normals, patches, patch_to_do, outliers, idx, flat_strictness )
    conv_str = patch_sdf( points, normals, patches, patch_to_do, outliers, idx, flat_strictness);
    sdf_norm = conv_str.sdf_norm;
    % create distance matrix
    patch_dists = earthMoversDistances(sdf_norm',sdf_norm');
    patch_dists2 = patch_dists / max(max(patch_dists));
    dist_matrix = patch_dists2 + eye(size(patch_dists2));
end