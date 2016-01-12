function dist_obj=gendist2object(dist)
if size(dist,1)==1
    dist={dist};
end
dist_obj=cell(length(dist),1);
for i=1:length(dist)
    dist_i=dist{i};
    if ~isa(dist_i, 'Distribution')
    switch dist_i{1}
        case 'beta'
            dist_obj{i}=BetaDistribution(dist_i{2}{:});
        case 'exponential'
            dist_obj{i}=ExponentialDistribution(dist_i{2}{:});
        case 'lognormal'
            dist_obj{i}=LogNormalDistribution(dist_i{2}{:});
        case 'normal'
            dist_obj{i}=NormalDistribution(dist_i{2}{:});
        case 'uniform'
            dist_obj{i}=UniformDistribution(dist_i{2}{:});
        %case 'translated'
        %    error(
        
    end
    else dist_obj{i}=dist_i;
    end
end
if length(dist_obj)==1
    dist_obj=dist_obj{1};
end