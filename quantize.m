function [xbar,misc] = quantize(x,varargin)
%------------------------------------------------------------------------
% Quantization of x. 
%
% INPUT
% -----
%
% Several options are available as optional parameters, listed below:
%
% - 'centroids' : centroids to use. Default is []. If non-empty, all other 
%                 parameters are ignored.
%
% If 'centroids' is not provided or empty, uniform quantization will be
% done. Parameters then are:
%
% - 'border_method': either QUANTILE (default) or ABSOLUTE
% - 'border_param': parameters for borders. depending on borders_method:
%       * border_method==ABSOLUTE: the actual min and max value to use
%       * border_method==QUANTILE: the quantiles to use for min and max.
%                                   default is [10,90]
% - 'control_method': either 'RATE' (default) or 'DISTORSION'
% - 'control_param': parameter for quantization step. depending on
%                    control_method:
%       * control_method==RATE: number of quantization cells
%       * control_method==DISTORSION: size of quantization cells.
% - 'display': 1 or 0, whether or not to display some figure
%
% OUTPUT
% ------
% - xbar: quantized signal
% - misc: a structure containing the following fields:
%       * centroids: the used centroids. min and max are the first and last
%                    ones, respectively
%       * group: for each data point, gives the index of the selected
%                centroid
%       * delta: only there if centroids were not provided: selected
%                quantization step
%--------------------------------------------------------------------------
%                                 Antoine Liutkus, Inria, 2016
%flattening
shape = size(x);
x=x(:);
%Parsing input and initialization
p = inputParser;
p.addParameter('centroids', [], @(x)isnumeric(x));
p.addParameter('borders_method', 'QUANTILE', @(x)ismember(upper(x),{'QUANTILE','ABSOLUTE'}));
p.addParameter('border_param', [0.1,0.9], @(x)isnumeric(x));
p.addParameter('control_method', 'RATE', @(x)ismember(upper(x),{'RATE','DISTORSION'}));
p.addParameter('control_param', 16, @(x)isnumeric(x));
p.addParameter('display', 0, @(x)isnumeric(x));
p.KeepUnmatched = true;
p.parse(varargin{:})
if isempty(p.Results.centroids)
    if strcmpi( p.Results.borders_method, 'QUANTILE')
        minmax = quantile(x,p.Results.border_param);
    else
        minmax = p.Results.border_param;
    end
    minimum = minmax(1);
    maximum = minmax(2);
    
    if strcmpi( p.Results.control_method, 'RATE')
        delta = (maximum-minimum)/p.Results.control_param;
    else
        delta = p.Results.control_param;
    end
    centroids = minimum:delta:maximum;
else
    centroids = p.Results.centroids;
    delta = nan
end
display   = p.Results.display;
[c,p] = sort(centroids);
[~,ic] = histc(x,[-inf,(c(1:end-1)+c(2:end))/2,inf]);
group = p(ic);
xbar = centroids(group);
if display
    nbins = 300;
    [HH,VV] = hist(x,nbins);
    [HHbar,VVbar] = hist(xbar,VV);
    figure
    bar(VV,HH/max(HH(:)),'b')
    hold on
    bar(VVbar,HHbar/max(HHbar(:)),'r')
    legend('original','quantized')
    grid on
    xlabel('value')
    ylabel('distribution')
    title('normalized histograms of original and quantized')
end
    
xbar=reshape(xbar,shape);
group = reshape(group,shape);
    
misc = struct;
misc.centroids = centroids;
misc.group = group;
if ~isnan(delta)
    misc.delta=delta;
end