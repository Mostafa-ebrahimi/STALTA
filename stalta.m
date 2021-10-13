function [ r ] = stalta(dtrace,slen,llen,varargin)
%stalta(trace,slen,llen)
%   calculates short term average-long term average ratio for trace
%   llen, slen are long and short windows in samples.
%   i.e., slen = seconds*samples per second.
%
%   If 4th argument is power, use amplitude squared.
%   Otherwise, default to absolute value of data trace

power = false;
if nargin > 3
    type = varargin{1};
    switch type
        case {'Power','power','POWER','Energy','energy'}
            power = true;
    end
end

if power
    dtrace = dtrace.^2;
else
    dtrace = abs(dtrace);
end
tmedian = median(dtrace);
npts = length(dtrace);
sta = dtrace*0 + tmedian;
lta = sta;
cumulative = cumsum(dtrace);

i = 1:(npts-slen);
sta(i) = (cumulative(i+slen) - cumulative(i))/slen;

i = (llen+1):npts;
lta(i) = (cumulative(i) - cumulative(i-llen))/llen;
r = sta./lta;

end

