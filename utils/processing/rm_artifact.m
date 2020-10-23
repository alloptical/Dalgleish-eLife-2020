function [out_traces] = rm_artifact(traces,artefact_start,artefact_end,varargin)
% Linearly interpolates over photostimulus artefacts intraces
% artefact_start and artefact_end are vectors of artefact start and stop
% times respectively

noise_flag = false;
if ~isempty(varargin)
    if strcmpi(varargin{1},'noise')
        noise_flag = true;
    end
end

out_traces = traces;
for i = 1:size(traces,1)
    for j = 1:numel(artefact_start)
        chunk                   = artefact_start(j):artefact_end(j);
        replace                 = linspace(traces(i,artefact_start(j)),traces(i,artefact_end(j)),numel(chunk));
        out_traces(i,chunk)     = replace;
    end
end

if noise_flag
    rng('default')
    noise_sd = std(diff(out_traces,[],2),[],2);
    n_traces = numel(noise_sd);
    for j = 1:numel(artefact_start)
        chunk = artefact_start(j):artefact_end(j);
        rand_chunk = 1 - (2 * rand(n_traces,numel(chunk)));
        out_traces(:,chunk) = out_traces(:,chunk) + (rand_chunk .* noise_sd);
    end
end


end

