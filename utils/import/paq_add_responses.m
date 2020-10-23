function pa = paq_add_responses(pa,varargin)

% Records whether there was a response after each stimulus (within a response
% window) and what the reaction time was (in none then NaN) and adds to
% sync data structure
%
% Input: pa = sync data structure (_paqanalysis.mat)
% Optional input: 'response_window' = response window duration in seconds

response_window = 1.2;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'response_window')
        response_window = varargin{v+1};
    end
end

num_stims = numel(pa.stims);
vars_per_stim = zeros(1,num_stims);
for i = 1:num_stims
    vars_per_stim(i) = numel(pa.stims(i).vars);
end
stims = find(vars_per_stim>0);

for s = stims
    [pa.stims(s).rxntimes,pa.stims(s).responses] = deal(cell(numel(pa.licks),vars_per_stim(s)));
    for v = 1:vars_per_stim(s)
        for t = 1:numel(pa.stims(s).in.sec{v})
            for l = 1:numel(pa.licks)
                lick_idx = find(pa.licks(l).sec>pa.stims(s).in.sec{v}(t),1,'First');
                if ~isempty(lick_idx)
                    rxntime = pa.licks(l).sec(lick_idx) - pa.stims(s).in.sec{v}(t);
                    if rxntime <= response_window
                        pa.stims(s).responses{l,v}(t) = true;
                        pa.stims(s).rxntimes{l,v}(t) = rxntime;
                    else
                        pa.stims(s).responses{l,v}(t) = false;
                        pa.stims(s).rxntimes{l,v}(t) = nan;
                    end
                else
                    pa.stims(s).responses{l,v}(t) = false;
                    pa.stims(s).rxntimes{l,v}(t) = nan;
                end
            end
        end
    end
end

end