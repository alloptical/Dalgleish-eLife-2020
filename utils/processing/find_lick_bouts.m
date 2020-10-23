function [bout_starts,bout_n] = find_lick_bouts(lick_idcs,min_d,pre_buffer)
% Finds bouts of licking (from an array of lick times), defined as licks
% that occur within some minimum time-window of each other (min_d) with the
% first lick beginning some minimum time period after any preceding licks
% (pre_buffer).
%
% Note times are internally consistent (i.e. lick_idcs, min_d and
% pre_buffer should all be samples OR seconds OR ms etc.)

lick_d = [inf diff(lick_idcs)];
bout_idcs = find(lick_d>pre_buffer);
bout_starts = lick_idcs(bout_idcs);
bout_n = 0*bout_starts;
for i = 1:numel(bout_idcs)
    n_licks = find(lick_d(bout_idcs(i)+1:end)>min_d,1,'First');
    if ~isempty(n_licks)
        bout_n(i) = n_licks;
    end
end
end

