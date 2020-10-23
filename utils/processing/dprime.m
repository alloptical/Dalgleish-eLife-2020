function [dp] = dprime(h_nr,h_nt,fa_nr,fa_nt)
% Calculates dprime
% Corrects 100%/0% according to Macmillan & Kaplan, 1985.

h_pc = h_nr./h_nt;
h_pc = adjust_pc(h_pc,h_nt);

fa_pc = fa_nr./fa_nt;
fa_pc = adjust_pc(fa_pc,fa_nt);

dp = norminv(h_pc) - norminv(fa_pc);

end

function pc = adjust_pc(pc,nt)

ceil_flag = pc==1;
floor_flag = pc==0;

pc(ceil_flag) = (nt(ceil_flag)-0.5) ./ nt(ceil_flag);
pc(floor_flag) = 0.5 ./ nt(floor_flag);

end

