function [st] = qstat(varargin)
% quick multi-purpose statistics function. Essentially pass in one or more
% columns or cells of data and enter stats info as keyword arguments (see
% below, e.g. 'paired' = true, 'normal' = true etc.)
%
% returns a stats structure with relevant info

% Initialise variables
paired = [];
normal = false;
stats = [];
c = [];
h = [];
p = [];
test = '';
tbl = [];
st.h = [];
st.p = [];
st.stats = [];
st.tbl = [];
st.c = [];
st.test = '';
st.normal = [];
st.paired = [];

% Parse inputs
idx = 0;
dat = {};
input_end = false;
for v = 1:numel(varargin)
    if input_end == false
        if iscell(varargin{v})
            for i = 1:numel(varargin{v})
                idx = idx+1;
                dat{idx} = varargin{v}{i}(:);
                dat{idx}(isnan(dat{idx})) = [];
            end
        elseif isnumeric(varargin{v})
            for i = 1:size(varargin{v},2)
                idx = idx+1;
                dat{idx} = varargin{v}(:,i);
                dat{idx}(isnan(dat{idx})) = [];
            end
        elseif ischar(varargin{v})
            input_end = true;
        end
    end
    if input_end == true
        if ischar(varargin{v})
            if strcmpi(varargin{v},'paired')
                if numel(unique(cellfun(@(x) numel(x),dat)))==1
                    paired = varargin{v+1}==1;
                else
                    paired = false;
                end
            end
        end
    end
end

% Deduce normality/paired
n = cellfun(@(x) numel(x),dat);
if isempty(paired)
    paired = numel(unique(n)) == 1;
end
if min(n) > 8
    normal = false(numel(dat),1);
    for i = 1:numel(dat)
        normal(i) = DagosPtest(dat{i});
    end
    normal = min(normal)==1;
else
    normal = false;
end

% Run appropriate test
if numel(dat)==1
    if normal
        test = 'ttest';
        [h,p] = ttest(dat{1});
    else
        test = 'signrank';
        [p,h] = signrank(dat{1});
    end
elseif numel(dat)==2
    if paired
        if normal
            test = 'ttest';
            [h,p] = ttest(dat{1},dat{2});
        else
            test = 'signrank';
            [p,h] = signrank(dat{1},dat{2});
        end
    else
        if normal
            test = 'ttest2';
            [h,p] = ttest2(dat{1},dat{2});
        else
            test = 'ranksum';
            [p,h] = ranksum(dat{1},dat{2});
        end
    end
else
    if paired
        dat_all = cell2mat(dat);
        test = 'friedman';
        [p,tbl,stats] = friedman(dat_all,1,'off');
    else
        group = [];
        dat_all = [];
        for i = 1:numel(dat)
            group = [group ; i*ones(numel(dat{i}),1)];
            dat_all = [dat_all ; dat{i}];
        end
        if normal
            test = 'anova1';
            [p,tbl,stats] = anova1(dat_all,group,'off');
        else
            test = 'kruskalwallis';
            [p,tbl,stats] = kruskalwallis(dat_all,group,'off');
        end
    end
    c = multcompare(stats,'CType','bonferroni','Display','off');
end

% Create statistics structure
st.h = h;
st.p = p;
st.stats = stats;
st.c = c;
st.normal = normal;
st.paired = paired;
st.test = test;
st.tbl = tbl;
st.mean = cellfun(@(x) mean(x(:)),dat);
st.sd = cellfun(@(x) std(x(:)),dat);
st.median = cellfun(@(x) median(x(:)),dat);
st.iqr = cellfun(@(x) iqr(x(:)),dat);

end

