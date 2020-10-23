function [varargout] = title_stats(st,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

plotFlag = true;
if ~isempty(varargin)
    plotFlag = varargin{1};
end

if isfield(st,'h')
    str = {[st.test ': h = ' num2str(st.h)] ['p = ' num2str(st.p)]};
elseif isfield(st,'rsq')
    str = {['R2 = ' num2str(st.rsq)] ['p = ' num2str(st.corr_p)]};
end
if plotFlag
    t = get(gca,'Title');
    if isempty(t.String)
        title(str)
    else
        str = [str{1} ', ' str{2}];
        if iscell(t.String)
            t.String = [t.String {str}];
        else
            t.String = {t.String str};
        end
    end
end
if nargout == 1
    str = [str{1} ', ' str{2}];
    varargout{1} = str;
end

end

