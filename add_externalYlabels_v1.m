function [a] = add_externalYlabels_v1(wg, hg, wb, hb, delta, str)

% delta = 0.05;
for i = 1:length(str)
    a(i) = annotation('textbox', [0 sum(hg(1:i))+sum(hb(1:i-1)) wg(1)-delta hb(i)], 'string', str{i});
end
set(a, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'right', ...
    'verticalalignment', 'middle')

