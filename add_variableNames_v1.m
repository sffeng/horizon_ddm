function [a, b] = add_variableNames_v1(wg, hg, wb, hb, delta, delta2)

% delta = 0.05;
a = annotation('textbox', [0 hg(1) wg(1)-delta hb(1)], 'string', 'non-decision time, T_0');
a(2) = annotation('textbox', [0 hg(1)+hg(2)+hb(1) wg(1)-delta hb(1)], 'string', 'bias, \alpha');
a(3) = annotation('textbox', [0 sum(hg(1:3))+sum(hb(1:2)) wg(1)-delta hb(1)], 'string', 'threshold, \beta');
a(4) = annotation('textbox', [0 sum(hg(1:4))+sum(hb(1:3)) wg(1)-delta hb(1)], 'string', 'drift rate, \mu');
set(a, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'right', ...
    'verticalalignment', 'middle')

% delta2 = 0.06;
b = annotation('textbox', [wg(1) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'baseline');
b(2) = annotation('textbox', [wg(1)+wb(1)+wg(2) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of \DeltaR');
b(3) = annotation('textbox', [wg(1)+wb(1)+wg(2)+wb(2)+wg(3) sum(hg(1:4))+sum(hb(1:4))+delta2 wb(1) hg(end)-delta2], 'string', 'effect of \DeltaI');
set(b, 'fontsize', 16, 'linestyle', 'none', ...
    'horizontalalignment', 'center', ...
    'verticalalignment', 'middle')
