function fig = alba9b_plot(R,titleText,fileStem,zoomHalfWidth)
%ALBA9B_PLOT Full-reference figure: the quantities ROM3 could not provide.
%   Middle row is the actuator: torque against the ASSUMED envelope, electrical
%   phase against the +-90 deg domain, and the ACTUAL rotor speed. Bottom row is
%   the closure close-up. No panel carries two y scales.
if nargin < 2 || isempty(titleText), titleText = 'ALBA Stage 9 B: full reference model'; end
if nargin < 3 || isempty(fileStem),  fileStem  = 'B_full'; end
if nargin < 4 || isempty(zoomHalfWidth), zoomHalfWidth = 1.6; end

ix = alba9b_log_index();
B  = alba9b_const();
colors = {'#0d8a72','#c97a10','#596bb2'};
ink = '#3b3b3a'; muted = '#6f6f6d';

fig = figure('Color','w','Position',[60 40 1500 1080],'Visible','off');
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,9);
for k = 1:9, ax(k) = nexttile(tl); hold(ax(k),'on'); end

labels = strings(1,numel(R));
hMain = gobjects(1,numel(R));
tc = arrayfun(@(r)r.diag.t_close_s,R);
zlo = min(tc)-zoomHalfWidth; zhi = max(tc)+zoomHalfWidth;
tgt = arrayfun(@(r)r.scenario.target_mL,R);
distinct = numel(unique(tgt)) == numel(R);

for i = 1:numel(R)
    r = R(i); t = r.t; y = r.y; c = colors{mod(i-1,numel(colors))+1};
    labels(i) = sprintf('%d mL | %s',r.scenario.target_mL,r.diag.mode_name);

    plot(ax(1),t,y(:,ix.Qref),'--','Color',[.62 .62 .61],'LineWidth',1);
    hMain(i) = plot(ax(1),t,y(:,ix.Q),'Color',c,'LineWidth',2);
    plot(ax(2),t,y(:,ix.discharged),'Color',c,'LineWidth',2);
    plot(ax(2),t,y(:,ix.collected),':','Color',c,'LineWidth',1.6);
    plot(ax(3),t,y(:,ix.n),'Color',c,'LineWidth',2);
    if distinct
        text(ax(2),t(end)+0.3,y(end,ix.collected),sprintf('%.0f mL',r.diag.collected_mL), ...
            'Color',ink,'FontSize',9,'VerticalAlignment','middle');
    end

    plot(ax(4),t,abs(y(:,ix.torque)),'Color',c,'LineWidth',2);
    plot(ax(5),t,abs(y(:,ix.phi))*180/pi,'Color',c,'LineWidth',2);
    plot(ax(6),t,y(:,ix.omega)*60/(2*pi),'Color',c,'LineWidth',2);

    m = t >= zlo & t <= zhi;
    plot(ax(7),t(m),y(m,ix.Q),'Color',c,'LineWidth',2);
    plot(ax(7),t(m),y(m,ix.q_outlet),'--','Color',c,'LineWidth',1.6);
    plot(ax(8),t(m),abs(y(m,ix.torque)),'Color',c,'LineWidth',2);
    plot(ax(9),t(m),y(m,ix.Fn),'Color',c,'LineWidth',2);
    for k = [7 8 9], xline(ax(k),r.diag.t_close_s,':','Color',c,'LineWidth',1); end
end

yline(ax(1),15,':','Color',muted,'LineWidth',1);
setax(ax(1),'Time (s)','Valve flow (mL/s)','Flow: solved from motor input, never imposed',muted);
setax(ax(2),'Time (s)','Volume (mL)','Bladder discharge solid, outlet collected dotted',muted);
yline(ax(3),3000,':','Color',muted,'LineWidth',1);
setax(ax(3),'Time (s)','Command (STEP/s)','Commanded step rate',muted);

yline(ax(4),B.torque_cap_Nm,':','Color',muted,'LineWidth',1.2);
text(ax(4),0.985,0.93,'0.29659 N\cdotm ASSUMED flat envelope, not a measured curve', ...
    'Units','normalized','HorizontalAlignment','right','Color',muted,'FontSize',9);
setax(ax(4),'Time (s)','|Torque| (N\cdotm)','Torque against the conditional envelope',muted);
yline(ax(5),90,':','Color',muted,'LineWidth',1.2);
text(ax(5),0.985,0.93,'90 deg: beyond it the model no longer represents the step', ...
    'Units','normalized','HorizontalAlignment','right','Color',muted,'FontSize',9);
setax(ax(5),'Time (s)','|Phase error| (deg)','Electrical phase vs the domain limit',muted);
setax(ax(6),'Time (s)','Rotor speed (rpm)','ACTUAL rotor speed, solved not commanded',muted);

setax(ax(7),'Time (s)','Flow (mL/s)','CLOSE-UP: valve flow solid, final outlet dashed',muted);
yline(ax(8),B.torque_cap_Nm,':','Color',muted,'LineWidth',1.2);
setax(ax(8),'Time (s)','|Torque| (N\cdotm)','CLOSE-UP: torque demanded AFTER closure',muted);
setax(ax(9),'Time (s)','Contact force (N)','CLOSE-UP: force against a sealed bladder',muted);
for k = 7:9, xlim(ax(k),[zlo zhi]); end

hV = plot(ax(1),NaN,NaN,'-', 'Color',ink,'LineWidth',2);
hO = plot(ax(1),NaN,NaN,'--','Color',ink,'LineWidth',1.6);
lg = legend(ax(1),[hMain hV hO], ...
    cellstr([labels,"solid = valve flow","dashed = final outlet incl. tail"]), ...
    'Box','off','FontSize',10,'TextColor',ink,'Orientation','horizontal','NumColumns',3);
lg.Layout.Tile = 'south';

title(tl,titleText,'FontWeight','bold','FontSize',15);
subtitle(tl,['Full five-state reference, ported from the frozen source and verified pointwise ' ...
    '(worst relative 4e-19 on the rhs). Torque envelope and mechanical parameters are ASSUMPTIONS.'], ...
    'FontSize',10,'Color',muted);

here = fileparts(mfilename('fullpath'));
exportgraphics(fig,fullfile(here,'..',[fileStem '.png']),'Resolution',150);
exportgraphics(fig,fullfile(here,'..',[fileStem '.pdf']),'ContentType','vector');
fprintf('Wrote %s.png and .pdf\n',fileStem);
end

function setax(a,xl,yl,ti,muted)
xlabel(a,xl); ylabel(a,yl);
title(a,ti,'FontWeight','bold','FontSize',10);
grid(a,'on'); a.GridAlpha = 0.12; a.Box = 'off';
a.XColor = muted; a.YColor = muted; a.FontSize = 9;
end
