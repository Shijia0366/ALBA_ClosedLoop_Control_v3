function fig = alba9_plot(R,titleText,fileStem,zoomHalfWidth)
%ALBA9_PLOT Stage-9 v2 figure: full run plus a closure close-up.
%
%   Nine panels, one colour per run held across all of them, so identity never
%   depends on which panel you are reading. No panel carries two y scales.
%
%   The bottom row is the closure close-up the report asks for. Actual rotor
%   speed and motor torque are NOT plotted because ROM3 does not have them:
%   the fast motor coordinates were eliminated when it was identified. Those
%   panels appear only in the full five-state reference replay.
%
%   Palette: the project figure palette re-stepped so every pair clears the
%   colour-vision separation, chroma and contrast checks.
if nargin < 2 || isempty(titleText),   titleText = 'ALBA Stage 9 v2'; end
if nargin < 3 || isempty(fileStem),    fileStem  = 'ALBA_Stage9_v2';  end
if nargin < 4 || isempty(zoomHalfWidth), zoomHalfWidth = 1.6; end

ix = alba9_log_index();
colors = {'#0d8a72','#c97a10','#596bb2'};
ink = '#3b3b3a'; muted = '#6f6f6d';

fig = figure('Color','w','Position',[60 40 1500 1080],'Visible','off');
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,9);
for k = 1:9, ax(k) = nexttile(tl); hold(ax(k),'on'); end

tgt = arrayfun(@(r)r.scenario.target_mL,R);
distinct = numel(unique(tgt)) == numel(R);
labels = strings(1,numel(R));
hMain = gobjects(1,numel(R));
hValve = gobjects(1); hOutlet = gobjects(1);

tc = arrayfun(@(r)r.diag.t_close_s,R);
zlo = min(tc)-zoomHalfWidth; zhi = max(tc)+zoomHalfWidth;

for i = 1:numel(R)
    r = R(i); t = r.t; y = r.y; c = colors{mod(i-1,numel(colors))+1};
    labels(i) = sprintf('%d mL | %s',r.scenario.target_mL,r.diag.mode_name);

    % ---- row 1: whole run
    plot(ax(1),t,y(:,ix.Qref),'--','Color',[.62 .62 .61],'LineWidth',1);
    hMain(i) = plot(ax(1),t,y(:,ix.Q),'Color',c,'LineWidth',2);
    plot(ax(2),t,y(:,ix.discharged),'Color',c,'LineWidth',2);
    plot(ax(2),t,y(:,ix.collected),':','Color',c,'LineWidth',1.6);
    plot(ax(3),t,y(:,ix.n),'Color',c,'LineWidth',2);
    if distinct
        text(ax(2),t(end)+0.3,y(end,ix.collected),sprintf('%.0f mL',r.diag.collected_mL), ...
            'Color',ink,'FontSize',9,'VerticalAlignment','middle');
    end

    % ---- row 2: whole run
    plot(ax(4),t,y(:,ix.F),'Color',c,'LineWidth',2);
    plot(ax(5),t,y(:,ix.valve_closed),'Color',c,'LineWidth',2);
    plot(ax(5),t,y(:,ix.latch),':','Color',c,'LineWidth',1.4);
    plot(ax(6),t,y(:,ix.e_track),'Color',c,'LineWidth',1.6);

    % ---- row 3: closure close-up
    m = t >= zlo & t <= zhi;
    plot(ax(7),t(m),y(m,ix.Q),'Color',c,'LineWidth',2);
    plot(ax(7),t(m),y(m,ix.q_outlet),'--','Color',c,'LineWidth',1.6);
    plot(ax(8),t(m),y(m,ix.n),'Color',c,'LineWidth',2);
    plot(ax(9),t(m),y(m,ix.F),'Color',c,'LineWidth',2);
    for k = [7 8 9]
        xline(ax(k),r.diag.t_close_s,':','Color',c,'LineWidth',1);
    end
end

yline(ax(1),15,':','Color',muted,'LineWidth',1);
setax(ax(1),'Time (s)','Valve flow (mL/s)','Flow: reference dashed grey, achieved solid',muted);
setax(ax(2),'Time (s)','Volume (mL)','Bladder discharge solid, outlet collected dotted',muted);
yline(ax(3),3000,':','Color',muted,'LineWidth',1);
setax(ax(3),'Time (s)','Command (STEP/s)','Applied command after slew limit and saturation',muted);
yline(ax(4),118,':','Color',muted,'LineWidth',1);
setax(ax(4),'Time (s)','Contact force (N)','Contact force vs the 118 N audited range',muted);
setax(ax(5),'Time (s)','State','Valve flag solid, monotone latch dotted',muted);
ylim(ax(5),[-0.08 1.15]);
setax(ax(6),'Time (s)','Qref - Q (mL/s)','Tracking error, logged even in open loop',muted);

setax(ax(7),'Time (s)','Flow (mL/s)','CLOSE-UP: valve flow solid, final outlet dashed',muted);
setax(ax(8),'Time (s)','Command (STEP/s)','CLOSE-UP: command keeps braking after closure',muted);
setax(ax(9),'Time (s)','Contact force (N)','CLOSE-UP: force rises again against a sealed bladder',muted);
for k = 7:9, xlim(ax(k),[zlo zhi]); end

% neutral proxies: these two entries describe LINE STYLE, not a series, so they
% must not wear any one run's colour
hValve  = plot(ax(1),NaN,NaN,'-', 'Color',ink,'LineWidth',2);
hOutlet = plot(ax(1),NaN,NaN,'--','Color',ink,'LineWidth',1.6);
lg = legend(ax(1),[hMain hValve hOutlet], ...
    cellstr([labels, "solid = valve flow", "dashed = final outlet incl. tail"]), ...
    'Box','off','FontSize',10,'TextColor',ink,'Orientation','horizontal','NumColumns',3);
lg.Layout.Tile = 'south';

title(tl,titleText,'FontWeight','bold','FontSize',15);
subtitle(tl,['Frozen three-state ROM3, ideal noise-free flow feedback. Rotor speed and motor ' ...
    'torque are absent by construction - ROM3 has no motor states. Not experimental validation.'], ...
    'FontSize',10,'Color',muted);

here = fileparts(mfilename('fullpath'));
exportgraphics(fig,fullfile(here,[fileStem '.png']),'Resolution',150);
exportgraphics(fig,fullfile(here,[fileStem '.pdf']),'ContentType','vector');
fprintf('Wrote %s.png and .pdf\n',fullfile(here,fileStem));
end

function setax(a,xl,yl,ti,muted)
xlabel(a,xl); ylabel(a,yl);
title(a,ti,'FontWeight','bold','FontSize',10);
grid(a,'on'); a.GridAlpha = 0.12; a.Box = 'off';
a.XColor = muted; a.YColor = muted; a.FontSize = 9;
end
