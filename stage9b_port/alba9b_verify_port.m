function rep = alba9b_verify_port(jsonfile)
%ALBA9B_VERIFY_PORT Pointwise verification of the MATLAB port of the FULL
%five-state reference model against the frozen Python source.
%
%   The MATLAB port is not allowed to be "close in spirit". Every observable,
%   every right-hand-side component, every Jacobian entry and every guard is
%   compared at 4800 states spanning water volume, flow, cable speed, contact
%   overlap, electrical phase and both valve states, plus a static loss sweep.
%
%   The tolerance is the same relative bound the Stage-8 package uses for its own
%   regression_checks: |a-b| / (1+|b|) < 1e-8 for the right-hand side, < 1e-7 for
%   the Jacobian. Anything worse is an error, not a warning.
%
%   Passing this means the MATLAB model IS the frozen model, expressed in another
%   language. It says nothing about whether the frozen model matches hardware.
here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(jsonfile)
    jsonfile = fullfile(here,'reference_points.json');
end
d = jsondecode(fileread(jsonfile));
P = alba9_const();
B = alba9b_const();
pts = d.points;

worst = struct('obs',0,'rhs',0,'jac',0,'guard',0,'loss',0);
whichWorst = struct('obs','','rhs','','jac','','guard','','loss','');
obsFields = {'phi','torque','Fm','Fn','fr','acc','tension','loss','slope', ...
             'Qdot','pe','pg','A','D','y'};

for i = 1:numel(pts)
    q = pts(i);
    o = alba9b_ref_obs(q.x,q.v,q.V,q.Q,q.count,q.sealed,P,B);
    for k = 1:numel(obsFields)
        f = obsFields{k};
        e = rel(o.(f),q.(f));
        if e > worst.obs, worst.obs = e; whichWorst.obs = f; end
    end
    dz = alba9b_ref_rhs(q.x,q.v,q.V,q.Q,q.count,q.sealed,P,B);
    e = max(rel(dz,q.rhs(:)));
    if e > worst.rhs, worst.rhs = e; whichWorst.rhs = sprintf('point %d',i); end

    J = alba9b_ref_jac(q.x,q.v,q.V,q.Q,q.count,q.sealed,P,B);
    e = max(max(rel(J,reshape(q.jac,6,6))));
    if e > worst.jac, worst.jac = e; whichWorst.jac = sprintf('point %d',i); end

    g = alba9b_guards(q.x,q.v,q.V,q.Q,q.count,q.sealed,P,B);
    e = max(rel(g,q.guards(:)));
    if e > worst.guard, worst.guard = e; whichWorst.guard = sprintf('point %d',i); end
end

for i = 1:numel(d.loss)
    L = d.loss(i);
    [lo,sl] = alba9b_loss(L.Q_m3_s,P);
    e = max(rel([lo;sl],[L.loss_Pa;L.slope]));
    if e > worst.loss, worst.loss = e; whichWorst.loss = sprintf('Q = %.4g',L.Q_m3_s); end
end

% the ported geometry and loss must also still agree with the ROM's own copies,
% otherwise the two halves of Stage 9 would be running different physics
gd = 0;
for VmL = [500 400 300 200 100 50]
    V = VmL*1e-6;
    [D1,A1,pe1,pg1] = alba9_geometry(V,P);
    [D2,A2,pe2,pg2] = alba9b_geometry(V,P);
    gd = max([gd, rel([D1;A1;pe1;pg1],[D2;A2;pe2;pg2]).']);
end
ld = 0;
for QmLs = [0 1 5 10 15 20]
    Q = QmLs*1e-6;
    ld = max(ld, rel(alba9_loss(Q,P), alba9b_loss(Q,P)));
end

rep = struct('points',numel(pts),'loss_points',numel(d.loss), ...
    'worst_relative_obs',worst.obs,'worst_obs_field',whichWorst.obs, ...
    'worst_relative_rhs',worst.rhs,'worst_relative_jac',worst.jac, ...
    'worst_relative_guard',worst.guard,'worst_relative_loss',worst.loss, ...
    'rom_vs_port_geometry',gd,'rom_vs_port_loss',ld, ...
    'tolerance_rhs',1e-8,'tolerance_jac',1e-7, ...
    'scope','numerical port verification only; says nothing about hardware');

fprintf('Port verification over %d states and %d loss points:\n',numel(pts),numel(d.loss));
fprintf('  observables  worst relative %.3e   (%s)\n',worst.obs,whichWorst.obs);
fprintf('  rhs          worst relative %.3e\n',worst.rhs);
fprintf('  jacobian     worst relative %.3e\n',worst.jac);
fprintf('  guards       worst relative %.3e\n',worst.guard);
fprintf('  loss curve   worst relative %.3e\n',worst.loss);
fprintf('  ROM vs port: geometry %.3e, loss %.3e\n',gd,ld);

assert(worst.obs   < 1e-8, 'ported observables exceed the 1e-8 relative bound');
assert(worst.rhs   < 1e-8, 'ported rhs exceeds the 1e-8 relative bound');
assert(worst.jac   < 1e-7, 'ported jacobian exceeds the 1e-7 relative bound');
assert(worst.guard < 1e-8, 'ported guards exceed the 1e-8 relative bound');
assert(worst.loss  < 1e-8, 'ported loss exceeds the 1e-8 relative bound');
assert(gd < 1e-12 && ld < 1e-12, 'the ROM and the port disagree on shared physics');

fid = fopen(fullfile(here,'port_verification.json'),'w');
c = onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(rep,'PrettyPrint',true),'char');
clear c
fprintf('PASSED. port_verification.json written.\n');
end

function e = rel(a,b)
e = abs(a-b)./(1+abs(b));
end
