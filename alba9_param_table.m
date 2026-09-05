function rep = alba9_param_table(stage8_mat)
%ALBA9_PARAM_TABLE Cross-check alba9_const.m against the frozen MAT and export
%the annotated parameter table.
%
%   Two jobs, both auditable:
%     1. CONSISTENCY. Every constant alba9_const.m carries is recomputed from
%        ALBA_Stage8_Results.mat by the same rule alba9_make_params used, and
%        compared bit-for-bit. Any mismatch is an error, not a warning.
%     2. TABLE. Every field of the MAT's plant and hydraulic blocks is joined
%        with its unit / source / assumption status. A field with no annotation
%        is reported as UNANNOTATED rather than silently dropped.
%
%   Integer-to-double compatibility is retained deliberately: the MAT stores
%   stiffness_powers and steps_per_rev as int64, and MATLAB integer arithmetic
%   would either error or silently round. See stage8_matlab_fix/.
here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(stage8_mat)
    stage8_mat = fullfile(here,'ALBA_Stage8_Results.mat');
end
S = load(stage8_mat);
p = S.plant_parameters;
h = S.hydraulic_parameters;
P = alba9_const();

% ---------------------------------------------------------------- 1. consistency
expected = struct( ...
    'Vsilicone',double(p.Vsilicone), 'Dref_inner',double(p.Dref_inner), ...
    'Dref_outer',double(p.Dref_outer), 'mu',double(p.mu), ...
    'cage_mass',double(p.cage_mass), 'gravity',double(p.gravity), ...
    'gravity_height_slope',double(p.gravity_height_slope), 'rho',double(p.rho), ...
    'water_viscosity',double(p.water_viscosity), 'pipe_diameter',double(p.pipe_diameter), ...
    'pipe_area',double(h.pipe_area_m2), 'pipe_length',double(h.pipe_length_m), ...
    'Rh',double(p.Rh), 'Rlam',double(h.R_pipe_laminar_Pa_s_m3), ...
    'Ih',double(h.Ih_Pa_s2_m3), 'Re_lam',double(h.Re_laminar), ...
    'Re_turb',double(h.Re_turbulent), 'roughness',double(h.pipe_absolute_roughness_m), ...
    'head_Pa',double(p.rho)*double(p.gravity)*double(p.head), ...
    'ku',2*pi*double(S.motor_radius_m)/double(S.steps_per_rev), ...
    'D0',double(p.D0), 'Vmin',double(p.Vmin), 'V0',double(p.V0), ...
    'contact_stiffness',double(p.contact_stiffness), ...
    'coef',reshape(double(S.stiffness_coefficients),1,[]), ...
    'powers',reshape(double(S.stiffness_powers),1,[]));

names = fieldnames(expected);
check = struct('name',{},'const_value',{},'mat_value',{},'exact_match',{},'abs_diff',{});
worst = 0;
for i = 1:numel(names)
    n = names{i};
    a = P.(n); b = expected.(n);
    same = isequal(a,b);
    d = max(abs(a(:)-b(:)));
    worst = max(worst,d);
    check(end+1,1) = struct('name',n,'const_value',a,'mat_value',b, ...
                            'exact_match',same,'abs_diff',d); %#ok<AGROW>
    if ~same
        error('alba9_param_table:mismatch', ...
            'alba9_const.%s does not match the frozen MAT (max abs diff %.3g)',n,d);
    end
end
extra = setdiff(fieldnames(P),names);
if ~isempty(extra)
    error('alba9_param_table:unchecked', ...
        'alba9_const carries unchecked fields: %s',strjoin(extra',', '));
end

% ---- the values the request calls out by name, re-asserted here -------------
must = { ...
    'V0',                 5e-4,                       0
    'Dref_outer',         0.07,                       0
    'Dref_inner',         0.066,                      0
    'Vsilicone',          2.9061826440808009e-5,      0
    'D0',                 0.10034658601199008,        0
    ... % Vmin is the WATER volume at a 50 mm outer diameter. The equivalent
    ... % sphere encloses water + silicone, so the shell volume comes back off.
    'Vmin',               pi/6*0.05^3-2.9061826440808009e-5, 1e-18
    'cage_mass',          0.270,                      0
    'head_Pa',            981,                        1e-9
    'Rh',                 1.6440083045256217e9,       0
    'pipe_diameter',      0.004,                      0
    'pipe_length',        0.36,                       1e-15
    'ku',                 4.9087385212340526e-6,      0};
for i = 1:size(must,1)
    n = must{i,1}; want = must{i,2}; tol = must{i,3};
    got = P.(n);
    if abs(got-want) > tol
        error('alba9_param_table:value','%s = %.17g, expected %.17g',n,got,want);
    end
end
% Values quoted in the request that live outside alba9_const
if abs(double(p.cable_gap)-0.003) > 0,  error('cable_gap must be 3 mm'); end
if abs(double(p.Dmin)-0.05) > 0,        error('Dmin must be 50 mm');     end
if abs(double(S.motor_radius_m)-0.0025) > 0, error('shaft radius must be 2.5 mm'); end
if double(S.steps_per_rev) ~= 3200,     error('steps_per_rev must be 3200'); end
% The mould wall follows from the two natural diameters; assert the reading.
if abs((double(p.Dref_outer)-double(p.Dref_inner))/2 - 0.002) > 1e-15
    error('the implied per-side MOULD wall must be 2 mm');
end

% ---------------------------------------------------------------- 2. table
rows = [tableFor(p,alba9_param_annotations(),'plant'); ...
        tableFor(h,alba9_hydraulic_annotations(),'hydraulic')];
rows = [rows; struct('block','stage8_rom','name','motor_radius_m','value_text', ...
            fmt(double(S.motor_radius_m)),'unit','m','status','given', ...
            'source','user / motor spec','note','Bare shaft winding radius.')];
rows = [rows; struct('block','stage8_rom','name','steps_per_rev','value_text', ...
            fmt(double(S.steps_per_rev)),'unit','STEP/rev','status','given', ...
            'source','simulation setting','note', ...
            'Simulation uses 3200 STEP/rev. Stored as int64 in the MAT; must be cast to double before arithmetic.')];
rows = [rows; struct('block','stage8_rom','name','ku','value_text',fmt(P.ku), ...
            'unit','m/STEP','status','derived','source','2*pi*r/steps_per_rev', ...
            'note','Cable taken up per commanded step.')];
rows = [rows; struct('block','stage8_rom','name','stiffness_coefficients','value_text', ...
            strjoin(arrayfun(@(x)fmt(x),P.coef,'UniformOutput',false),'; '), ...
            'unit','-','status','learned','source','frozen Stage-8 SINDy', ...
            'note','Read from the MAT at full precision. Never retyped, never refitted in Stage 9.')];
rows = [rows; struct('block','stage8_rom','name','stiffness_powers','value_text', ...
            strjoin(arrayfun(@(x)fmt(x),P.powers,'UniformOutput',false),'; '), ...
            'unit','-','status','learned','source','frozen Stage-8 SINDy', ...
            'note','Stored as int64 in the MAT; must be cast to double before use as an exponent.')];
rows = [rows; struct('block','stage8_rom','name','audited_force_range','value_text','0 to 118', ...
            'unit','N','status','learned','source','frozen Stage-8 audit', ...
            'note',['The range over which the LEARNED stiffness was checked positive. ' ...
                    'It is NOT a motor force limit and NOT a force saturation in the model. ' ...
                    'Exceeding it means the ROM is being used outside its audited domain, not that the hardware saturated.'])];

unann = {rows(strcmp({rows.status},'UNANNOTATED')).name};

% ---------------------------------------------------------------- 3. export
T = struct2table(rows);
writetable(T,fullfile(here,'parameter_table.csv'),'Encoding','UTF-8');
writeMarkdown(fullfile(here,'parameter_table.md'),rows,worst,numel(check),unann);

rep = struct('checked_constants',numel(check),'worst_abs_diff',worst, ...
             'all_exact',all([check.exact_match]),'rows',numel(rows), ...
             'unannotated',{unann},'source_mat',stage8_mat, ...
             'source_mat_sha256',sha256file(stage8_mat));
fid = fopen(fullfile(here,'parameter_consistency.json'),'w');
c = onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(struct('consistency',rep,'checks',check),'PrettyPrint',true),'char');
clear c

fprintf('Parameter check: %d constants, all exact (worst |diff| = %.3g).\n',numel(check),worst);
fprintf('Parameter table: %d rows -> parameter_table.csv / .md\n',numel(rows));
if ~isempty(unann)
    fprintf(2,'UNANNOTATED fields: %s\n',strjoin(unann,', '));
end
end

% -------------------------------------------------------------------------
function rows = tableFor(blk,ann,blockName)
rows = struct('block',{},'name',{},'value_text',{},'unit',{},'status',{},'source',{},'note',{});
names = fieldnames(blk);
for i = 1:numel(names)
    n = names{i};
    k = find(strcmp({ann.name},n),1);
    if isempty(k)
        a = struct('unit','?','status','UNANNOTATED','source','?','note', ...
                   'Present in the frozen MAT but not described in the annotation table.');
    else
        a = ann(k);
    end
    rows(end+1,1) = struct('block',blockName,'name',n,'value_text',valueText(blk.(n)), ...
        'unit',a.unit,'status',a.status,'source',a.source,'note',a.note); %#ok<AGROW>
end
end

function s = valueText(v)
if ischar(v) || isstring(v)
    s = char(v);
elseif islogical(v)
    s = string(v);  s = char(s);
elseif isnumeric(v)
    v = double(v);
    if isscalar(v)
        s = fmt(v);
    else
        s = strjoin(arrayfun(@(x)fmt(x),v(:).','UniformOutput',false),'; ');
    end
else
    s = '<non-scalar>';
end
end

function s = fmt(x)
s = sprintf('%.17g',x);
end

function writeMarkdown(file,rows,worst,nchk,unann)
fid = fopen(file,'w','n','UTF-8');
c = onCleanup(@() fclose(fid));
fprintf(fid,'# ALBA Stage 9 v2 - frozen parameter table\n\n');
fprintf(fid,'Generated by `alba9_param_table.m` from `ALBA_Stage8_Results.mat`. Nothing here is hand-copied.\n\n');
fprintf(fid,'`alba9_const.m` was cross-checked against the MAT: **%d constants, all exact, worst |difference| = %.3g**.\n\n',nchk,worst);
if isempty(unann)
    fprintf(fid,'Every field in the MAT carries an annotation.\n\n');
else
    fprintf(fid,'**UNANNOTATED fields:** %s\n\n',strjoin(unann,', '));
end
fprintf(fid,['**Status vocabulary.** `given` = design / nameplate / user-supplied, taken as given, ' ...
    'not independently measured here. `derived` = computed from other entries. `assumption` = chosen by ' ...
    'the modeller or inherited as a conditional prior, not identified. `learned` = fitted by the frozen ' ...
    'Stage-8 SINDy run against synthetic noise-free data, not a measured constant. `solver` = numerical ' ...
    'setting. `legacy` = present but unused by Stage 9, read the note. `metadata` = descriptive string.\n\n']);
fprintf(fid,'**No row in this table is an experimentally validated hardware value.**\n\n');
fprintf(fid,'| block | parameter | value | unit | status | source | note |\n');
fprintf(fid,'|---|---|---|---|---|---|---|\n');
for i = 1:numel(rows)
    r = rows(i);
    fprintf(fid,'| %s | `%s` | %s | %s | **%s** | %s | %s |\n', ...
        r.block,r.name,r.value_text,r.unit,r.status,r.source,strrep(r.note,'|','\|'));
end
end

function d = sha256file(f)
fid = fopen(f,'r');
c = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8');
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes);
d = sprintf('%02x',typecast(md.digest(),'uint8'));
end
