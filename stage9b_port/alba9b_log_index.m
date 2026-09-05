function ix = alba9b_log_index()
%ALBA9B_LOG_INDEX Column map of the full-reference closed-loop log (32 channels).
%   Channels 24-32 are the ones ROM3 cannot provide at all: electrical phase,
%   torque, ACTUAL rotor speed, cable tension, the model-validity guards and the
%   independently integrated flow.
ix = struct( ...
    'Qref',1, ...          % mL/s, commanded flow reference
    'Q',2, ...             % mL/s, flow through the valve (a solved state)
    'V',3, ...             % mL, bladder water
    'discharged',4, ...    % mL
    'Fn',5, ...            % N, cable-coordinate contact generalized force, NOT one button's normal force
    'n',6, ...             % STEP/s, applied command
    'n_raw',7, ...         % STEP/s, controller output before shaping
    'x',8, ...             % m, cable displacement
    'v',9, ...             % m/s, ACTUAL cable speed (solved, never imposed)
    'valve_closed',10, ... % latched valve flag
    'q_tail',11, ...       % mL/s
    'drained',12, ...      % mL
    'collected',13, ...    % mL
    'e_track',14, ...      % mL/s
    'n_dot',15, ...        % STEP/s^2, applied command acceleration
    'ni',16, ...           % STEP/s, PI integral state
    'q0',17, ...           % mL/s, captured pre-closure flow
    'tclose',18, ...       % s
    'q_outlet',19, ...     % mL/s, final outlet: valve flow, then tail
    'latch',20, ...        % monotone valve latch
    'cmp_raw',21, ...      % raw threshold comparator
    'sat_flag',22, ...
    'slew_flag',23, ...
    'phi',24, ...          % rad, ELECTRICAL PHASE ERROR. |phi| < pi/2 is the
    ...                    % declared domain; beyond it the torque relation turns
    ...                    % over and the frozen model does not represent stall.
    'torque',25, ...       % N*m, torque_cap*sin(phi)
    'omega',26, ...        % rad/s, ACTUAL rotor speed = v/r
    'tension',27, ...      % N, cable tension. A cable cannot push.
    'guard_violations',28, ... % count of failed guards, NOT mixed-unit minimum
    'guard_first_failed',29, ... % source-order index, NOT closest-boundary ranking
    'intQ',30, ...         % mL, flow integrated as its OWN state, for an
    ...                    % independent mass balance against the bladder
    'count',31, ...        % steps, commanded position = integral of n
    'acc',32, ...          % m/s^2, actual cable acceleration
    'guards',33:39, ...    % rad, m, m^3, m, rad/s, N, m^3/s respectively
    'overlap_raw',40, ...  % m, signed x-y (negative means gap)
    'pressure',41, ...     % Pa, bladder outlet gauge pressure
    'V_est',42, ...        % mL, V0 - independently integrated ideal flow
    'quiet_dwell',43, ...  % s, consecutive mechanically quiet interval
    'done',44, ...         % explicit normal-completion condition
    'terminal_stop',45, ...% separately labelled optional early-stop request
    'quiet',46);           % actual velocity AND acceleration AND command small
end
