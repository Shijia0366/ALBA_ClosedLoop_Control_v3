function ix = alba9_log_index()
%ALBA9_LOG_INDEX Column map of the alba9_log signal matrix (v2, 22 channels).
%   Kept next to build_alba9_model so the Mux wiring and every reader agree.
ix = struct( ...
    'Qref',1, ...         % mL/s, commanded flow reference
    'Q',2, ...            % mL/s, flow THROUGH THE VALVE (an ODE state)
    'V',3, ...            % mL, water remaining in the bladder
    'discharged',4, ...   % mL, displaced from the bladder
    'F',5, ...            % N, contact generalised force in cable coordinates
    'n',6, ...            % STEP/s, applied command after slew limit and saturation
    'n_raw',7, ...        % STEP/s, controller output before slew and saturation
    'keff',8, ...         % N/m, learned effective stiffness at the current F
    'Ax',9, ...           % m^2, |dV/d(cable coordinate)|
    'valve_closed',10,... % 0 open, 1 closed. LATCHED: monotone by construction
    'q_tail',11, ...      % mL/s, downstream residual-water drain at the outlet
    'drained',12, ...     % mL, tail volume delivered so far
    'collected',13, ...   % mL, total collected at the final outlet
    'e_track',14, ...     % mL/s, Qref - Q. Always the TRUE error, including in
    ...                   % open-loop fixed-step-rate mode where it is not fed back
    'n_dot',15, ...       % STEP/s^2, the APPLIED command acceleration. Logged
    ...                   % rather than finite-differenced: differencing a
    ...                   % variable-step log across the slew limit reports
    ...                   % spurious values just above the envelope.
    'ni',16, ...          % STEP/s, PI integral state
    'q0',17, ...          % mL/s, tracked pre-closure flow, frozen at closure.
    ...                   % A tau_q0-lagged value, NOT the exact left limit of Q.
    'tclose',18, ...      % s, closing time latch (equals t while open, then frozen)
    'q_outlet',19, ...    % mL/s, flow at the FINAL outlet: valve flow while open,
    ...                   % tail drain after closure. Distinct from channel 2.
    'latch',20, ...       % monotone valve latch state in [0 1]
    'cmp_raw',21, ...     % the RAW comparator before latching. Kept deliberately:
    ...                   % it still chatters at the threshold, and logging it is
    ...                   % how we show the latch absorbs the chatter rather than
    ...                   % the chatter having been tuned away.
    'sat_flag',22, ...    % +1 command riding n_max, -1 riding 0, 0 otherwise
    'slew_flag',23);      % 1 while the command acceleration is being clipped
end
