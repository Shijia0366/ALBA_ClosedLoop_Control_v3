function [S,vec] = alba9_scenario(target_mL,mode,varargin)
%ALBA9_SCENARIO Stage-9 v2 design scenario and its packed parameter vector.
%   Every number here is a DESIGN value. None is a measurement, a hardware
%   guarantee or a medical criterion. The actuator limits are the Stage-8
%   conditional envelope, not measured capability.
%
%   v2 keeps every v1 control value unchanged so that the effect of the
%   PROGRAM fixes can be read separately from any later control change.
%
%   Override any field:  alba9_scenario(300,3,'Kp',80,'Ki',300)
S = struct();
S.target_mL    = target_mL;   % collected at the FINAL outlet: 200 / 300 / 400 mL
S.tail_comp_mL = 2.0;         % tail compensated when choosing the closing threshold
S.tail_vol_mL  = 2.0;         % tail assumed to be delivered downstream
S.q_plateau    = 15.0;        % mL/s design plateau
S.q_close      = 3.0;         % mL/s design flow held into closure
S.rise_s       = 2.0;         % startup raised-cosine duration, s
S.R0_mL        = 28.5;        % remaining volume where deceleration begins
S.r_hold_mL    = 1.5;         % remaining volume held at q_close before closing
S.Kp           = 50.0;        % STEP/s per (mL/s)   - error is in mL/s, output in
S.Ki           = 600.0;       % STEP/s^2 per (mL/s)   STEP/s. The gains carry that
S.Kaw          = 10.0;        % 1/s                   conversion; do NOT also
                              %                       convert the error to SI.
S.n_max        = 3000.0;      % conditional STEP/s ceiling
S.a_max        = 6000.0;      % conditional STEP/s^2 ceiling, valve OPEN
S.a_brake      = 2000.0;      % conditional STEP/s^2 ceiling, valve CLOSED
S.tau_slew     = 0.005;       % s. COMMAND-SHAPING parameter that realises the
                              % slew limit smoothly. It is NOT a measured motor
                              % electrical or mechanical time constant.
S.mode         = mode;        % 1 fixed step rate | 2 pure PI | 3 feedforward+PI
S.n_fixed      = 606.94;      % STEP/s, mode 1 only: the geometric feedforward
                              % evaluated once at the initial 500 mL, so the flow
                              % then sags as Ax shrinks
S.latch_gain   = 1e5;         % 1/s. v2 valve latch rate. The valve flag trips at
                              % latch >= 0.5, i.e. 0.5/latch_gain = 5 us after the
                              % volume threshold. Deterministic, reported, and the
                              % price of guaranteeing no closed->open chatter.
S.tau_q0       = 1e-3;        % s. Time constant of the pre-closure flow tracker
                              % that is frozen at closure. The captured q0 is
                              % therefore a 1 ms lagged value, NOT the exact
                              % left limit of Q at the closing instant.

S.terminal_manager = 0;    % v3 optimisation only; zero preserves v2 control
S.stop_remaining_mL = 0;   % stop winding before closure, never fabricate water
S.stop_n_steps_s = 1e-3;   % numerical command-stop threshold, NOT power-off
S.stop_omega_rad_s = 1e-3; % actual mechanical settling criterion
S.stop_acc_m_s2 = 1e-3;    % reject instantaneous zero-speed crossings
S.stop_dwell_s = 0.1;      % criteria must hold continuously this long
for k = 1:2:numel(varargin)
    if ~isfield(S,varargin{k})
        error('alba9_scenario:unknownField','Unknown scenario field %s',varargin{k});
    end
    S.(varargin{k}) = varargin{k+1};
end

vec = [S.target_mL; S.tail_comp_mL; S.tail_vol_mL; S.q_plateau; S.q_close; ...
       S.rise_s; S.R0_mL; S.r_hold_mL; S.Kp; S.Ki; S.Kaw; S.n_max; S.a_max; ...
       S.a_brake; S.tau_slew; S.mode; S.n_fixed; S.latch_gain; S.tau_q0; ...
       S.terminal_manager; S.stop_remaining_mL; S.stop_n_steps_s; ...
       S.stop_omega_rad_s; S.stop_acc_m_s2; S.stop_dwell_s];
end
