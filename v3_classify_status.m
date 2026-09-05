function status=v3_classify_status(done,guard,atHorizon,executionEvent)
% Explicit completion precedence; a closed valve alone is never success.
if nargin<4,executionEvent='';end
if strcmp(executionEvent,'solver_failure'),status='solver_failure';
elseif strcmp(executionEvent,'manual_interruption'),status='manual_interruption';
elseif guard,status='boundary_triggered';
elseif done,status='normal_completion';
elseif atHorizon,status='timeout';
else,status='interrupted_or_unclassified_stop';end
end
