function meshadapt_step(indicator,refmax,reftol,crstol)
%MESHADAPT_STEP Single mesh adaptation step
%
%   MESHADAPT_STEP(INDICATOR,REFMAX,REFTOL,CRSTOL)
% 
%   performs a single mesh adaptation step based on the given element-wise
%   INDICATOR array and the prescribed tolerances REFTOL and CRSTOL.
%   Elements are only refined up the the maximum refinement level REFMAX.
%
% Input arguments:
%
%    INDICATOR     element-wise refinement/coarsening indicator
%    REFMAX        maximum admissible refinement level of the initial mesh
%    REFTOL        refinement tolerance
%    CRSTOL        coarsening tolerance
%
% Author: M. Moller, TU Delft, 2014.
    
    global lp_meshadapt
    
    % Perform one step of mesh adaptation
    calllib('meshadapt', 'madapt_step_dble2', lp_meshadapt, ...
            length(indicator), indicator(:), refmax, reftol, crstol);