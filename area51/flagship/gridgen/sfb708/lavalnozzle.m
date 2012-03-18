function f = lavalnozzle(area, dist, segment)
%LAVALNOZZLE optimization of Laval nozzle for SFB 708
%
% area - 
    
    % Calculate angle of converging throat
    phi = calc_phi(1e-12);
    
    % Calculate radius of converging throat
    r = 8.85./(1+cos(phi));
    
    % Calculate angle of diverging throat
    ang = atan(0.75/48);
    
    % Generate structure with initial nozzle design
    S = struct('width1',   10.0,   ... % length of inlet chamber
               'width2',    0.0,   ... % length of inlet chamber
               'width3',   10.0,   ... % length of outlet chamber
               'height0',  11.15,  ... % center of inlet
               'height1',   0.0,   ... % height of inlet chamber
               'height2',   1.0,   ... % height of inlet chamber
               'height3',  10.0,   ... % height of outlet chamber
               'angle1',    phi,   ... % angle of circle segment
               'angle2',    ang,   ... % angle of diverging throat
               'radius0',   4.075, ... % radius of inlet
               'radius1',   r,     ... % radius of circle segment
               'x1',        0.0,   ... % x-coordinate of circle segment
               'y1',       13.0,   ... % y-coordinate of circle segment
               'filename','nozzle',... % name of temporal files
               'length',   56.85,  ... % length of nozzle
               'area',     area,   ...
               'dist',     dist,   ...
               'segment',  segment,...
               'pid',       0);
    
    % Generate initial nozzle design
    cmd = cmd_gridgen(S(1));
    s   = system(cmd);
        
    function phi = calc_phi(tol)
    %CALC_PHI calculates the angle of the converging throat
        
        f    = @(phi) 11.75./sin(phi).*(1+cos(phi))+8.85;
        df   = @(phi) -(11.75.*(1+cos(phi)).*cos(phi))./sin(phi).^2-11.75;
        phi0 = 3/2*pi;
        
        while 1
            phi  = phi0-f(phi0)./df(phi0);
            err  = abs(phi-phi0);
            phi0 = phi;
            if err < tol, return, end
        end        
    end
    
    function cmd = cmd_gridgen(s)
    %CMD_GRIDGEN generates the command string for grid generation

        cmd = ['gridgenlaval2d' ...
               ' -w1 '   num2str(s.width1) ...
               ' -w2 '   num2str(s.width2) ...
               ' -w3 '   num2str(s.width3) ...
               ' -h0 '   num2str(s.height0) ...
               ' -h1 '   num2str(s.height1) ...
               ' -h2 '   num2str(s.height2) ...
               ' -h3 '   num2str(s.height3) ...
               ' -ang1 ' num2str(s.angle1) ...
               ' -ang2 ' num2str(s.angle2) ...
               ' -r0 '   num2str(s.radius0) ...
               ' -r1 '   num2str(s.radius1) ...
               ' -x1 '   num2str(s.x1) ...
               ' -y1 '   num2str(s.y1) ...
               ' -l '    num2str(s.length) ...
               ' -f '    s.filename ...
               ' -a '    num2str(s.area) ...
               ' -d '    num2str(s.dist) ...
               ' -s '    num2str(s.segment) ...
               ' && rm -f ' s.filename '.1.* ' s.filename '.poly']
end
