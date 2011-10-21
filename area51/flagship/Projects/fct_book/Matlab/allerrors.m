function allerrors(dirname,extname)

files=dir([dirname '/*.' extname]);
fid=0;

for i=1:size(files),
    
    % Check if we start with a new sequence
    if ~isempty(findstr(files(i).name,'Nlev2')),
        % Close previous file?
        if fid~=0,
            % Write Latex table footer
            fprintf(fid, '\\hline\n');
            fprintf(fid, '\\end{tabular}\n');
            
            % Close output file
            fclose(fid);
        end
        % Open output file
        fid=fopen([dirname '/' strrep(files(i).name, '-Nlev2.csv', '') '.tex' ...
                  ], 'w');
        
        % Clear arrays
        Rho_err1=[];
        P_err1  =[];
        U_err1  =[];
        
        % Write Latex table header
        fprintf(fid, '\\begin{tabular}{llllll}\n');
        fprintf(fid, '\\hline\n');
        fprintf(fid, ['$E_1(\\rho)$ & $p(\\rho)$ & $E_1(u)$ & $p(u)$ & ' ...
                      '$E_1(p)$ & $p(p)$\\\\\n']);
        fprintf(fid, '\\hline\n');
    end
    
    % Compute L1-errors
    [rho_err1, p_err1, u_err1] = errors([dirname '/' files(i).name]);
    
    % File export
    if isempty(Rho_err1),
        % Write errors
        fprintf(fid, '%10.4e & %10.4e & %10.4e & %10.4e & %10.4e & %10.4e\\\\\n', ...
            rho_err1, 0, u_err1, 0, p_err1, 0);
    else
        % Write errors and convergence rate
        fprintf(fid, '%10.4e & %10.4e & %10.4e & %10.4e & %10.4e & %10.4e\\\\\n', ...
            rho_err1, log2(Rho_err1(end)/rho_err1), ...
            u_err1,   log2(U_err1(end)/u_err1), ...
            p_err1,   log2(P_err1(end)/p_err1));
    end
    
    % Append L1-erros
    Rho_err1 = [Rho_err1 rho_err1];
    P_err1   = [P_err1 p_err1];
    U_err1   = [U_err1 u_err1];
end

% Write Latex table footer
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');

% Close output file
fclose(fid);
