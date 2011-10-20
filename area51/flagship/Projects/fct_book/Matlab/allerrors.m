function allerrors(dirname,extname)

files=dir([dirname '/*.' extname]);
fid=fopen([dirname '/convstudy.table'], 'w');

for i=1:size(files),
    
    % Check if we start with a new sequence
    if ~isempty(findstr(files(i).name,'Nlev2')),
        fprintf(fid, '\n\n%s\n', files(i).name);
        
        % Clear arrays
        Rho_err1=[];
        P_err1  =[];
        U_err1  =[];
    end
    
    % Compute L1-errors
    [rho_err1, p_err1, u_err1] = errors([dirname '/' files(i).name]);
    
    % File export
    if isempty(Rho_err1),
        % Write errors
        fprintf(fid, '%10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e\n', ...
            rho_err1, 0, u_err1, 0, p_err1, 0);
    else
        % Write errors and convergence rate
        fprintf(fid, '%10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e\n', ...
            rho_err1, log2(Rho_err1(end)/rho_err1), ...
            u_err1,   log2(U_err1(end)/u_err1), ...
            p_err1,   log2(P_err1(end)/p_err1));
    end
    
    % Append L1-erros
    Rho_err1 = [Rho_err1 rho_err1];
    P_err1   = [P_err1 p_err1];
    U_err1   = [U_err1 u_err1];
end

fclose(fid);