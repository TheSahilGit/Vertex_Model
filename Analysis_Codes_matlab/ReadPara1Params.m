function p = ReadPara1Params(filename)
% READPARA1PARAMS  Robustly parse para1_in.dat into a named struct.
%
%   p = ReadPara1Params("../para1_in.dat")
%   dt = p.dt; it_dump = p.it_dump;
%
% Why this exists instead of readtable(): several existing scripts do
%   para1 = readtable("../para1_in.dat");
%   dt = table2array(para1(8,1));
% This "magic row number" convention is fragile and non-obvious -- verified
% directly that MATLAB's readtable() silently drops the first two lines of
% this file via its header/description auto-detection, so table row N
% actually corresponds to file line N+2 (e.g. row 8 = line 10 = dt). That
% happens to be why "para1(8,1)" gives the right answer today, but it is
% not a mapping anyone should have to re-derive or trust across MATLAB
% versions or file-format tweaks.
%
% This function instead reads para1_in.dat line-by-line, in the EXACT
% order Fortran's allocation.f90::read_input consumes them (the two lists
% must be kept in sync if that subroutine's read order ever changes), and
% returns a struct indexed by parameter name -- e.g. p.dt, p.it_dump,
% p.totT, p.if_Shear_tissue -- so callers never depend on a row number.

names = { ...
    'nrun','nrun2_initialTime','Ao','Co','lambda','beta','gamm','eta', ...
    'totT','dt','if_Do_T1','min_d_T1','if_Do_T2','min_area_T2', ...
    'if_Fixed_boundary','if_bottom_borders_fixed','if_top_borders_fixed', ...
    'it_dump','T1_time_interval','T2_time_interval', ...
    'if_motility','if_motility_gradient','etas_max','etas_min','mot_Lc', ...
    'if_motility_decay','motility_decay_timeScale','if_motility_hotspot', ...
    'number_of_hotspot','sigma_hotspot', ...
    'if_Shear_tissue','if_Sudden_Shearing','sudden_shearStrength','sudden_shearWhen', ...
    'if_Oscillatory_Shearing','Oscl_shearStrength','Oscl_shearWhen','Oscl_freq_wo', ...
    'if_Perturb_tissue','if_sin_perturb','sin_perturb_when','sin_perturb_strength', ...
    'sin_perturb_waveNumber','if_dirac_comb','comb_onPeriod','comb_offPeriod', ...
    'if_limb_force','limb_force_strength','if_cell_division','area_0', ...
    'if_squeeze_tissue','squeeze_when','percent_squeeze', ...
    'if_active_contractility','tau_contr','active_contr_strength', ...
    'if_ABP','vo','Dr','if_RhoROCK','if_RK4','nhill','K_hill','A_Rho','A_ROCK', ...
    'A_Myosin','D_rho','D_ROCK','D_Myosin','Myosin_Coupling_Strength', ...
    'if_myosin_noise','if_gaussian_noise','Myosin_noise_strength', ...
    'if_coupling_noise','coupling_noise_strength', ...
    'if_polar_motility','polar_motility_strength' ...
    };

fid = fopen(filename, 'r');
if fid == -1
    error('ReadPara1Params:fileNotFound', 'Could not open %s', filename);
end

p = struct();
for i = 1:numel(names)
    line = fgetl(fid);
    if ~ischar(line)
        fclose(fid);
        error('ReadPara1Params:tooFewLines', ...
            ['para1_in.dat has fewer lines (%d) than expected parameters ' ...
             '(%d) -- the file format may have changed; check it against ' ...
             'allocation.f90::read_input''s read order and update the ' ...
             '''names'' list in this function to match.'], i-1, numel(names));
    end
    p.(names{i}) = parseToken(line);
end
fclose(fid);

end


function val = parseToken(line)
% First whitespace-delimited token on the line (comment-safe: para1_in.dat
% lines look like "1.0d-3           ! dt", so the first token is the value).
tok = strtrim(strtok(line));

lowtok = lower(tok);
if strcmp(lowtok, '.true.') || strcmp(lowtok, 't')
    val = true;
elseif strcmp(lowtok, '.false.') || strcmp(lowtok, 'f')
    val = false;
else
    % Fortran double-precision literals use 'd'/'D' for the exponent
    % (e.g. 1.0d-3) -- str2double doesn't understand that, so normalize
    % to 'e'/'E' first.
    numtok = strrep(strrep(tok, 'd', 'e'), 'D', 'E');
    val = str2double(numtok);
    if isnan(val)
        error('ReadPara1Params:parseError', ...
            'Could not parse "%s" as a number or logical (from line: "%s")', ...
            tok, line);
    end
end
end
