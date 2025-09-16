%% main.m — Lambert + Propagation Comparison Report
clear; clc;

%% 1) Problem definition (consistent units)
mu = 398600.4418;                % [km^3/s^2] Earth
r1 = [  5000; 10000; 2100 ];     % [km]
r2 = [ -14600;  2500; 7000 ];    % [km]
dt = 1*60*60;                    % [s] 1 hour

%% 2) Methods to run (add more rows as you implement others)
% Each row: {name, function_handle, name_value_args{:}}
methods = {
    'UV', @lambert_uv, 'z0', 1.5, 'verbose', true
    % 'Battin',  @lambert_battin, ...
    % 'Gooding', @lambert_gooding, ...
};

%% 3) Run all methods, compute COE, propagate, collect metrics
results = struct('method',{},'v1',{},'v2',{},'info',{},'coe',{}, ...
                 'r2p',{},'v2p',{},'pos_err_km',{},'vel_err_kms',{});

for k = 1:size(methods,1)
    name = methods{k,1};
    fun  = methods{k,2};
    args = methods(k,3:end);

    try
        % --- Solve Lambert ---
        [v1, v2, info] = fun(r1, r2, dt, mu, args{:});

        % --- COEs from (r1,v1) ---
        coe = coe_from_sv(r1, v1, mu);   % [h e RAAN i w TA a]

        % --- Propagate (r1,v1) to t=dt and compare to (r2,v2) ---
        [r2p, v2p] = kepler_universal(r1, v1, dt, mu);  %#ok<ASGLU> (keep v2p too)
        pos_err = norm(r2p - r2);       % [km]
        vel_err = norm(v2p - v2);       % [km/s]

        % --- Store ---
        results(end+1) = struct( ...
            'method', name, ...
            'v1', v1, ...
            'v2', v2, ...
            'info', info, ...
            'coe', coe, ...
            'r2p', r2p, ...
            'v2p', v2p, ...
            'pos_err_km', pos_err, ...
            'vel_err_kms', vel_err);
    catch ME
        warning('Method "%s" failed: %s', name, ME.message);
        results(end+1) = struct( ...
            'method', name, 'v1', nan(3,1), 'v2', nan(3,1), ...
            'info', struct('converged',false,'iterations',NaN,'tof_err',NaN,'z',NaN,'message','ERROR'), ...
            'coe', nan(1,7), ...
            'r2p', nan(3,1), 'v2p', nan(3,1), ...
            'pos_err_km', NaN, 'vel_err_kms', NaN);
    end
end

%% 4) Summary header
fprintf('\n================ LAMBERT + PROPAGATION REPORT ================\n');
fprintf('r1 = [%10.3f %10.3f %10.3f] km\n', r1);
fprintf('r2 = [%10.3f %10.3f %10.3f] km\n', r2);
fprintf('Δt = %g s (%.3f h)   μ = %.6f km^3/s^2\n\n', dt, dt/3600, mu);

%% 5) Compact summary table (one line per method)
% Columns: Method | Conv | Iters | |F| [s] | z | |r2p-r2| [km] | |v2p-v2| [km/s]
fmtHeader = '%-10s  %-5s  %5s  %12s  %12s  %14s  %12s\n';
fmtRow    = '%-10s  %-5s  %5d  %12.3e  %12.5g  %14.6e  %12.6e\n';
fprintf(fmtHeader, 'Method','Conv','Iters','|F| [s]','z','|Δr| [km]','|Δv| [km/s]');
for k = 1:numel(results)
    R = results(k);
    fprintf(fmtRow, ...
        R.method, tf(getfielddef(R.info,'converged',false)), ...
        getfielddef(R.info,'iterations',NaN), ...
        getfielddef(R.info,'tof_err',NaN), ...
        getfielddef(R.info,'z',NaN), ...
        R.pos_err_km, R.vel_err_kms);
end
fprintf('\n');

%% 6) Detailed per-method section
for k = 1:numel(results)
    R = results(k);
    fprintf('--- %s details ---\n', R.method);

    % Convergence info
    fprintf('  converged: %s   iters: %d   |F|: %.3e s   z: %.6g\n', ...
        tf(getfielddef(R.info,'converged',false)), ...
        getfielddef(R.info,'iterations',NaN), ...
        getfielddef(R.info,'tof_err',NaN), ...
        getfielddef(R.info,'z',NaN));
    if isfield(R.info,'message')
        fprintf('  note: %s\n', R.info.message);
    end

    % Velocities
    fprintf('  v1 (km/s): [%12.6f %12.6f %12.6f]\n', R.v1);
    fprintf('  v2 (km/s): [%12.6f %12.6f %12.6f]\n', R.v2);

    % COE (angles in degrees)
    if all(isfinite(R.coe))
        h     = R.coe(1);
        e     = R.coe(2);
        RAAN  = rad2deg(R.coe(3));
        inc   = rad2deg(R.coe(4));
        argp  = rad2deg(R.coe(5));
        nu    = rad2deg(R.coe(6));
        a     = R.coe(7);

        fprintf('  COE:\n');
        fprintf('    h      = % .6f km^2/s\n', h);
        fprintf('    e      = % .8f [-]\n',   e);
        fprintf('    RAAN   = % .6f deg\n',   RAAN);
        fprintf('    i      = % .6f deg\n',   inc);
        fprintf('    ω      = % .6f deg\n',   argp);
        fprintf('    ν      = % .6f deg\n',   nu);
        fprintf('    a      = % .6f km\n',    a);
    else
        fprintf('  COE: n/a (solver failed)\n');
    end

    % Propagation check
    fprintf('  Propagate (r1,v1) → t=Δt   |Δr| = %.6e km,  |Δv| = %.6e km/s\n\n', ...
        R.pos_err_km, R.vel_err_kms);
end

%% 7) (Optional) Plot a successful transfer
% Pick first converged method to plot (once you have a plotter)
%ok = find(arrayfun(@(S)isfield(S,'info') && S.info.converged, results), 1, 'first');
%if ~isempty(ok)
    %plot_transfer3d(r1, r2, results(ok).v1, mu, dt, ...
        %'N', 400, 'BodyRadius', 6378, 'BodyName', 'Earth', ...
       % 'ShowVel', true, 'ShowPlane', true, 'NewFigure', true);
%end

%% ----------------- helpers (local) -----------------
function s = tf(b)
% true/false → 'yes'/'no' for neat printing
    if islogical(b) && isscalar(b)
        s = ternary(b,'yes','no');
    else
        s = 'n/a';
    end
end

function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end

function v = getfielddef(S, fld, default)
% get field with default if missing
    if isstruct(S) && isfield(S, fld), v = S.(fld); else, v = default; end
end

