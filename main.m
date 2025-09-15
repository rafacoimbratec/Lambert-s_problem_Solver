%% main.m — Lambert comparison driver (extensible)
clear; clc;

%% 1) Problem definition (consistent units)
mu = 398600.4418;                % [km^3/s^2] Earth
r1 = [  5000; 10000; 2100 ];     % [km]
r2 = [ -14600;  2500; 7000 ];    % [km]
dt = 1*60*60;                    % [s] 1 hour

%% 2) Methods to run (add more later)
% Each cell is {name, function_handle, name_value_args{:}}
methods = {
    'UV', @lambert_uv, 'z0', 1.5, 'verbose', true   % Universal Variables Iteractive method
    % 'Battin', @lambert_battin, ...                 % (when you add it)
    % 'Izzo',   @lambert_izzo,   ...
    % 'Gooding',@lambert_gooding,...
};

%% 3) Run all methods, collect results
results = struct('method',{},'v1',{},'v2',{},'info',{},'coe',{});
for k = 1:size(methods,1)
    name = methods{k,1};
    fun  = methods{k,2};
    args = methods(k,3:end);

    try
        % Call solver
        [v1, v2, info] = fun(r1, r2, dt, mu, args{:});

        % Compute COE from (r1, v1)
        coe = coe_from_sv(r1, v1, mu);   % [h e RAAN i w TA a]

        % Store
        results(end+1) = struct( ...
            'method', name, ...
            'v1', v1, ...
            'v2', v2, ...
            'info', info, ...
            'coe', coe);
    catch ME
        warning('Method "%s" failed: %s', name, ME.message);
        results(end+1) = struct( ...
            'method', name, 'v1', nan(3,1), 'v2', nan(3,1), ...
            'info', struct('converged',false,'message','ERROR'), ...
            'coe', nan(1,7));
    end
end

%% 4) Pretty print per-method report
fprintf('\n================ LAMBERT SOLVER REPORT ================\n');
fprintf('r1 = [%10.3f %10.3f %10.3f] km\n', r1);
fprintf('r2 = [%10.3f %10.3f %10.3f] km\n', r2);
fprintf('Δt = %g s (%.3f h)   μ = %.6f km^3/s^2\n\n', dt, dt/3600, mu);

for k = 1:numel(results)
    R = results(k);
    fprintf('--- Method: %s ---\n', R.method);

    % Convergence info
    if isfield(R,'info') && isfield(R.info,'converged')
        fprintf('  converged: %s   iters: %d   |F|: %.3e s   z: %.6g\n', ...
            tf(R.info.converged), getfielddef(R.info,'iterations',NaN), ...
            getfielddef(R.info,'tof_err',NaN), getfielddef(R.info,'z',NaN));
        if isfield(R.info,'message')
            fprintf('  note: %s\n', R.info.message);
        end
    end

    % Velocities
    fprintf('  v1 (km/s): [%12.6f %12.6f %12.6f]\n', R.v1);
    fprintf('  v2 (km/s): [%12.6f %12.6f %12.6f]\n', R.v2);

    % COE table (angles in deg)
    if all(isfinite(R.coe))
        coe = R.coe;
        h     = coe(1);
        e     = coe(2);
        RAAN  = rad2deg(coe(3));
        inc   = rad2deg(coe(4));
        argp  = rad2deg(coe(5));
        nu    = rad2deg(coe(6));
        a     = coe(7);

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
    fprintf('\n');
end

%% 5) (Optional) Cross-check by propagating r1,v1 to t=dt
% Pick first successful solution
%ok = find(arrayfun(@(S)isfield(S,'info') && isfield(S.info,'converged') && S.info.converged, results), 1, 'first');
%if ~isempty(ok)
%    v1_ok = results(ok).v1;
    % If you have kepler_universal.m, uncomment to validate:
    % [r2p,~] = kepler_universal(r1, v1_ok, dt, mu);
    % fprintf('Endpoint error |r2p-r2| = %.6e km\n', norm(r2p - r2));
%end

%% 6) (Optional) Plotting
% plots.plot_transfer(...) or your own plot once you write it.

% ----------------- small helpers -----------------

