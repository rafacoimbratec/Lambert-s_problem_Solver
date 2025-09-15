function [v1, v2, info] = lambert_uv(r1, r2, dt, mu, varargin)
%LAMBERT_UV  Universal-Variables iterative Lambert solver (0-rev).
%
%   [v1,v2,info] = lambert_uv(r1, r2, dt, mu, 'Name',Value, ...)
%
% Inputs
%   r1,r2 : 3x1 position vectors [km]
%   dt    : time of flight [s], dt>0
%   mu    : gravitational parameter [km^3/s^2]
%
% Name–Value options
%   'longway'  (false)  -> if true, solve the long-way (Δθ in (π,2π))
%   'tol'      (1e-10)  -> tolerance on TOF residual |T(z)-dt| [s]
%   'maxIter'  (50)     -> maximum iterations
%   'verbose'  (false)  -> print iteration trace
%   'zMin'     (-4)     -> lower bound for universal variable z (hyperbolic)
%   'zMax'     ( +4)    -> upper bound for z (elliptic)
%
% Outputs
%   v1,v2 : 3x1 velocities at r1 and r2 [km/s]. NaNs if not converged.
%   info  : struct with fields:
%       .converged, .iterations, .z, .tof_err, .theta, .A, .y
%       .f, .g, .gdot
%       .message
%
% Notes
%   - Implements the classic UV Lambert with iteration on z.
%   - Targets 0-revolution solutions. (Multi-rev requires separate bracketing.)
%   - Requires helper STUMPFF.M on MATLAB path.

    %----------- parse & validate inputs (robust for older MATLAB) ----------
    ip = inputParser;
    ip.FunctionName = 'lambert_uv';
    addRequired(ip,'r1',@(x)isnumeric(x)&&isvector(x)&&numel(x)==3);
    addRequired(ip,'r2',@(x)isnumeric(x)&&isvector(x)&&numel(x)==3);
    addRequired(ip,'dt',@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addRequired(ip,'mu',@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'longway',false,@(x)islogical(x)||ismember(x,[0,1]));
    addParameter(ip,'tol',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'maxIter',50,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
    addParameter(ip,'verbose',false,@(x)islogical(x)||ismember(x,[0,1]));
    addParameter(ip,'zMin',-4,@(x)isnumeric(x)&&isscalar(x)&&x<0);
    addParameter(ip,'zMax',+4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(ip,r1,r2,dt,mu,varargin{:});
    opt = ip.Results;

    r1 = r1(:); r2 = r2(:);
    r1n = norm(r1); r2n = norm(r2);

    if r1n==0 || r2n==0
        error('lambert_uv:badInput','r1 and r2 must be nonzero.');
    end
    if opt.zMin >= opt.zMax
        error('lambert_uv:badBounds','zMin must be < zMax.');
    end

    %------------------------------ geometry --------------------------------
    % Transfer angle Δθ using atan2 with cross product
    cr = cross(r1,r2);
    sin_d = norm(cr) / (r1n*r2n);                 % >=0
    cos_d = max(-1,min(1, dot(r1,r2)/(r1n*r2n))); % clamp for safety
    theta = atan2(sin_d, cos_d);                  % (0, π]
    if opt.longway
        theta = 2*pi - theta;                     % (π, 2π)
    end

    den = 1 - cos(theta);
    if den <= eps
        % Δθ ~ 0 (r1≈r2 and aligned). UV Lambert becomes ill-conditioned.
        v1 = nan(3,1); v2 = nan(3,1);
        info = packInfo(false,0,NaN,Inf,theta,NaN,NaN,NaN,NaN,NaN,'Δθ≈0: degenerate geometry');
        return;
    end

    % A parameter (captures short/long via sin(theta) sign)
    A = sin(theta) * sqrt( (r1n*r2n) / den );
    if abs(A) < 1e-14
        v1 = nan(3,1); v2 = nan(3,1);
        info = packInfo(false,0,NaN,Inf,theta,A,NaN,NaN,NaN,NaN,'A≈0: degenerate geometry');
        return;
    end

    %------------------------ residual function F(z) ------------------------
    % F(z) = T(z) - dt, where T(z) is UV time-of-flight
    function [F, y, C, S] = tof_residual(z)
        [C, S] = stumpff(z);

        % Guard for tiny/negative numerical C
        if (z > 0) && (C <= 0)
            C = max(C, realmin);
        end

        sqrtC = sqrt(C);
        % y(z) formula (Vallado/Lancaster-Blanchard UV Lambert)
        y = r1n + r2n + A * (z*S - 1) / sqrtC;

        if ~(isfinite(y)) || y <= 0
            F = +Inf;   % invalid branch for this z
            return;
        end
        chi = sqrt(y/C);
        Tz  = chi^3 * S + A * sqrt(y);   % note: units sqrt(mu)*seconds
        F   = Tz/sqrt(mu) - dt;          % seconds
    end

    %------------------------- initial guess for z --------------------------
    z = 0.0;                              % parabolic-like start
    [F, y] = tof_residual(z);
    if ~isfinite(F)
        % Probe a few points within bounds to get a finite residual
        probes = linspace(opt.zMin,opt.zMax,9);
        got = false;
        for p = probes
            [Fp, yp] = tof_residual(p);
            if isfinite(Fp)
                z = p; F = Fp; y = yp; got = true; break;
            end
        end
        if ~got
            v1 = nan(3,1); v2 = nan(3,1);
            info = packInfo(false,0,NaN,Inf,theta,A,NaN,NaN,NaN,NaN,'Failed to find valid initial z');
            return;
        end
    end

    %------------------ iterate (safeguarded Newton/secant) -----------------
    iter = 0;
    converged = false;
    last_z = z; last_F = F;

    while iter < opt.maxIter
        iter = iter + 1;

        % numerical derivative dF/dz
        dz = max(1e-8, 1e-6*(1+abs(z)));
        [F2, ~] = tof_residual(z + dz);
        if ~isfinite(F2)
            [F2, ~] = tof_residual(z - dz);
            dF = (F - F2)/dz;
        else
            dF = (F2 - F)/dz;
        end

        if ~isfinite(dF) || dF == 0
            % fallback to secant using last step if possible
            if iter > 1 && (z ~= last_z) && isfinite(last_F)
                dF = (F - last_F) / (z - last_z);
                if ~isfinite(dF) || dF == 0
                    dF = sign(F)*max(1,abs(F)); % crude but safe
                end
            else
                dF = sign(F)*max(1,abs(F));
            end
        end

        % Newton step
        z_new = z - F/dF;

        % Safeguard: keep within [zMin,zMax] and avoid huge jumps
        if z_new < opt.zMin, z_new = 0.5*(z + opt.zMin); end
        if z_new > opt.zMax, z_new = 0.5*(z + opt.zMax); end
        if ~isfinite(z_new)
            z_new = 0.5*(opt.zMin + opt.zMax);
        end

        % evaluate at new z
        last_z = z; last_F = F;
        z = z_new;
        [F, y] = tof_residual(z);

        if opt.verbose
            fprintf('it=%2d  z=% .6e  F=% .3e\n', iter, z, F);
        end

        if isfinite(F) && abs(F) < opt.tol
            converged = true;
            break;
        end
    end

    if ~converged || ~isfinite(F) || y <= 0
        v1 = nan(3,1); v2 = nan(3,1);
        info = packInfo(false,iter,z,abs(F),theta,A,y,NaN,NaN,NaN,'Did not converge or invalid y');
        return;
    end

    %-------------------- recover Lagrange f,g,gdot and v -------------------
    f    = 1 - y/r1n;
    g    = A * sqrt(y/mu);
    gdot = 1 - y/r2n;

    if abs(g) < 1e-14 || ~isfinite(g)
        v1 = nan(3,1); v2 = nan(3,1);
        info = packInfo(false,iter,z,abs(F),theta,A,y,f,g,gdot,'g≈0 (singular)');
        return;
    end

    v1 = (r2 - f*r1) / g;
    v2 = (gdot*r2 - r1) / g;

    info = packInfo(true,iter,z,abs(F),theta,A,y,f,g,gdot,'OK');
end

%============================= helpers =====================================
function s = packInfo(conv,it,z,toferr,theta,A,y,f,g,gdot,msg)
    s = struct('converged',logical(conv), ...
               'iterations',it, ...
               'z',z, ...
               'tof_err',toferr, ...
               'theta',theta, ...
               'A',A, ...
               'y',y, ...
               'f',f, ...
               'g',g, ...
               'gdot',gdot, ...
               'message',msg);
end

