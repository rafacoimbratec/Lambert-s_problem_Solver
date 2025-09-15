function [C,S] = stumpff(z)
%STUMPFF  Compute Stumpff functions C(z), S(z) robustly.
% Handles z ~ 0 with series expansions.

    epsZ = 1e-8;

    if abs(z) < epsZ
        % Series about z=0
        z2 = z*z;
        C = 1/2 - z/24 + z2/720 - z2*z/40320;
        S = 1/6 - z/120 + z2/5040 - z2*z/362880;
        return;
    end

    if z > 0
        s = sqrt(z);
        C = (1 - cos(s)) / z;
        S = (s - sin(s)) / (s^3);
    else
        s = sqrt(-z);
        C = (1 - cosh(s)) / z;     % note: z<0
        S = (sinh(s) - s) / (s^3);
    end
end
