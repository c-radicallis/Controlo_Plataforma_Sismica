function d2y = secondDerivativeTime7(y, dt)
% secondDerivativeTime7  Compute second derivative of a 1D signal in time
%                       using a seven‐point central‐difference stencil.
%
%   d2y = secondDerivativeTime7(y, dt) returns a vector d2y of the same
%   length as y, where each element is an approximation of d^2y/dt^2.
%   y must be a row or column vector; dt is the constant timestep.
%
%   Interior (i = 4..N-3), sixth-order accurate:
%     d2y(i) ≈ [  2*y(i-3) - 27*y(i-2) + 270*y(i-1)
%               - 490*y(i) + 270*y(i+1) - 27*y(i+2)
%               +   2*y(i+3) ] / (180 * dt^2)
%
%   Boundaries (i = 1,2,3 and N-2,N-1,N): fallback to three-point,
%   second-order:
%     d2y(i) ≈ ( y(i+1) - 2*y(i) + y(i-1) ) / dt^2
%   (for i=1 and i=N you use the one‐sided forward/backward versions)

    % reshape and check
    origIsRow = isrow(y);
    y = y(:);
    N = numel(y);
    if N < 7
        error('Need at least 7 samples for a seven-point stencil.');
    end

    d2y = zeros(N,1);

    %--- second-order boundaries ---
    % i = 1 forward
    d2y(1)   = (y(3) - 2*y(2) + y(1)) / dt^2;
    % i = 2,3 and i = N-2,N-1 central
    for i = 2:3
        d2y(i)     = (y(i+1) - 2*y(i) + y(i-1)) / dt^2;
        d2y(N+1-i) = (y(N+2-i) - 2*y(N+1-i) + y(N-i)) / dt^2;
    end
    % i = N backward
    d2y(N)   = (y(N) - 2*y(N-1) + y(N-2)) / dt^2;

    %--- interior seventh‐point, sixth-order stencil ---
    for i = 4:(N-3)
        d2y(i) = (  2*y(i-3) ...
                  -27*y(i-2) ...
                  +270*y(i-1) ...
                  -490*y(i)   ...
                  +270*y(i+1) ...
                  -27*y(i+2) ...
                  +  2*y(i+3) ) ...
                / (180 * dt^2);
    end

    % restore shape
    if origIsRow
        d2y = d2y.';
    end
end
