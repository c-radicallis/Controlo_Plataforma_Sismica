function d2y = secondDerivativeTime3(y, dt)
% secondDerivativeTime3  Compute second derivative of a 1D signal in time
%                       using a five-point central difference interior.
%
%   d2y = secondDerivativeTime3(y, dt) returns a vector d2y of the same
%   length as y, where each element is an approximation of d^2y/dt^2 at
%   the corresponding time point. y must be a row or column vector; dt is
%   the constant time step between samples.
%
%   Interior (i = 3..N-2), fourth-order accurate:
%       d2y(i) ≈ [ -y(i+2) + 16*y(i+1) - 30*y(i) + 16*y(i-1) - y(i-2) ] 
%                 / (12 * dt^2)
%
%   Boundaries:
%     • i = 1: forward, second-order:
%         d2y(1)   ≈ ( y(3) - 2*y(2) +   y(1) ) / dt^2
%     • i = 2: central, second-order:
%         d2y(2)   ≈ ( y(3) - 2*y(2) +   y(1) ) / dt^2
%     • i = N-1: central, second-order:
%         d2y(N-1) ≈ ( y(N) - 2*y(N-1) + y(N-2) ) / dt^2
%     • i = N: backward, second-order:
%         d2y(N)   ≈ ( y(N) - 2*y(N-1) + y(N-2) ) / dt^2
%
%   Example:
%       t  = 0:0.01:2;                   
%       y  = sin(2*pi*5*t);             
%       dt = t(2) - t(1);                
%       d2y = secondDerivativeTime3(y, dt);
%       plot(t, d2y, 'b', t, -(2*pi*5)^2*sin(2*pi*5*t), 'r--');
%       legend('Numeric','Analytic','Location','Best');

    % Ensure column for indexing
    origIsRow = isrow(y);
    y = y(:);
    N = numel(y);
    if N < 5
        error('Need at least 5 samples for a five-point stencil.');
    end

    d2y = zeros(N,1);

    %--- boundaries (second-order) ---
    % i = 1: forward
    d2y(1)   = ( y(3)   - 2*y(2)   +   y(1) ) / dt^2;
    % i = 2: simple central
    d2y(2)   = ( y(3)   - 2*y(2)   +   y(1) ) / dt^2;

    %--- interior (fourth-order five-point) ---
    for i = 3:(N-2)
        d2y(i) = ( -y(i+2) ...
                   +16*y(i+1) ...
                   -30*y(i) ...
                   +16*y(i-1) ...
                   -y(i-2) ) ...
                 / (12 * dt^2);
    end

    %--- boundaries (second-order) ---
    % i = N-1: central
    d2y(N-1) = ( y(N)   - 2*y(N-1) + y(N-2) ) / dt^2;
    % i = N: backward
    d2y(N)   = ( y(N)   - 2*y(N-1) + y(N-2) ) / dt^2;

    % restore shape
    if origIsRow
        d2y = d2y.';
    end
end
