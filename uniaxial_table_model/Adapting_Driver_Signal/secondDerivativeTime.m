function d2y = secondDerivativeTime(y, dt)
% secondDerivativeTime  Compute second derivative of a 1D signal in time.
%
%   d2y = secondDerivativeTime(y, dt) returns a vector d2y of the same
%   length as y, where each element is an approximation of d^2y/dt^2 at the
%   corresponding time point.  y must be a row or column vector, and dt is
%   the constant time step between samples.
%
%   Uses a second‐order central‐difference for interior points:
%       d2y(i) ≈ (y(i+1) – 2*y(i) + y(i-1)) / dt^2
%   For the first and last samples, it falls back to a one‐sided (forward/backward)
%   difference of second order:
%       d2y(1)   ≈ (y(3) – 2*y(2) + y(1)) / dt^2
%       d2y(end) ≈ (y(end) – 2*y(end-1) + y(end-2)) / dt^2
%
%   Example:
%       t  = 0:0.01:2;                   % time vector 0 to 2 s in 0.01 s steps
%       y  = sin(2*pi*5*t);             % a 5 Hz sine wave
%       dt = t(2) - t(1);                % sampling interval = 0.01 s
%       d2y = secondDerivativeTime(y, dt);
%
%   Plot to verify (analytically, d2(sin(2π5t))/dt^2 = - (2π5)^2 * sin(2π5t) ):
%       figure;
%       plot(t, d2y, 'b', t, -(2*pi*5)^2*sin(2*pi*5*t), 'r--');
%       legend('Numeric 2nd‐derivative','Analytic 2nd‐derivative','Location','Best');
%       xlabel('Time (s)'); ylabel('d^2y/dt^2');

    % Ensure y is a column vector for indexing convenience
    y = y(:);
    N = numel(y);
    if N < 3
        error('Input signal must have at least 3 samples.');
    end

    % Preallocate output
    d2y = zeros(N,1);

    % Forward‐difference (second order) for the first point:
    %   d2y(1) ≈ (y(3) - 2*y(2) + y(1)) / dt^2
    d2y(1) = (y(3) - 2*y(2) + y(1)) / dt^2;

    % Central‐difference (second order) for interior points 2..N-1:
    %   d2y(i) = (y(i+1) - 2*y(i) + y(i-1)) / dt^2
    for i = 2:(N-1)
        d2y(i) = (y(i+1) - 2*y(i) + y(i-1)) / dt^2;
    end

    % Backward‐difference (second order) for the last point:
    %   d2y(N) ≈ (y(N) - 2*y(N-1) + y(N-2)) / dt^2
    d2y(N) = (y(N) - 2*y(N-1) + y(N-2)) / dt^2;

    % If the original input was a row vector, return a row vector
    if isrow(y)
        d2y = d2y.';
    end
end
