function C = checkerboard_2(n)
  % case of even rows
  if mod(n, 2) == 0
    C = zeros(n + 1, n);    % create a matrix having a plus row
    C(2 : 2 : end) = 1;     % select every second element
    C(end, :) = [];         % remove the last row from the matrix
  else
    C = zeros(n, n);        % create the matrix
    C(2 : 2 : end) = 1;     % select every second element
  end
end