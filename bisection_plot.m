function bisection_plot(f, a, b, tol)
  %BISECTION_PLOT Plot the function f and the evolution of the bisection method
  %
  %   BISECTION_PLOT(f, a, b, x_iter)
  %
  %   Inputs  : f     = function handle to the function f(x)
  %             a     = left bound
  %             b     = right bound
  %             tol   = requested tolerance

  figure
  hold on, box on
  x_plot = linspace(a, b, 1000);
  plot(x_plot, f(x_plot), 'LineWidth', 2)
  plot(x_plot, 0*x_plot, 'k-', 'LineWidth', 1)
  xlabel('x','FontSize', 16)
  ylabel('f(x)','FontSize', 16)
  set(gca,'FontSize', 16)
  set(gca,'LineWidth', 1.5)
  pause

  left_color = 'rx-';
  right_color = 'cx-';
  old_color = 'kx-';
  new_color = 'gx-';
  plot([a a], [0 f(a)], left_color, 'LineWidth', 2, 'MarkerSize', 8)
  plot([b b], [0 f(b)], right_color, 'LineWidth', 2, 'MarkerSize', 8)
  pause

  % Number of necessary iterations to match the requested tolerance
  max_iterb = ceil( log((b - a)/tol)/log(2) );
  
  for (iter = 1:max_iterb)
    x_iter(iter) = a + (b-a)/2;
    f_iter(iter) = f(x_iter(iter));
    
    if ( f(b)*f_iter(iter) < 0 )
      plot([x_iter(iter) x_iter(iter)], [0 f_iter(iter)], new_color, 'LineWidth', 2, 'MarkerSize', 8)
      pause
      plot([a a], [0 f(a)], old_color, 'LineWidth', 2, 'MarkerSize', 8)
      a = x_iter(iter);
      plot([x_iter(iter) x_iter(iter)], [0 f_iter(iter)], left_color, 'LineWidth', 2, 'MarkerSize', 8)
      pause
    elseif ( f(a)*f_iter(iter) < 0 )
      plot([x_iter(iter) x_iter(iter)], [0 f_iter(iter)], new_color, 'LineWidth', 2, 'MarkerSize', 8)
      pause
      plot([b b], [0 f(b)], old_color, 'LineWidth', 2, 'MarkerSize', 8)
      b = x_iter(iter);
      plot([x_iter(iter) x_iter(iter)], [0 f_iter(iter)], right_color, 'LineWidth', 2, 'MarkerSize', 8)
      pause
    else % f(x) = 0
      break;
    end
  end

end
