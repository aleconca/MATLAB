%9.1
%a
weights_1 = 1/4*[1 2 1];
nodes_1   = [0 1/2 1];
f = @(x) x.^0; sum(weights_1 .* f(nodes_1))
f = @(x) x.^1; sum(weights_1 .* f(nodes_1))% inexact

%b
%yes, midpoint



%9.2
x = [0 0 0];
nodes_1   = [0 1 0];

f = @(x) x.^0; sum(x.* f(nodes_1))=1
f = @(x) x.^1; sum(x .* f(nodes_1))=0.5



%9.3
