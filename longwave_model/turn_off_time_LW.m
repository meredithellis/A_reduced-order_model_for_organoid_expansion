function [rW,f] = turn_off_time_LW(t, W, W_maxtol)

if nargin<7
    W_maxtol=0.7;
end
j=1;
%W_maxtol=0.7;


for i=1:length(t)

% r =  roots(W(:,i)-W_maxtol);
%  
%   if isempty(r)==0   % if there is a root of W-0.7=0
%   %
%   elseif max(W(:,i))> W_maxtol % if there is no root, but the max(W)>Wmaxtol, then have W>Wmaxtol
% %
%   else % no root and max W<Wmaxtol
% 
%     
%    
%  end
r =  roots(W(:,i)-W_maxtol);
 if isempty(r)~=0 && max(W(:,i)) < W_maxtol
            j=j+1;
    rW(j,1) = t(i);
    rW(j,2) = 1;
 elseif isempty(r)==0 && max(W(:,i))< W_maxtol
                 j=j+1;
    rW(j,1) = t(i);
    rW(j,2) = min(r);
 end
 
 
f = max(rW(j,1));






end
end