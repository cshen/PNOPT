function s = mergestruct(s1,s2)
% mergestruct : Merge two structures
% 
%   $Revision: 0.1.0 $  $Date: 2012/07/15 $
% 
  s = s1;
  
  names = fieldnames(s2);
  for k = 1:length(names)
    if isfield(s,names{k})
      if isstruct(s.(names{k}))
        s.(names{k}) = mergestruct(s.(names{k}),s2.(names{k}));
      else
        s.(names{k}) = s2.(names{k});
      end
    else
      s.(names{k}) = s2.(names{k});
    end
  end