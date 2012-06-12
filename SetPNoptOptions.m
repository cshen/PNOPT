function options = SetPNoptOptions(varargin)
% SetProxNewtonOptions : Set options for PNopt
%
% options = SetPNoptOptions('param1',val1,'param2',val2,...) creates an options
%   structure in which the named parameters have the specified values. Default
%   values are used for parameters that are not specified.
%   
% SetPNoptOptions with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
%   SetPNoptOptions.m
%   $Revision: 0.1.2 $  $Date: 2012/06/12 $

  % Print out possible values of options.
  if nargin == 0 && nargout == 0
    fprintf('  Default parameters for PNopt.m:\n');
    fprintf('          CheckOpt : [          logical |         1 ] \n');
    fprintf('           Display : [ positive integer |         1 ] \n');
    fprintf('  LbfgsCorrections : [ positive integer |        50 ] \n');
    fprintf('  LineSearchMemory : [ positive integer |        10 ] \n');
    fprintf('  LineSearchMethod : [ string           |     Curvy ] \n');
    fprintf('       MaxFunEvals : [ positive integer |      5000 ] \n');
    fprintf('           MaxIter : [ positive integer |       500 ] \n');
    fprintf('            Method : [ string           |     Lbfgs ] \n');
    fprintf('        PrintEvery : [ positive integer |        10 ] \n');
    fprintf('    SmallFirstStep : [ logical          |         1 ] \n');
    fprintf('  SubproblemMethod : [ string           |     Tfocs ] \n');
    fprintf('         TfocsOpts : [ struct           | TfocsOpts ] \n');
    fprintf('            TolFun : [ positive scalar  |      1e-9 ] \n');
    fprintf('            TolOpt : [ positive scalar  |      1e-6 ] \n');
    fprintf('              TolX : [ positive scalar  |      1e-9 ] \n');
    fprintf('           UseHess : [ logical          |         0 ] \n');
    fprintf('            UseMex : [ logical          |         1 ]\n');
    fprintf('\n');
    return;
  end

  Names = [
    'SpgOptions       '
    'CheckOpt         '
    'Display          '
    'Hess             '
    'LbfgsCorrections '
    'LineSearchMemory '
    'LineSearchMethod '
    'MaxFunEvals      '
    'MaxIter          '
    'Method           '
    'PrintEvery       '
    'SmallFirstStep   '
    'SubproblemMethod '
    'TfocsOpts        '
    'TolFun           '
    'TolOpt           '
    'TolX             '
    'UseHess          '
    'UseMex           '
    ];
  [m,n] = size(Names);
  names = lower(Names);

  % Combine all leading options structures o1, o2, ... in l1Set(o1,o2,...).
  options = [];
  for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
  end
  i = 1;
  while i <= nargin
    arg = varargin{i};
    if ischar(arg), break; end
    if ~isempty(arg)                      % [] is a valid options argument
      if ~isa(arg,'struct')
        error(sprintf(['Expected argument %d to be a string parameter name ' ...
          'or an options structure\ncreated with SetPNoptOptions.'], i));
      end
      for j = 1:m
        if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
          eval(['val = arg.' Names(j,:) ';']);
        else
          val = [];
        end
        if ~isempty(val)
          eval(['options.' Names(j,:) '= val;']);
        end
      end
    end
    i = i + 1;
  end

  % A finite state machine to parse name-value pairs.
  if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
  end
  expectval = 0;                          % start expecting a name, not a value
  while i <= nargin
    arg = varargin{i};

    if ~expectval
      if ~ischar(arg)
         error(sprintf('Expected argument %d to be a string parameter name.', i));
      end

      lowArg = lower(arg);
      j = strmatch(lowArg,names);
      if isempty(j)                       % if no matches
        error(sprintf('Unrecognized parameter name ''%s''.', arg));
      elseif length(j) > 1                % if more than one match
        % Check for any exact matches (in case any names are subsets of others)
        k = strmatch(lowArg,names,'exact');
        if length(k) == 1
          j = k;
        else
          msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
          msg = [msg '(' deblank(Names(j(1),:))];
          for k = j(2:length(j))'
            msg = [msg ', ' deblank(Names(k,:))];
          end
          msg = sprintf('%s).', msg);
          error(msg);
        end
      end
      expectval = 1;                      % we expect a value next

    else
      eval(['options.' Names(j,:) '= arg;']);
      expectval = 0;

    end
    i = i + 1;
  end

  if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
  end
end
