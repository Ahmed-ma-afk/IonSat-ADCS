%QNORM_TEST runs timing tests for Quaternion Normalization algorithms.

% $Source: /home/stpierre/cvsroot/matlab/simulink/qlib/test/qnorm_test.m,v $
% $Revision: 1.3 $
% $Date: 2009-07-21 03:23:26 $

% Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.

if ~exist('timesim','file')
  error(['You do not appear to have the test_tools toolbox installed,', 10,...
         'which is necessary to run this test.  You should be able to', 10,...
         'get the test_tools toolbox from the same place you obtained', 10,...
         'this toolbox.']);
end

stop_time=12e5;
q=[1 2 3 4];

save qnorm_test_vars stop_time q

qnorm_test1
qnorm_test2
qnorm_test3
disp(' ')

clear
disp('qnorm_test1:')
disp('========')
load qnorm_test_vars
ts=timesim('qnorm_test1');
if simout == qnorm(q)
  disp('simout == qnorm(q)  ***PASSED***')
else
  disp('simout ~= qnorm(q)  ***FAILED***')
  simout
  qnorm(q)
end
disp(['Time: ', num2str(ts)])
disp(' ')

clear
disp('qnorm_test2:')
disp('========')
load qnorm_test_vars
ts=timesim('qnorm_test2');
if simout == qnorm(q)
  disp('simout == qnorm(q)  ***PASSED***')
else
  disp('simout ~= qnorm(q)  ***FAILED***')
  simout
  qnorm(q)
end
disp(['Time: ', num2str(ts)])
disp(' ')

clear
disp('qnorm_test3:')
disp('========')
load qnorm_test_vars
ts=timesim('qnorm_test3');
if simout == qnorm(q)
  disp('simout == qnorm(q)  ***PASSED***')
else
  disp('simout ~= qnorm(q)  ***FAILED***')
  simout
  qnorm(q)
end
disp(['Time: ', num2str(ts)])
disp(' ')

