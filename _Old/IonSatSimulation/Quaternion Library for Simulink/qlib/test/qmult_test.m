%QMULT_TEST runs timing tests for Quaternion Multiplication algorithms.

% $Source: /home/stpierre/cvsroot/matlab/simulink/qlib/test/qmult_test.m,v $
% $Revision: 1.3 $
% $Date: 2009-07-21 03:23:26 $

% Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.

if ~exist('check_value','file')
  error(['You do not appear to have the test_tools toolbox installed,', 10,...
         'which is necessary to run this test.  You should be able to', 10,...
         'get the test_tools toolbox from the same place you obtained', 10,...
         'this toolbox.']);
end

stop_time=12e5;
q1=[ 1 2 3  4];
q2=[-1 2 1 -3];

save qmult_test_vars stop_time q1 q2

qmult_test1
qmult_test2
disp(' ')

clear
disp('qmult_test1:')
disp('============')
load qmult_test_vars
ts=timesim('qmult_test1');
truth_value=('qmult(q1,q2)');
test_value=('simout');
check_value(truth_value, test_value);
disp(['Time: ', num2str(ts)])
disp(' ')

clear
disp('qmult_test2:')
disp('============')
load qmult_test_vars
ts=timesim('qmult_test2');
truth_value=('qmult(q1,q2)');
test_value=('simout');
check_value(truth_value, test_value);
disp(['Time: ', num2str(ts)])
disp(' ')
