%UNIT_TEST runs unit tests for the Quaternion SIMULINK Library.

% $Source: /home/stpierre/cvsroot/matlab/simulink/qlib/test/unit_test.m,v $
% $Revision: 1.5 $
% $Date: 2009-12-12 20:13:00 $

% Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.

if ~exist('check_value','file')
  error(['You do not appear to have the test_tools toolbox installed,', 10,...
         'which is necessary to run these unit tests.  You should be', 10,...
         'able to get the test_tools toolbox from the same place you', 10,...
         'obtained this toolbox.']);
end

unix('rm -f test.out')
diary test.out

test_title = 'Quaternion SIMULINK Library blocks';
disp_test_title(test_title);

%%% Quaternion

q=[1 2 3 4];
sim_q=q;

%%% Vector

v=[1 2 3];
sim_v=v;

%%% Run first test model

qtest1
sim('qtest1')

%%% Check results

failures=0;

disp_test_name('Quaternion Normalize');
truth_value = 'qnorm(sim_q)';
test_value  = 'sim_qnorm';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Quaternion Conjugate');
truth_value = 'qconj(sim_q)';
test_value  = 'sim_qconj';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Quaternion Multiply');
truth_value = 'qmult(sim_q,sim_q)';
test_value  = 'sim_qmult';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Quaternion Decomposition');
[v,phi]=qdecomp(qnorm(sim_q)) %#ok<NOPTS>
truth_value = '[v phi]';
test_value  = '[sim_vector sim_phi]';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Quaternion Vector Transform');
truth_value = 'qvxform(qnorm(sim_q),sim_v)';
test_value  = 'sim_qvxform';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Quaternion Vector Rotation');
truth_value = 'qvrot(qnorm(sim_q),sim_v)';
test_value  = 'sim_qvrot';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Quaternion to DCM');
dcm='q2dcm(qnorm(sim_q))',disp(' ') %#ok<NOPTS>
dcm=eval(dcm);
truth_value = '[dcm(1,:) dcm(2,:) dcm(3,:)]';
test_value  = 'sim_q2dcm';
failures=failures+check_value(truth_value, test_value);

disp_test_name('DCM to Quaternion');
truth_value = 'dcm2q(dcm)';
test_value  = 'sim_dcm2q';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Row Major to Matrix');
truth_value = 'reshape(sim_q2dcm,3,3)''';
test_value  = 'sim_rm2mat';
failures=failures+check_value(truth_value, test_value);

disp_test_name('Matrix to Row Major');
truth_value = 'sim_q2dcm';
test_value  = 'sim_mat2rm';
failures=failures+check_value(truth_value, test_value);

%%% Run second test model

disp_test_name('Quaternion Integration');

step_size=0.01, disp(' ') %#ok<NOPTS>
stop_time=10, disp(' ') %#ok<NOPTS>
decimation=100, disp(' ') %#ok<NOPTS>
set_val('t', '0:step_size*decimation:stop_time')
set_val('zero_col', 'zeros(length(t), 1)')

qtest2
sim('qtest2')

%%% Check results

disp_test_name('Quaternion Integration: X Rotation');
truth_value = '[sin(t/2).'' zero_col zero_col cos(t/2).'']';
test_value  = 'xrot';
failures=failures+check_value(truth_value, test_value, 5e-15);

disp_test_name('Quaternion Integration: Y Rotation');
truth_value = '[zero_col sin(t/2).'' zero_col cos(t/2).'']';
test_value  = 'yrot';
failures=failures+check_value(truth_value, test_value, 5e-15);

disp_test_name('Quaternion Integration: Z Rotation');
truth_value = '[zero_col zero_col sin(t/2).'' cos(t/2).'']';
test_value  = 'zrot';
failures=failures+check_value(truth_value, test_value, 5e-15);

disp_num_failures(test_title, failures)

diary off
