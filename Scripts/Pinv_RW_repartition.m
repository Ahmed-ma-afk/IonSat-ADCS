%Original by X19 students
% alpha = pi / 6;
% sat.wheel.repartition_matrix_4RW = [[sin(alpha), sin(alpha), sin(alpha), sin(alpha)],
%     [cos(alpha), 0, -cos(alpha), 0],
%     [0, cos(alpha), 0, -cos(alpha)]];
% sat.wheel.repartition_metrix_4RW_inverse = pinv(sat.wheel.repartition_matrix_4RW);

%Up to date with CubeSpace RWA
alpha = 26.6 * pi/180;      %from deg to rad
sat.wheel.repartition_matrix_4RW = [[cos(alpha),      0,    -cos(alpha),      0    ];...
                                    [     0,     cos(alpha),     0,     -cos(alpha)];...
                                    [sin(alpha), sin(alpha), sin(alpha), sin(alpha)]];
sat.wheel.repartition_metrix_4RW_inverse = pinv(sat.wheel.repartition_matrix_4RW);