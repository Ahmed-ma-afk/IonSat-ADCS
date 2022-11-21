alpha = pi / 6;
sat.wheel.repartition_matrix_4RW = [[sin(alpha), sin(alpha), sin(alpha), sin(alpha)],
    [cos(alpha), 0, -cos(alpha), 0],
    [0, cos(alpha), 0, -cos(alpha)]];
sat.wheel.repartition_metrix_4RW_inverse = pinv(sat.wheel.repartition_matrix_4RW);
