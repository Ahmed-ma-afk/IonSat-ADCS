function [D] = distance_on_sphere(t1, p1, t2, p2)
D = acos(cos(t1)*cos(t2)+sin(t1)*sin(t2)*cos(p1-p2))
end

