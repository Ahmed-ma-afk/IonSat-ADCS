function [R] = interpol(t,p)
    norm = 0 ;
    sum = 0;
    s=0.1;
    for i=1:5
        t2 = Ts(i,0);
        p2 = Ps(i,0);
        c = Cs(i,0);
        dist = distance_on_sphere(t,p,t2,p2);
        if dist<1e-10
            R = c;
        end  
        coef = (1/dist)^s;
        sum = sum + coef*c;
        R = sum/norm;
    end
end