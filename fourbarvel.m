%Four bar velocity anslysis 
function [w_vec, VA, VBA, VB] = fourbarvel(l, th_vec_out, w_2)
    l1 = l(1);
    l2 = l(2);
    l3 = l(3);
    l4 = l(4);

    th_1 = th_vec_out(1);
    th_2 = th_vec_out(2);
    th_3 = th_vec_out(3);
    th_4 = th_vec_out(4);
    
    w_1 = 0;
    w_3 = (l2*w_2/l3) * (sind(th_4 - th_2)/sind(th_3 - th_4));
    w_4 = (l2*w_2/l4) * (sind(th_2 - th_3)/sind(th_4 - th_3));

    VA = l2*w_2 *(i*cosd(th_2) - sind(th_2));
    VAx = l2*w_2 * (-sind(th_2));
    VAy = l2*w_2 * cosd(th_2);

    VBA = l3*w_3 *(i*cosd(th_3) - sind(th_3));
    VB = l4*w_4 *(i*cosd(th_4) - sind(th_4));
    
    w_vec = [w_1, w_2, w_3, w_4];
    
end
