function [alpha_vec, AA, ABA, AB] = fourbaraccel(l, th_vec_out, w_vec_out, alpha_2)
    alpha_1 = 0;

    l1 = l(1);
    l2 = l(2);
    l3 = l(3);
    l4 = l(4);

    th_1 = th_vec_out(1);
    th_2 = th_vec_out(2);
    th_3 = th_vec_out(3);
    th_4 = th_vec_out(4);

    w_1 = w_vec_out(1);
    w_2 = w_vec_out(2);
    w_3 = w_vec_out(3);
    w_4 = w_vec_out(4);
    
    A = l4*sind(th_4);
    B = l3*sind(th_3);
    C = l2*alpha_2*sind(th_2) + l2*w_2*w_2*cosd(th_2) + l3*w_3*w_3*cosd(th_3) - l4*w_4*w_4*cosd(th_4);
    D = l4*cosd(th_4);
    E = l3*cosd(th_3);
    F = l2*alpha_2*cosd(th_2) - l2*w_2*w_2*sind(th_2) - l3*w_3*w_3*sind(th_3) + l4*w_4*w_4*sind(th_4);

    alpha_4 = (C*E - B*F) /(A*E - B*D);
    alpha_3 = (C*D - A*F) / (A*E - B*D); 
    alpha_vec = [alpha_1, alpha_2, alpha_3, alpha_4];
    
    Ax = -l2*alpha_2*sind(th_2) - l2*w_2*w_2*cosd(th_2);
    Ay = l2*alpha_2*cosd(th_2) - l2*w_2*w_2*sind(th_2);
    AA = Ax + j*Ay; 

    ABAx = -l3*alpha_3*sind(th_3) - l3*w_3*w_3*cosd(th_3);
    ABAy = l3*alpha_3*cosd(th_3) - l3*w_3*w_3*sind(th_3);
    ABA = ABAx + j*ABAy; 

    ABx = -l4*alpha_4*sind(th_4) - l4*w_4*w_4*cosd(th_4);
    ABy = l4*alpha_4*cosd(th_4) - l4*w_4*w_4*sind(th_4);
    AB = ABx + j*ABy;
    
    alpha_vec = [alpha_1, alpha_2, alpha_3, alpha_4];
end