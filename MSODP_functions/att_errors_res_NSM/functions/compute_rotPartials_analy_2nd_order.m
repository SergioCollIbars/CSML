function [H] = compute_rotPartials_analy_2nd_order(T)
    S1 = [0,0,0;0,0,1;0,-1,0];
    S2 = [0,0,-1;0,0,0;1,0,0];
    S3 = [0,1,0;-1,0,0;0,0,0];
    
    % square partials terms
    dT_d1d1 = S1 * S1 * T + S1 * T * S1' + S1 * T * S1' + T * (S1 * S1)';
    dT_d2d2 = S2 * S2 * T + S2 * T * S2' + S2 * T * S2' + T * (S2 * S2)';
    dT_d3d3 = S3 * S3 * T + S3 * T * S3' + S3 * T * S3' + T * (S3 * S3)';
    
    % symetric partials terms
    dT_d3d1 = S1 * S3 * T + S3 * T * S1' + S1 * T * S3' + T * (S1 * S3)';
    dT_d3d2 = S2 * S3 * T + S3 * T * S2' + S2 * T * S3' + T * (S2 * S3)';
    dT_d2d1 = S1 * S2 * T + S2 * T * S1' + S1 * T * S2' + T * (S1 * S2)';

    h_d1d1 = [dT_d1d1(1,1);dT_d1d1(1,2);dT_d1d1(1,3);dT_d1d1(2,1);...
        dT_d1d1(2,2);dT_d1d1(2,3);dT_d1d1(3,1);dT_d1d1(3,2);dT_d1d1(3,3)];
    h_d2d2 = [dT_d2d2(1,1);dT_d2d2(1,2);dT_d2d2(1,3);dT_d2d2(2,1);...
        dT_d2d2(2,2);dT_d2d2(2,3);dT_d2d2(3,1);dT_d2d2(3,2);dT_d2d2(3,3)];
    h_d3d3 = [dT_d3d3(1,1);dT_d3d3(1,2);dT_d3d3(1,3);dT_d3d3(2,1);...
        dT_d3d3(2,2);dT_d3d3(2,3);dT_d3d3(3,1);dT_d3d3(3,2);dT_d3d3(3,3)];

    h_d3d1 = [dT_d3d1(1,1);dT_d3d1(1,2);dT_d3d1(1,3);dT_d3d1(2,1);...
        dT_d3d1(2,2);dT_d3d1(2,3);dT_d3d1(3,1);dT_d3d1(3,2);dT_d3d1(3,3)];
    h_d3d2 = [dT_d3d2(1,1);dT_d3d2(1,2);dT_d3d2(1,3);dT_d3d2(2,1);...
        dT_d3d2(2,2);dT_d3d2(2,3);dT_d3d2(3,1);dT_d3d2(3,2);dT_d3d2(3,3)];
    h_d2d1 = [dT_d2d1(1,1);dT_d2d1(1,2);dT_d2d1(1,3);dT_d2d1(2,1);...
        dT_d2d1(2,2);dT_d2d1(2,3);dT_d2d1(3,1);dT_d2d1(3,2);dT_d2d1(3,3)];

    H = 0.5 * [h_d3d3, 2.*h_d3d2, 2.*h_d3d1, h_d2d2, 2.*h_d2d1, h_d1d1];
end

