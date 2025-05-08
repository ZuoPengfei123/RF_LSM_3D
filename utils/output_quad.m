function [quadrant] = output_quad(back_azimuth)
% X1,Y1 : any point
% X0,Y0 : referrence point

    if back_azimuth > 0 && back_azimuth <= 90
        quadrant = 1;
    elseif back_azimuth > 90 && back_azimuth <= 180
        quadrant = 2;
    elseif back_azimuth > 180 && back_azimuth <= 270
        quadrant = 3;
    elseif back_azimuth > 270 && back_azimuth <= 360
        quadrant = 4;
    end

end
