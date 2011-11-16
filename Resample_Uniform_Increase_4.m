function [shape] = Resample_Uniform_Increase_4(polygon, NUM_SAMPLES)
% resample from a set of points for open and closed curve
% output:
%  - shape: 2 * m, row 1: x coordinate, row 2: y coordinate.
% 
% 
% history:
%  - Feb, 23, 2010: version 1, only for open curve
%  - Feb, 23, 2010: version 2, to be for both the open and closed curves
%  - Mar, 2, 2010: version 3, the bug will generate more than the NUM_SAMPLES points is removed.
%  - Mar, 25, 2010: version 4, the bug of the open curve that miss the
%  last point is solved. And the open or closed curve is determined before
%  resampling. Also the one line to tell 'residual < interval' is changed
%  into 'residual <= interval'. But it need further verification.
% Note that whatever the open or closed curve, the head and rear point are
% different. It means, if the closed curve need to make the head and rear
% the same, the user should add the head point up to the rear of the
% returned shape.
%
% author: Hongquan Sun


isShow = 0;
if isempty(polygon)
    shape = [];
    return;
elseif size(polygon, 2) == 1
    shape = polygon;
    return;
end

if isShow
    figure, plot(polygon(1,:), polygon(2,:), 'r.'); hold on;
    plot(polygon(1,:), polygon(2,:), 'b');
    plot(polygon(1,:), polygon(2,:), 'c');
end
% polygon(1,:) = Xnew;
% polygon(2,:) = Ynew;

num = size(polygon, 2);
polyLength = 0;
for i = 2:num
    polyLength = polyLength + sqrt((polygon(1,i) - polygon(1,i-1))^2 + (polygon(2,i) - polygon(2,i-1))^2);
end

if isClosedCurve(polygon)
    unit = polyLength / NUM_SAMPLES;
else
    unit = polyLength / (NUM_SAMPLES - 1);
end


% two issue, the ratio determine the direction
% the unit determine the points to pass
shape(:,1) = polygon(:,1);
start = polygon(:,1);
if isShow
    plot(shape(1,1), shape(2,1), 'g*');
end
k = 1;
% t = 0;
i = 1;
residual = unit;
while k < NUM_SAMPLES
    if i < num % to exeption implementation, i+1 could not exeed num
        interval = sqrt((polygon(1,i+1) - polygon(1,i))^2 + (polygon(2,i+1) - polygon(2,i))^2);
    else
        if k < NUM_SAMPLES % specifically, the residual and the interval could not exactly equal, but with a very slight difference, and it could not distinct from residual <= interval, such that we nned to add the last point for the sampling.
            shape(:,end + 1) = polygon(:, end);
            if isShow
                plot(shape(1,end),shape(2,end), 'g*');
            end
            break;
        else
            break;
        end        
    end
    if residual <= interval
        cost = (polygon(1,i+1) - polygon(1,i)) / interval;
        sint = (polygon(2,i+1) - polygon(2,i)) / interval;        
        while residual <= interval
            if polygon(1, i+1) ~= polygon(1,i)
                x = start(1) + residual * cost;
                y = start(2) + residual * sint;
            else
                x = start(1);
                y = start(2) + residual * sign(polygon(2,i+1) - polygon(2,i));
            end
            k = k + 1;
            if k <= NUM_SAMPLES
                shape(1,k) = x;
                shape(2,k) = y;
            else
                break;
            end
            if isShow
                plot(x,y, 'g*');
            end
            residual = residual + unit;
%             pause;
        end
    else
        i = i + 1;
        start = polygon(:,i);
        residual = residual - interval;
    end
end

