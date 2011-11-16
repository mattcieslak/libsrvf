function flag = isClosedCurve(curve)
% to determine the curve is closed or not. 
% for the closed curve, the first and the end point are the common;
% otherwise, it is open for our work. 
% history: 
%  - Mar, 18, 2010
% author: Hongquan Sun

head = curve(:, 1);
rear = curve(:, end); 
if all(head == rear) == 1
    flag = 1;
else 
    flag = 0;
end