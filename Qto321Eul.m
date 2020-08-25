function [psi theta phi] = Qto321Eul(q);

q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

yaw_num = (2*(q0*q1 + q3*q2));
yaw_denum = (q3^2 - q2^2 - q1^2 + q0^2);
phi = atan2(yaw_num,yaw_denum);

pitch_num = (-2*(q0*q2 - q1*q3));
theta = asin(pitch_num);

roll_num = (2*(q1*q2 + q3*q0));
roll_denum = (q3^2 + q2^2 - q1^2 - q0^2);
psi = atan2(roll_num,roll_denum);