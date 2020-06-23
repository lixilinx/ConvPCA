function w = coswin( F )
% return a cosine window with length F
w = [1:F] - 0.5*(F+1);
w = 0.5 + 0.5*cos(2*pi*w/F);