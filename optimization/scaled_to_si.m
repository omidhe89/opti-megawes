function w = scaled_to_si(scaled_w,lbw,ubw)
    w = scaled_w .*(ubw-lbw) + lbw;
end