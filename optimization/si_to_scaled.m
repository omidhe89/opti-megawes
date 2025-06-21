function scaled_w = si_to_scaled(w,lbw,ubw)
    scaled_w = (w-lbw)./(ubw-lbw);
end