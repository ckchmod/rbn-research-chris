function lambda = lyapunov( F )

    tempA = bnActivity(F);
    temp = sum(tempA,2);
    lambda = log(mean(temp))

end

