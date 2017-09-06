function [C,S]=transcendentalfunc(y)

    if y>0
        sqrty=sqrt(y);
        C=(1-cos( sqrty ))/y;
        S=(sqrt(y)-sin(sqrt(y)))/sqrt(y^3);
        
    elseif y<0
        sqrtmy=sqrt(-y);
            C=(cosh(sqrtmy)-1)/(-y);
            S=(sinh(sqrtmy)-sqrtmy)/sqrt(-y^3);
    else
       C=0.5;
       S=0.1667;
        
    end


end