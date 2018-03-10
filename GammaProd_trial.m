function fnl = GammaProd_trial(p,x)
   fnl  = 1;
    
    for i = 1:length(p)
       fnl  = gamma(p(i) - x(i))*fnl;
    end
end 
    
    