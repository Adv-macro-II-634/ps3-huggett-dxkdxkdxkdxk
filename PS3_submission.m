

beta = .994; 
sigma = 1.5; 
b = 0.5; 
y_s = [1, b]; 
Pimx = [.97 .03; .5 .5]; 
Pimx1 = Pimx ^ 100;
a_lo = -2; 
a_hi = 5;
num_a = 101;

a = linspace(a_lo, a_hi, num_a); 

q_min = 0; 
q_max = 1;
q_guess = (q_min + q_max) / 2;  % initial guess

n=1;
aggsav = 1 ;
while abs(aggsav) >= 0.01 ;
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); 
						
ret(cons < 0) = -Inf;

v_guess = zeros(2, num_a);
   
v_tol = 1;
i=1;

while v_tol >.0001;
    
    value_fn = ret + beta * ...
        repmat(permute((Pimx * v_guess), [3 2 1]), [num_a 1 1]);
    
    [vfn, pol_indx] = max(value_fn, [], 2);
    vfn = permute(vfn,[3 1 2]);
    
        
    v_tol = max(abs(v_guess(:) - vfn(:)));
    v_guess = vfn;
    i=i+1;
 end

pol_indx = permute(pol_indx, [3 1 2]);
g = a(pol_indx); 
Mu = (1/(2*num_a))*ones(2,num_a);
Mu_tol   = 1;

while Mu_tol > 1e-8;   
    
    [emp_ind, a_ind] = find(Mu > 0); 
        
          MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... 
            (Pimx(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)))';
    end
        
        Mu_tol   = max(abs(MuNew(:)-Mu(:)));
        Mu= MuNew;
        
    end
   
    MarketClearing=sum(sum(Mu.*g));
    
   
    if MarketClearing>= 0.01% price is too low ; too much savinf
        q_min = (q_min + q_max) / 2; %new q_min must be higher than before   
    else
        q_max = (q_min + q_max) / 2;   
    end 
    q_guess = (q_min + q_max) / 2;
    aggsav = MarketClearing;
    n=n+1;
    q=q_guess;
    
    disp('market-clearing; iteration; price')
    disp(full([MarketClearing n q]))
end


%Lorenz Curve

 income_E = a+1;
 income_U = a+b;
 Income = [income_E income_U].*Mu(:)';
 TotalIncome=sum(Income(:));
 [Income sort] = sort(Income);

 perc = Income/TotalIncome;
 pop = Mu(:)';
 pop = pop(sort);

 for i = 2:length(Income)
     perc(i) = perc(i)+perc(i-1);
     pop(i) = pop(i)+pop(i-1);
 end
 
 plot([0 1],[0 1],'r',pop,perc,'--','Linewidth',1);
 title('Lorenz Curve')
 xlabel('Cumulative Share from Lowest to Highest Income') 
 ylabel('Cumulative Share of Income')
 legend({'Equality line','Lorenz Curve'},'Location','southeast')