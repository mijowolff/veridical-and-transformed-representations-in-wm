% simple p-value calculation from a previously obtained null-distribution

% written by Michael J. Wolff, Oxford, 2020

function p = FastPvalue(obs,null_distr,tail)

dim_sum=length(size(null_distr));

p=sum(null_distr>obs,dim_sum)/length(null_distr);

if tail==-1
    p=1-p;
elseif tail==2
    p(p>=.5) = 1-p(p>=.5);
    p = p*2;
end
p(p == 1) = (length(null_distr) - 1)/length(null_distr);
p(p == 0) = 1/length(null_distr);
end