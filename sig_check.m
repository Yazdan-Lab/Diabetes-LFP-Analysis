function [hA, hB] =sig_values = (age, db)
age_sig = sig_check(age);
db_sig = sig_check(db);
if age <= 0.05
hA = suplabel(['age effect: ' age_sig],'t',[0.8 0.08 0.01 0.8]);
set(A,'FontSize',12)
end
if db <= 0.05
hB = suplabel(['db effect: ' db_sig],'t',[0.8 0.08 0.01 0.75]);
set(B,'FontSize',12)
end

function sig = sig_check(P)
if P > 0.05
    sig = 'n.s.';
elseif P <= 0.05 && P >0.01
    sig = '*';
elseif P <= 0.01 && P >0.001
    sig = '**';
elseif P <= 0.001
    sig = '***';
end
end