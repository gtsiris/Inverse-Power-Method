%Georgios Tsiris, 1115201700173
clear; clc;
format long g;

% Εφαρμογή της μεθόδου των δυνάμεων για την προσέγγιση της 
% μεγαλύτερης κατά μέτρο ιδιοτιμής και του αντίστοιχου ιδιοδιανύσματος του πίνακα Α


A=[16 -8 2 1; 2 -12 1 0; -1 1 -4 1; 0 -1 2 3];
n=4;nn=4; 
[n nn]=size(A);

y=[1 1 1 1]';
nz=1;n1=4;
[n1 nz]=size(y);

maxiter=100;
tol=(1/2)*10^-6;  %1e-8;

if (n~=nn)
    disp('Λάθος της μεθόδου δυνάμεων: Ο πίνακας Α πρέπει να είναι τετραγωνικός');return;
elseif (n~=n1)
    disp('Λάθος της μεθόδου δυνάμεων: οι διαστάσεις του πίνακα Α και του διανύσματος y δεν είναι συμβατές');return;
end;

%Ακριβείς ιδιοτιμές και ιδιοδιανύσματα
[eigvectors, D]=eig(A);
eigvalues=diag(D)
eigvectors
% disp('idiotimes tou A'); disp(eigvalues);
% disp('idiodianusmata tou A'); disp(eigvectors);

eigval_max=max(abs(eigvalues))
pos=min(find(abs(eigvalues-eigval_max)<tol))
eigvec_max=eigvectors(:,pos)
eigval_min=min(abs(eigvalues))
pos=min(find(abs(eigvalues-eigval_min)<tol))
eigvec_min=eigvectors(:,pos)


disp("\n\n");
disp('******Find lambda_max and z_max******');
% προσέγγιση της μεγαλύτερης κατά μέτρο ιδιοτιμής και του αντίστοιχου
% ιδιοδιανύσματος
for q=[(eigval_max-0.3):-0.1:(eigval_max-0.9)]
    disp('*************************************');
    q
    [lambda_temp,z_max] = eig_power(inv(A-q*eye(n)),y,tol,maxiter);
    lambda_max = (1/lambda_temp)+q
    error_at_lambda_max = abs(eigval_max-lambda_max)
    z_max
    error_at_z_max = norm(eigvec_max - z_max)
end


disp("\n\n");
disp('******Find lambda_min and z_min******');
% προσέγγιση της μικρότερης κατά μέτρο ιδιοτιμής και του αντίστοιχου
% ιδιοδιανύσματος

for q=[(eigval_min+0.3):0.1:(eigval_min+0.9)]
    disp('*************************************');
    q
    [lambda_temp,z_min] = eig_power(inv(A-q*eye(n)),y,tol,maxiter);
    lambda_min = (1/lambda_temp)+q
    error_at_lambda_min = abs(eigval_min-lambda_min)
    z_min
    error_at_z_min = norm(eigvec_min - z_min)
end
