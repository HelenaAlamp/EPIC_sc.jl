fileID = fopen("test_20MeV_0phi_pbeam_500turns.CORRECT.txt",'r');

formatSpec = '%f';
sizeA = [100000 100000];

A = fscanf(fileID,formatSpec,sizeA);

X=A(:,1);
Y=A(:,2);

plot(X,Y,'b','Marker','.','LineStyle','none');
xlabel('z(m)');
ylabel('dp/p_0')

grid on

% %%%%%%%%%%
% 
% fileID = fopen("test_500turns_20MeV_pbeam_correct.txt",'r');
% 
% formatSpec = '%f';
% sizeA = [100000 100000];
% 
% A = fscanf(fileID,formatSpec,sizeA);
% 
% X=A(:,1);
% Y=A(:,2);
% 
% hold on
% plot(X,Y,'r','Marker','.','LineStyle','none');
% xlabel('z(m)');
% ylabel('dp/p_0')
% 
% grid on
% 
% %%%%
% 
% fileID = fopen("test_500turns_30MeV.txt",'r');
% 
% formatSpec = '%f';
% sizeA = [100000 100000];
% 
% A = fscanf(fileID,formatSpec,sizeA);
% 
% X=A(:,1);
% Y=A(:,2);
% 
% hold on
% plot(X,Y,'g','Marker','.','LineStyle','none');
% xlabel('z(m)');
% ylabel('dp/p_0')


%%%%%%%%%%%%%%

% fileID = fopen("test_500turns_25MeV_10ph.txt",'r');
% 
% formatSpec = '%f';
% sizeA = [100000 100000];
% 
% A = fscanf(fileID,formatSpec,sizeA);
% 
% X=A(:,1);
% Y=A(:,2);
% 
% hold on
% plot(X,Y,'d','Marker','.','LineStyle','none');
% xlabel('z(m)');
% ylabel('dp/p_0')

Position = {cursor_info(1).Position}.';
Position1 = {cursor_info(2).Position}.';
v = cell2mat(Position(1));
v1 = cell2mat(Position1(1));
sigma_l = v-v1

