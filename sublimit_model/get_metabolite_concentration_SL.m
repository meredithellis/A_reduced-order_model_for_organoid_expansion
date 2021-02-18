%% Get metabolite concentration profile data

%prolif_possibilities=cell(4,3);
prolif_possibilities{1,1}='low';
prolif_possibilities{2,1}='mid';
prolif_possibilities{3,1}=2;
prolif_possibilities{4,1}='high';

%uptake_possibilities=cell{3,1};
uptake_possibilities{1,1}='low';
uptake_possibilities{2,1}='mid';
uptake_possibilities{3,1}='high';

%C_functions=cell(4,3);
%W_functions=cell(4,3);
%t_functions=cell(4,3);
%%
for i=1:4
    proliferation_rate=prolif_possibilities{i,1};
    
    for j=1:3
        uptake_rate=uptake_possibilities{j,1};
        
        [t,X, c,w]=metabolite_concentration_SL(proliferation_rate, uptake_rate);
        
       
        C_functions{i,j} = c;
        W_functions{i,j} = w;
        t_functions{i,j} = t;
        X_functions{i,j} = X;
    end
end

%%
C_SL_LL = C_functions{1,1};
C_SL_ML = C_functions{2,1};
C_SL_2L = C_functions{3,1};
C_SL_HL = C_functions{4,1};
C_SL_LM = C_functions{1,2};
C_SL_MM = C_functions{2,2};
C_SL_2M = C_functions{3,2};
C_SL_HM = C_functions{4,2};
C_SL_LH = C_functions{1,3};
C_SL_MH = C_functions{2,3};
C_SL_2H = C_functions{3,3};
C_SL_HH = C_functions{4,3};

W_SL_LL = W_functions{1,1};
W_SL_ML = W_functions{2,1};
W_SL_2L = W_functions{3,1};
W_SL_HL = W_functions{4,1};
W_SL_LM = W_functions{1,2};
W_SL_MM = W_functions{2,2};
W_SL_2M = W_functions{3,2};
W_SL_HM = W_functions{4,2};
W_SL_LH = W_functions{1,3};
W_SL_MH = W_functions{2,3};
W_SL_2H = W_functions{3,3};
W_SL_HH = W_functions{4,3};

t_SL = t_functions{1,1};
X_SL = X_functions{1,1};
% t_SL_LL = t_functions{1,1};
% t_SL_ML = t_functions{2,1};
% t_SL_2L = t_functions{3,1};
% t_SL_HL = t_functions{4,1};
% t_SL_LM = t_functions{1,2};
% t_SL_MM = t_functions{2,2};
% t_SL_2M = t_functions{3,2};
% t_SL_HM = t_functions{4,2};
% t_SL_LH = t_functions{1,3};
% t_SL_MH = t_functions{2,3};
% t_SL_2H = t_functions{3,3};
% t_SL_HH = t_functions{4,3};
% 
% X_SL_LL = X_functions{1,1};
% X_SL_ML = X_functions{2,1};
% X_SL_2L = X_functions{3,1};
% X_SL_HL = X_functions{4,1};
% X_SL_LM = X_functions{1,2};
% X_SL_MM = X_functions{2,2};
% X_SL_2M = X_functions{3,2};
% X_SL_HM = X_functions{4,2};
% X_SL_LH = X_functions{1,3};
% X_SL_MH = X_functions{2,3};
% X_SL_2H = X_functions{3,3};
% X_SL_HH = X_functions{4,3};


%% Save data
% 
%     filename=sprintf('metabolite_concentration_data.mat');
%     save(filename, 'C_functions','W_functions','t_functions', 't')
    
    %%
   clearvars -except C* W* t* X* prolif_possibilities uptake_possibilities
   %%
   filename=sprintf('metabolite_concentration_SL_data.mat');
    save(filename)