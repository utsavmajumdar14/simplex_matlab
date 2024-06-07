clear all
clc
format short;

% Initialisation %

n=0; % No of variables
m=0; % No of constraints
Objin=""; % Max or Min
Cin=[]; % Cost vector
Ain=[]; % A matrix
ineq=[]; % Inequality matrix
bin=[]; % RHS matrix
C=[]; % Transposed and adjusted Cost vector (generalised to max)
flag=0; % To test for infeasibility, unbounded, optimal feasible solutions
chkminrat =[]; % Used to shift common min ratio in case of degeneracy
called = 0; % Checking number of times the Bland's rule function is called
BFS = []; % Storing Basic Feasible Solution


fprintf(" ********* SIMPLEX ********* \n");
fprintf(" This code can be run in 2 methods : Test mode and User Input \n");
fprintf(" Type 4 for user input or 0,1,2,3 for any of the test modes \n\n");

% Default inputs for test cases or user-input %

prompt0 = " Enter mode (Test for Infeasible : 0 / Feasible : 1 / \n Degeneracy : 2 / Unbounded : 3 / Example1 : 4 /\n Example2 : 5 / Example3 : 6 / Production Problem : 7 / User input : 8): \n";
mode = input(prompt0);
if (mode == 0)
    % Test case for infeasibility %
    n = 2;
    m = 5;
    Objin = "Min";
    Cin = [10;15];
    Ain = [2,3;3,-5;2,10;1,-2;3,-2];
    ineq = ['>';'>';'=';'<';'<'];
    bin = [7;-8;27;-6;3];

elseif (mode == 1)
    % Test case for Feasible solution %
    n = 3;
    m = 3;
    Objin = "Max";
    Cin = [3;-1;2];
    Ain = [1,3,1;2,-1,1;4,3,-2];
    ineq = ['<';'>';'='];
    bin = [5;2;5];

elseif (mode == 2)
    % Test case for degeneracy %
    n = 4;
    m = 3;
    Objin = "Max";
    Cin = [10;-57;-9;-24];
    Ain = [0.5,-5.5,-2.5,9;0.5,-1.5,-0.5,1;1,1,1,1];
    ineq = ['<';'<';'<'];
    bin = [0;0;1];

elseif (mode == 3)
    % Test case for unboundedness %
    n = 2;
    m = 2;
    Objin = "Max";
    Cin = [2;1];
    Ain = [1,-1;2,-1];
    ineq = ['<';'<'];
    bin = [10;40];

elseif (mode == 4)
    % Example1 %
    n = 2;
    m = 2;
    Objin = "Min";
    Cin = [-5;-7];
    Ain = [-3,2;-2,1];
    ineq = ['<';'<'];
    bin = [30;12];

elseif (mode == 5)
    % Example2 %
    n = 5;
    m = 10;
    Objin = "Max";
    Cin = [60;40;30;30;15];
    Ain = [1,1,1,1,1;4,2,2,2,1;0,1,0,1,0;1,0,0,0,0;0,0,1,0,0;1,1,1,0,0;0,0,0,1,1;0,1,0,0,0;0,0,0,1,0;0,0,0,0,1];
    ineq = ['<';'<';'<';'<';'<';'<';'<';'>';'>';'>'];
    bin = [7;8;3;1.8;0.3;3.8;3.2;0.5;0.5;0.4];

elseif (mode == 6)
    % Example3 %
    n = 4;
    m = 3;
    Objin = "Min";
    Cin = [-0.75;20;-0.5;6];
    Ain = [0.25,-8,-1,9;0.5,-12,-0.5,3;0,0,1,0];
    ineq = ['<';'<';'<'];
    bin = [0;0;1];

elseif (mode == 7)
    % Production Problem %
    fprintf("\n\n **** PRODUCTION PROBLEM **** \n\n");
    fprintf("\n We have 2 assemblies made out of 10 child parts : \n");
    fprintf(" Assembly 1 requires the following parts with respective quantities : \n");
    fprintf(" P01 : 2 \n P02 : 1 \n P03 : 0 \n P04 : 20 \n P05 : 1 \n P06 : 1 \n P07 : 3 \n P08 : 12 \n P09 : 1 \n P10 : 2 \n ");
    fprintf(" Assembly 2 requires the following parts with respective quantities : \n");
    fprintf(" P01 : 1 \n P02 : 3 \n P03 : 2 \n P04 : 15 \n P05 : 0 \n P06 : 2 \n P07 : 4 \n P08 : 5 \n P09 : 3 \n P10 : 1 \n ");
    fprintf(" Selling price for each Assembly 1 : 600$ \n");
    fprintf(" Selling price for each Assembly 2 : 700$ \n");
    fprintf(" We have the following stock for each part \n");
    fprintf(" P01 : 800 \n P02 : 500 \n P03 : 200 \n P04 : 6000 \n P05 : 500 \n P06 : 700 \n P07 : 2000 \n P08 : 3000 \n P09 : 900 \n P10 : 600 \n ");
    fprintf(" Our objective is to maximize sales revenue of both assemblies while adhering to inventory constraints \n");

    n = 2;
    m = 10;
    Objin = "Max";
    Cin = [600;700];
    Ain = [2,1;1,3;0,2;20,15;1,0;1,2;3,4;12,5;1,3;2,1];
    ineq = ['<';'<';'<';'<';'<';'<';'<';'<';'<';'<'];
    bin = [800;500;200;6000;500;700;2000;3000;900;600];

else
    % User input % 
    prompt = "Enter number of variables : \n";
    n = input(prompt); % Taking number of variables from user
    prompt2 = "Enter number of constraints : \n";
    m = input(prompt2); % Taking number of constraints from user
    
    prompt3 = "Enter Max or Min : \n";
    Objin = input(prompt3,"s"); % Taking Objective Max / Min input from user
    
    fprintf("Enter cost coefficients : \n");
    % Taking cost vector from user (as a column which will be transposed)
    for i=1:n
        fprintf("c%d : ",i);
        Cin(i,1) = input('');
    end
    % Taking A matrix coefficients from user
    for i=1:m
        fprintf('Enter variable coefficients for constraint no : %d \n',i);
        for j=1:n
            fprintf("a%d%d : ",i,j);
            Ain(i,j) = input('');
        end
        % Taking inequalities from user (</>/=)
        fprintf("Inequality %d (< / = / >) : ",i);
        ineq(i,1) = input('',"s");
        fprintf("b%d : ",i);
        bin(i,1) = input('');
    end
end
fprintf("Objective : ",Objin);
fprintf("\n");

fprintf("C = \n");
disp(Cin);

fprintf("A = \n");
disp(Ain);

fprintf("Inequalities : \n");
disp(char(ineq));
fprintf("\n");

fprintf("b = \n");
disp(bin);
    
% Standardising into Maximization %
if(Objin == "Min")
    for i=1:n
        C(1,i) = -Cin(i,1); % Transposing and reversing sign %
    end
else
    for i=1:n
        C(1,i) = Cin(i,1); % Transposing %
    end
end


% Standardising RHS negativities %
for i=1:m
    % Flipping inequalities %
    if(bin(i,1)<0)
        bin(i,1) = -bin(i,1);
        if(ineq(i,1)=='<')
            ineq2(i,1) = '>';
        end
        if(ineq(i,1)=='>')
            ineq2(i,1) = '<';
        end
        % Multiplying by -1 on LHS and RHS %
        for j = 1:n
            Ain(i,j) = -Ain(i,j);
        end
    else
        ineq2(i,1) = ineq(i,1);
    end
end

fprintf("New A = \n");
disp(Ain);

fprintf("New Inequalities : \n");
disp(char(ineq2));
fprintf("\n");

fprintf("New B = \n");
disp(bin);

% Creating Matrix for slack, surplus and artificial variables %
counter = 0;
for i=1:m
    if(ineq2(i,1)=='>')
        surpluscolumn(i,1) = -1; 
        counter = counter + 1;
    else
        surpluscolumn(i,1) = 0;
    end
end
fprintf("Number of > inequalities : %d \n",counter);
fprintf("Surplus column : \n");
disp(surpluscolumn);
surpluscolumn2 = zeros(m,counter);
for i=1:m
    surpluscolumn2(i,i) = surpluscolumn(i,1);
end
fprintf("Surplus column2 : \n");
disp(surpluscolumn2);

surpluscolumn2(:,all(surpluscolumn2 == 0))=[];
fprintf("Surplus column2 Changed : \n");
disp(surpluscolumn2);

I = eye(m);
A = [Ain surpluscolumn2 I];
AB = [A bin];

% Printing the re-adjusted problem %

fprintf("C' = \n");
disp(C);

fprintf("New Inequalities : \n");
disp(char(ineq2));
fprintf("\n");

fprintf("New b = \n");
disp(bin);

fprintf("After slack surplus artificial : \n A = \n");
disp(A);

fprintf("A:B = \n");
disp(AB);

% Locating columns with Basic Variables %
BVcol = [];
for j = 1:size(I,2)
    for i = 1:size(A,2)
        if A(:,i) == I(:,j)
            BVcol = [BVcol i];
        end
    end
end

fprintf("BV Columns : \n");
disp(BVcol);
B = A(: , BVcol);
fprintf('B = \n');
disp(B);

% Starting Simplex %
% Phase 1 %
% Identifying the rows with artificial %
artifrow = [];
countart = 0;
counteq = 0;
for i=1:m
    if(ineq2(i,1) == '>')
        artifrow = [artifrow i];
        countart = countart + 1;
    end
    if(ineq2(i,1) == '=')
        artifrow = [artifrow i];
        countart = countart + 1;
        counteq = counteq + 1;
    end
end

fprintf("Indicating rows that have artificial variables \n");
disp(artifrow);

Costphase1 = zeros(1,m+n+countart-counteq+1); % +1 to denote an extra 0 in the end to accomodate for RHS

fprintf("Phase1 Cost initialization : \n");
disp(Costphase1);

fprintf("No of surplus variables : \n");
disp(counter);

fprintf("No of artificial variables : \n");
disp(countart);

fprintf("No of equalities : \n");
disp(counteq);

for j=1:size(artifrow,2)
        Costphase1(1,n+counter+artifrow(j))=-1;
end

fprintf("Phase1 Cost : \n");
disp(Costphase1);

artifcol=[];
for i=1:size(artifrow,2)
    artifcol = [artifcol A(:,n+counter+artifrow(i))];
end

fprintf("Columns that have artificial variables : \n");
disp(artifcol);

OrigCost = [C Costphase1(:,n+1:size(Costphase1,2))];
fprintf("Original Cost with slack surplus and artificial : \n");
disp(OrigCost);

ArtificialColumns = artifrow + n + counter;
fprintf("Columns containing artificial variables in new cost vector : \n");
disp(ArtificialColumns);

% Finding Zj-Cj %

ZjCj = Costphase1(BVcol)*AB - Costphase1;

fprintf("ZjCj = \n");
disp(ZjCj);

disp(AB);

X = [ZjCj; AB];
fprintf("X = Zj-Cj and AB : \n");
disp(X);

InitialTab = array2table(X);
InitialTab
StartBV = BVcol;
fprintf('\n **** Phase1 starts **** \n')
checkedMinRatioRow = [];

% Iteration start %

RUN = true;
while RUN

    ZC = ZjCj(1,1:end-1); % Not considering the objective RHS so omitting last column
    fprintf("ZC = \n");
    disp(ZC);
    

    if any(ZC<0) % Checking if any negative value of reduced cost exists
        fprintf("Current solution is Not Optimal \n");
        % Entering Variable %
        [EnterCol, pvt_col] = min(ZC);
        fprintf("Entering Column = \n");
        disp(pvt_col);

        % Leaving Variable %
        sol = AB(:,end); 
        fprintf("Intermediate RHS : \n");
        disp(sol);
        Column = AB(:,pvt_col); 
        fprintf("Leaving column : \n");
        disp(Column);

        % Finding Minimum Ratio %

        for i=1:size(AB,1)
            if Column(i)>0
                ratio(i) = sol(i)./Column(i);
            else
                ratio(i) = inf;
            end
        end

        fprintf("Ratios (Transposed) : \n");
        disp(ratio);
        
        % Handling Degeneracy %
        countRatio = 0;

        [MinRatio, pvt_row] = min(ratio);
        
        for i=1:size(ratio,2)
            if(ratio(1,i)==MinRatio)
                countRatio = countRatio+1;
                checkedMinRatioRow = [checkedMinRatioRow, i];
            end
        end

        if(countRatio > 1) % If multiple min ratio exists
            fprintf("Degeneracy detected at rows : \n");
            disp(checkedMinRatioRow);
            % Record the min ratio rows and shift with each
            % iteration as per Bland's rule
            checkedMinRatioRow2 = blandrule(checkedMinRatioRow,called);
            disp(checkedMinRatioRow2);
            pvt_row = checkedMinRatioRow2(1,1);
            fprintf("Entering Row No : \n");
            disp(pvt_row);
            disp(checkedMinRatioRow2);
            checkedMinRatioRow = [];
                        
        end

        
        fprintf("Min Ratio Value = \n");
        disp(MinRatio);
        fprintf("Entering Row No : \n");
        disp(pvt_row);

        % Updating BFS %

        fprintf("Basic Variable columns before :\n");
        disp(BVcol);
        BVcol(pvt_row) = pvt_col;
    
        fprintf("Basic Variable columns after :\n");
        disp(BVcol);

        % Pivot Key %
        pvt_key = AB(pvt_row, pvt_col);
        fprintf("Pivot Key :\n");
        disp(pvt_key);

        fprintf("Old A:B Aug Matrix : \n");
        disp(AB);
        fprintf("Old Zj-Cj : \n");
        disp(ZjCj);

        % Updating Entries %
        AB(pvt_row,:) = AB(pvt_row,:)./pvt_key;
        for i=1:size(AB,1)
            if i~=pvt_row
                AB(i,:) = AB(i,:) - AB(i,pvt_col).*AB(pvt_row,:);
            end
        end
        ZjCj = ZjCj - ZjCj(pvt_col).*AB(pvt_row,:);

        fprintf("Old A:B Aug Matrix : \n");
        disp(AB);
        fprintf("New Zj-Cj : \n");
        disp(ZjCj);

        % Printing Table %
        ZCj = [ZjCj;AB];
        TABLE = array2table(ZCj);
        TABLE

        BFS(BVcol) = AB(:,end);
        fprintf("New BFS (Updated RHS) : \n");
        disp(BFS(BVcol));
    
    else
        RUN = false;
        if any (ismember(BVcol,ArtificialColumns)) == 1
            BFS = BVcol;
            fprintf("Final BFS at Phase 1 optimality : \n");
            disp(BFS);
            fprintf("Infeasible, since artificial variable lingering in system at optimality \n");
            flag = -1;
        else
            flag = 0;
            fprintf("Optimal Solution for Phase1 Reached \n");
            fprintf(" **** Phase1 End **** \n");
            BFS = BVcol;
            fprintf(" BFS columns for Phase 1 Optimality : \n");
            disp(BFS);
        end
    end
end

% Phase 2 starts %

if flag == 0 % Only enters if feasibility is achieved
    fprintf("\n ***** Phase 2 Starts ***** \n");
    countBV = 0;
   
    
    for i=1:size(BFS,2)
        for j=1:size(ArtificialColumns,2)
            if BFS(1,i)>ArtificialColumns(1,j)
                countBV = countBV + 1;
            end
        end
        BFS(1,i) = BFS(1,i)-countBV;
        countBV=0;
    end
    fprintf("New BFS Columns :\n");
    disp(BFS);

    AB(:,ArtificialColumns)=[];
    OrigCost(:,ArtificialColumns) = [];
    fprintf ("Phase2 Aug A:b Matrix with removed artificial columns : \n");
    disp(AB);
    fprintf("Cost vector after removing artificial variables : \n");
    disp(OrigCost);
end

% Phase 2 %
BVcol = BFS;
checkedMinRatioRow = [];
if flag==0
    ZjCj = OrigCost(BVcol)*AB - OrigCost;

    fprintf("ZjCj = \n");
    disp(ZjCj);

    disp(AB);

    X = [ZjCj; AB];
    fprintf("X = Zj-Cj and AB : \n");
    disp(X);

    InitialTab = array2table(X);
    InitialTab
    StartBV = BVcol;
    if all(ZjCj>=0)
        fprintf("Optimality Reached even before starting Phase2 \n");
        
        fprintf(" Optimal objective solution : \n");
        disp(ZjCj(1,end));
    else
        RUN = true;
        while RUN
        
            ZC = ZjCj(1,1:end-1);
            fprintf("ZC = \n");
            disp(ZC);

            if any(ZC<0)
                fprintf("Current solution is Not Optimal \n");
                % Entering Variable %
                [EnterCol, pvt_col] = min(ZC);
                fprintf("Entering Column = \n");
                disp(pvt_col);

                % Leaving Variable %
                sol = AB(:,end); 
                fprintf("Intermediate RHS : \n");
                disp(sol);
                Column = AB(:,pvt_col); 
                fprintf("Leaving column : \n");
                disp(Column);
                if Column<0
                    flag = 1;
                    fprintf("Unbounded Solution \n");
                    if(Objin=="Max")
                        fprintf("Optimal solution is %d \n",inf);
                        fprintf("Unbounded direction \n");
                        disp(C');
                    else
                        fprintf("Optimal solution is %d \n",-inf);
                        fprintf("Unbounded direction \n");
                        disp(-C');
                    end
                    RUN = false;
                    
                else
                    % Finding Minimum Ratio %

                    for i=1:size(AB,1)
                        if Column(i)>0
                            ratio(i) = sol(i)./Column(i);
                        else
                            ratio(i) = inf;
                        end
                    end

                    % Handling Degeneracy %
                    countRatio = 0;

                    [MinRatio, pvt_row] = min(ratio);
                    fprintf("Ratios (Transposed) : \n");
                    disp(ratio);
                    fprintf("Min Ratio Value = \n");
                    disp(MinRatio);
                    fprintf("Entering Row No : \n");
                    disp(pvt_row);

                    % Count number of times minimum ratio appears %
                    for i=1:size(ratio,2)
                        if(ratio(1,i)==MinRatio)
                            countRatio = countRatio+1;
                            checkedMinRatioRow = [checkedMinRatioRow, i];
                        end
                    end

                    if(countRatio > 1) % If multiple min ratio exists
                        fprintf("Degeneracy detected at rows : \n");
                        disp(checkedMinRatioRow);
                        % Record the min ratio rows and shift with each
                        % iteration as per Bland's rule
                        checkedMinRatioRow2 = blandrule(checkedMinRatioRow,called);
                        disp(checkedMinRatioRow2);
                        pvt_row = checkedMinRatioRow2(1,1);
                        fprintf("Entering Row No : \n");
                        disp(pvt_row);
                        disp(checkedMinRatioRow2);
                        checkedMinRatioRow = [];
                        
                    end

                    
                    fprintf("Ratios (Transposed) : \n");
                    disp(ratio);
                    fprintf("Min Ratio Value = \n");
                    disp(MinRatio);
                    fprintf("Entering Row No : \n");
                    disp(pvt_row);

                    % Updating BFS %

                    fprintf("Basic Variable columns before :\n");
                    disp(BVcol);
                    BVcol(pvt_row) = pvt_col;
    
                    fprintf("Basic Variable columns after :\n");
                    disp(BVcol);

                    % Pivot Key %
                    pvt_key = AB(pvt_row, pvt_col);
                    fprintf("Pivot Key :\n");
                    disp(pvt_key);

                    fprintf("Old A:B Aug Matrix : \n");
                    disp(AB);
                    fprintf("Old Zj-Cj : \n");
                    disp(ZjCj);

                    % Updating Entries %
                    AB(pvt_row,:) = AB(pvt_row,:)./pvt_key;
                    for i=1:size(AB,1)
                        if i~=pvt_row
                            AB(i,:) = AB(i,:) - AB(i,pvt_col).*AB(pvt_row,:);
                        end
                    end
                    ZjCj = ZjCj - ZjCj(pvt_col).*AB(pvt_row,:);

                    fprintf("Old A:B Aug Matrix : \n");
                    disp(AB);
                    fprintf("New Zj-Cj : \n");
                    disp(ZjCj);

                    % Printing Table %
                    ZCj = [ZjCj;AB];
                    TABLE = array2table(ZCj);
                    TABLE

                    BFS(BVcol) = AB(:,end);
                    fprintf("New BFS (Updated RHS) : \n");
                    disp(BFS(BVcol));

                end
            else
                RUN = false;
                flag = 2;
                fprintf("Optimal Solution for Phase2 Reached \n");
                fprintf(" **** Phase2 End **** \n");
                fprintf(" Basic Variable Columns at Phase 2 Optimality : \n");
                disp(BVcol);
                if(ZjCj(1,end)== inf)
                    fprintf(" BFS at Phase 2 Optimality : \n");
                    disp(BFS(BVcol));
                    fprintf(" Optimal objective solution : \n");
                    disp(ZjCj(1,end));
                    fprintf(" Unbounded Solution ");
                elseif (mode == 7)
                    fprintf(" BFS at Phase 2 Optimality : \n");
                    disp(round(BFS(BVcol)));
                    fprintf(" Optimal objective solution : \n");
                    disp(round(ZjCj(1,end)));
                else
                    fprintf(" BFS at Phase 2 Optimality : \n");
                    disp(BFS(BVcol));
                    fprintf(" Optimal objective solution : \n");
                    disp(ZjCj(1,end));
                end
            end
        end
    end
end

% Shifting min ratio row by Bland's rule %

function [NewRatioList] = blandrule(chkminrat,called)
    NewRatioList = chkminrat;
    if(called == 0)
        NewRatioList;
    else
        for i=0:called+1
            NewRatioList(:,1)=[] ;
        end
    end
    called = called + 1;
end