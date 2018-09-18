%% setup running
H2O2 = [0.005,0.1,0.25,0.5,1,2];
for i = 1:length(H2O2)
    H2O2_Conc = H2O2(i); %mM
    OPD_Conc = 0.5; %mM
    Diff_Desc = [1,1];
    [Conc_record_1_col,Conc_record_1_row,Conc_record_2_col,Conc_record_2_row,H2O2_field,OPD_field,time,DAP_Coa]...
    = DiffusionSimulation2D(H2O2_Conc,OPD_Conc,Diff_Desc);
    save(['DiffBoth_Result_H2O2',num2str(H2O2_Conc*1e3),'uM_OPD_Conc_',num2str(OPD_Conc*1e3),'Complicated']);
end
%exit;

Diff_DM = [1 1;1 0;0 1;0 0];
for i = 1:length(Diff_DM)
    H2O2_Conc = 1; %mM
    OPD_Conc = 0.5; %mM
    Diff_Desc = Diff_DM(i,:);
    [Conc_record_1_col,Conc_record_1_row,Conc_record_2_col,Conc_record_2_row,H2O2_field,OPD_field,time2,time1,DAP_Coa]...
    = DiffusionSimulation2D(H2O2_Conc,OPD_Conc,Diff_Desc);
    save(['Diff',num2str(Diff_Desc(1)),num2str(Diff_Desc(2)),'_H2O2_',num2str(H2O2_Conc*1e3),'uM_OPD_Conc_',num2str(OPD_Conc*1e3),'FullKinetics']);
    disp('file saved');
end
disp('Done')
exit;