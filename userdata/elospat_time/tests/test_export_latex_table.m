function test_export_latex_table()
% test_export_latex_table  Focused validation for Task 34 helper.

fprintf('\n========================================\n');
fprintf('TEST: export_latex_table (Task 34)\n');
fprintf('========================================\n');

tbl = table();
tbl.Method = ["Dense Arc"; "RMS-ELO-SPAT"];
tbl.dose_targetD95Min = [73.1; 72.5];
tbl.dose_targetD98Min = [71.0; 70.4];
tbl.dose_oarDmaxMax = [75.2; 74.3];
tbl.delivery_nActiveLayers = [40; 18];
tbl.delivery_nActiveSpots = [2000; 900];
tbl.delivery_nSwitchUp = [12; 5];
tbl.delivery_T_energy = [50.4; 18.2];
tbl.delivery_T_spot = [4.3; 2.1];
tbl.delivery_T_gantry = [12.0; 12.0];
tbl.delivery_T_total = [80.0; 36.5];

latexText = export_latex_table(tbl);

assert(contains(latexText, '\begin{tabular}'));
assert(contains(latexText, 'Dense Arc'));
assert(contains(latexText, 'RMS-ELO-SPAT'));
assert(contains(latexText, 'CTV D95'));
assert(contains(latexText, 'T\_total'));

fprintf('[PASS] LaTeX table export rendered expected content.\n');
end
