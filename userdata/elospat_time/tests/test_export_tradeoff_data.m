function test_export_tradeoff_data()
% test_export_tradeoff_data  Focused validation for Task 35 helper.

fprintf('\n========================================\n');
fprintf('TEST: export_tradeoff_data (Task 35)\n');
fprintf('========================================\n');

tbl = table();
tbl.Method = ["Dense Arc"; "RMS-ELO-SPAT"];
tbl.dose_targetD95Min = [73.1; 72.5];
tbl.delivery_T_total = [80.0; 36.5];
tbl.dose_oarDmaxMax = [75.2; 74.3];
tbl.delivery_nActiveLayers = [40; 18];
tbl.delivery_nActiveSpots = [2000; 900];
tbl.solver_iterations = [500; 180];

outDir = fullfile(tempdir, 'elospat_tradeoff_test');
if exist(outDir, 'dir') == 7
    rmdir(outDir, 's');
end
mkdir(outDir);

outputs = export_tradeoff_data(tbl, 'unit_tradeoff', outDir);

assert(istable(outputs.tradeoffTable));
assert(height(outputs.tradeoffTable) == 2);
assert(all(outputs.tradeoffTable.estimatedBDT == tbl.delivery_T_total));
assert(exist(outputs.csvFile, 'file') == 2);
assert(exist(outputs.matFile, 'file') == 2);

fprintf('[PASS] Tradeoff data exported to MAT and CSV.\n');
end
