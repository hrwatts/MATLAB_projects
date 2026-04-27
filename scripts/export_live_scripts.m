function export_live_scripts()
%EXPORT_LIVE_SCRIPTS Export selected live scripts to HTML.

targets = {
    'examples/tutorial/matlab_tutorial.mlx'
    'examples/tutorial/matlab_principles.mlx'
    'examples/bifurcation/bifurcation.mlx'
    'examples/epidemiology/SIR_fit_AL.mlx'
};

for idx = 1:numel(targets)
    input_file = targets{idx};
    [folder, name] = fileparts(input_file);
    output_file = fullfile(folder, [name '.html']);
    matlab.internal.liveeditor.openAndConvert(input_file, output_file);
end
end
