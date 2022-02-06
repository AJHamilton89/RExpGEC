%no. Iridis Arrays
no_arrays=10;

%%iridis HPC functionality

array = str2double(getenv('SLURM_ARRAY_TASK_ID'));

if isempty(array)
    array = randi([0 no_arrays-1]);
end

fn1=sprintf('Figures/test=%i',array);
save(fn1)