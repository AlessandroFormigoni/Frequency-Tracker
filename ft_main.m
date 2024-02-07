clear all;
clc;
sel = menu('Scegli la modalità', 'Bittanti Discreto', 'Bittanti Campionato');

if sel == 1
    ft_run_disc
elseif sel == 2
    ft_run_sampl
end
