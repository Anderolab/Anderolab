eps = [300*30, 15*30, 30*30, 180*30];
g_frec = [10, 20, 40, 60, 80];
t_rat = [4, 2, 1, .5, .1];
pt_rat = .5;
Neurons = {};
for i = 1:10
    for g = 1:length(g_frec)
        basal = g_frec(g);
        for t = 1:length(t_rat)
            ep_f = [basal, basal.*pt_rat, basal.*t_rat(t), basal];
            Neurons{g, t, i} = F_SyntNeur(eps, ep_f);
        end
    end
end
save("SampleNeurons.mat", "Neurons")