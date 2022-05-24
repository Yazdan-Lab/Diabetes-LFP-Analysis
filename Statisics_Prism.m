% Manually copying data from Prism files to create figures.
    
age_labs = [200 200 200 200 200 200 200 200 200 200 200 200 200 400 400 400 400 400 400 400 400 400 400]';
db_labs = [{'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'}]';

alpha_age_labs = [200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 ]';
alpha_db_labs = [{'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'}]';

beta_age_labs = [200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400]';
beta_db_labs = [{'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} ]';

mi_age_labs = [200 200 200 200 200 200 200 200 200 200 200 200 200 400 400 400 400 400 400 400 400 400]';
mi_db_labs = [{'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'DB'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'ctrl'} {'DB'} {'DB'} {'DB'} {'DB'}]';

%% CA1 Signal Power
hgamma_vals = [1.45 2.17 1.95 2.53 2.61 2.47 2.35 2.2 2.474356 2.519189 2.5 2.52 2.22 1.98 2.18 2.3 1.75 2.04 1.51 2.06 1.75 2.37 2.23]';

[hg_P, hg_T, hg_Stats] = anovan(hgamma_vals,{db_labs, age_labs}, 'model', 'interaction','display','off');
[hg_C,hg_M,~,hg_N] = multcompare(hg_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(hg_M(:,2), hg_M(:,1), hg_C);
sig_values(hg_P(2), hg_P(1));
title('High Gamma')
ylabel('Signal power')


gamma_vals = [1.7 1.91 1.95 2.53 2.61 2.47 2.35 2.2 2.474356 2.519189 2.27 2.52 2.22 2.11 2.08 2.29 1.75 2.04 1.56 1.62 1.75 2.37 2.23]';

[g_P, g_T, g_Stats] = anovan(gamma_vals,{db_labs, age_labs}, 'model', 'interaction','display','off');
[g_C,g_M,~,g_N] = multcompare(g_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(g_M(:,2), g_M(:,1), g_C);
sig_values(g_P(2), g_P(1));
title('Gamma')
ylabel('Signal power')


beta_vals = [1.34950482 1.55844095 1.74690263 1.81025 1.99365971 1.98514479 1.94581509 1.58091928 1.58761893 1.57640586 2.09187298 2.16557667 2.17884002 2.0363426 1.87704979 1.54918426 1.41560703 1.90251585 1.98142855 1.58822292 1.2177576 1.93655496 2.2634429 1.69193156 1.73124135 1.96524546 2.08606997 1.75531448 1.76712172 1.70379848 1.82157706 1.29212255 1.31918659 1.30799712 1.68278633 1.66416467 1.35539404 1.55342301 1.40265164 1.67008064]';

[b_P, b_T, b_Stats] = anovan(beta_vals,{beta_db_labs, beta_age_labs}, 'model', 'interaction','display','off');
[b_C,b_M,~,b_N] = multcompare(b_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(b_M(:,2), b_M(:,1), b_C);
sig_values(b_P(2), b_P(1));
title('Beta')
ylabel('Signal power')


alpha_vals = [1.27309623 1.59222674 1.60791267 1.88288355 1.97323087 2.05837211 1.90369757	1.37201488 1.4943383 1.5657109 2.03637096 1.96390394 1.52369705 1.74371682 1.80092794 1.77118182 1.90571273 2.02184647 1.75802724 1.42493225 1.89729162 2.29471306 1.90265819 0.58575392 1.83079646 1.80382843 1.65021716 1.69384474 1.82683254 1.64497891 1.85646243 1.22806396 1.39269144 1.35052349 1.6957178 1.6488733 1.34416033 1.57941886 1.36890209]';

[a_P, a_T, a_Stats] = anovan(alpha_vals,{alpha_db_labs, alpha_age_labs}, 'model', 'interaction','display','off');
[a_C,a_M,~,a_N] = multcompare(a_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(a_M(:,2), a_M(:,1), a_C);
sig_values(a_P(2), a_P(1));
title('Alpha')
ylabel('Signal power')

theta_vals = [1.45 2.17 1.95 2.53 2.61 2.47 2.35 2.2 2.474356 2.519189 2.5 2.52 2.22 1.98 2.18 2.3 1.75 1.99 1.47 2.06 1.75 2.37 2.23]';

[t_P, t_T, t_Stats] = anovan(theta_vals,{db_labs, age_labs}, 'model', 'interaction','display','off');
[t_C,t_M,~,t_N] = multcompare(t_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(t_M(:,2), t_M(:,1), t_C);
sig_values(t_P(2), t_P(1));
title('Theta')
ylabel('Signal power')

%% Modulation index
CA1_ctx_vals = [-3.41 -2.93 -3 -3.43 -3.77 -3.02 -3.35 -2.91 -2.248 -2.95731 -3.1 -3.53 -2.42 -3.02 -3.05 -2.51 -3.19 -2.89 -2.76 -1.92 -2.1 -2.62]';

[mp_P, mp_T, mp_Stats] = anovan(CA1_ctx_vals,{mi_db_labs, mi_age_labs}, 'model', 'interaction','display','off');
[mp_C,mp_M,~,mp_N] = multcompare(mp_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(mp_M(:,2), mp_M(:,1), mp_C);
sig_values(mp_P(2), mp_P(1));
title('Ctx-CA1 Modulation')
ylabel('Ctx Amplitude')
xlabel('CA1 Phase')

SLM_ctx_vals = [-3.3 -3.25 -2.96 -3.43 -3.68 -3.02 -3.33 -3.09 -2.27091 -3.05403 -2.75 -3.61 -1.95 -2.83 -2.96 -2.79 -2.38 -3.2 -2.62 -1.63 -2.08 -2.36]';

[ms_P, ms_T, ms_Stats] = anovan(SLM_ctx_vals,{mi_db_labs, mi_age_labs}, 'model', 'interaction','display','off');
[ms_C,ms_M,~,ms_N] = multcompare(ms_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(ms_M(:,2), ms_M(:,1), ms_C);
sig_values(ms_P(2), ms_P(1));
title('Ctx-SLM Modulation')
ylabel('Ctx Amplitude')
xlabel('SLM Phase')

%% TD ratio
td_vals = [0.602 0.31 1.11 0.774 0.969 0.506 0.553 0.531 0.923932 0.628 0.633 0.78 0.392 0.453 0.666 0.395 0.552 0.3 0.292 0.306 0.273 0.553 0.569]';

[td_P, td_T, td_Stats] = anovan(td_vals,{db_labs, age_labs}, 'model', 'interaction','display','off');
[td_C,td_M,~,td_N] = multcompare(td_Stats, 'Dimension', [1 2], 'CType', 'bonferroni','display','off');
figure
create_bar_figure(td_M(:,2), td_M(:,1), td_C);
sig_values(td_P(2), td_P(1));
title('Time Spent in HTD')
ylabel('Proportion')








