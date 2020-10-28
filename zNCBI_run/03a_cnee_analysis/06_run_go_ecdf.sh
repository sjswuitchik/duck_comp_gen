## in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/postPhyloAcc/go_perms

module load R/4.0.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.0.2

Rscript go_ecdf.R > std.Rout 2> std.Rerr

Rscript analyze_go_ecdf_top1.R > std.Rout 2> std.Rerr
Rscript analyze_go_ecdf_top2.R > std.Rout 2> std.Rerr
Rscript analyze_go_ecdf_top3.R > std.Rout 2> std.Rerr
