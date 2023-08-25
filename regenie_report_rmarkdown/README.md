Run the rendering job in the cluster. This code takes around 5-10 minutes for initial calculations, then ~ 2 min per trait.
To choose what flags/parameters to add, type `Rscript render_summarize_mbGWAS.R -h` for the help page. **NB:** The plotly option is no longer available, it will be removed in the following commits.

```
sbatch -c16 --mem=32G -pslow ~/bin/R/Rscript-4.2-BioC3.15 render_summarize_mbGWAS.R -t 16
```
If you prefer to skip the single-trait report part (partly done by Regenie, too), you can flag -f after calling the rendering 

```
sbatch -c16 --mem=32G -pslow /shared/bioinf/R/bin/Rscript-4.3-BioC3.17 render_summarize_mbGWAS.R -f -t 16
```

My personal preference is not to have interactive plots, to keep the report smaller. To keep it also complete, I choose to leave out the -f|--fast flag. For an eye on the progress, you can use the flag `-v|--verbose`, that will print when it is doing which plot in the per-trait section, at the end

```
sbatch -c16 --mem=32G -pslow /shared/bioinf/R/bin/Rscript-4.3-BioC3.17 render_summarize_mbGWAS.R -t 16 -v
```


 
