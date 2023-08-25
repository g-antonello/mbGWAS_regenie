0. This folder should be copied automatically by my setup, but in case, copy it manually into your folder
```
git clone https://gitlab.gm.eurac.edu/gantonello/regenie_report_rmarkdown
```

if it asks for **username** and **password**, put your gitlab login credentials

It is still important so far that the destination folder is "report_markdown". After this, the rest is quite simple

1. Move to this newly generated directory
	
		cd path/to/regenie/project/report_markdown/

2. Check that the main settings are good inside the launching script. in case they are not, edit them with nano/vim/notepad/whatever.exe

```
cat render_summarize_mbGWAS.R
```
Eg:
	* You MUST make sure the working directory is correct
	* you CAN chance the p value threshold (default is set by Michele's regenie "[...].config" file in the conf directory

3. Run the rendering job in the cluster. these setting takes ~40 min to complete for 87 traits

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


 
