# Step 5: Perform GO Enrichment Analysis using g:Profiler

After generating the gene list for GO analysis in Step 4, you can use g:Profiler to identify enriched biological functions and pathways.

---

## Input File

Use the gene list output from Step 4:

```
results/04_<prefix>_shinygo_gene_list.txt
```

It should look like this
```
AT1G12780
AT1G63180
AT5G17420
AT4G23920
AT1G08260
AT5G44480
AT5G21060
```


Open the enrichment tool
## g:Profiler Web Server

https://biit.cs.ut.ee/gprofiler/gost

---

## Instructions

1. Open the g:Profiler web interface in your browser.
2. Paste or upload the gene list file from Step 4.
3. Set Organism to: Arabidopsis thaliana
4. Adjust optional parameters:
   - Multiple testing correction method (e.g., g:SCS, Bonferroni)
   - Significance threshold (e.g., p-value < 0.05)
5. Click "Run Query"

---

## Export Results
Examples output
https://biit.cs.ut.ee/gplink/l/a2Uyjp24XTu

You can export the results as:

- Enrichment table (CSV)
- Publication-ready plots (PDF or PNG)

---

## Automation Tip

If needed, g:Profiler also supports programmatic access through:

- R package: gprofiler2
- Python client: gprofiler-official
