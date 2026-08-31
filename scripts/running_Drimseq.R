---
title: "Running_Drimseq"
author: "Selam Mehreteab(smehrete@ucsc.edu/selammeh2004@gmail.com)"
date: "`r format(Sys.time(), '%B %d, %Y')`"
purpose: "Run Drimseq on my Juncbase results"
---

Rscript /private/groups/brookslab/cafelton/git-flair/flair/src/flair/diffExp_drimseq.R --threads 12 --group1 wtdmso --group2 wtcsc --batch b1 --matrix /private/groups/brookslab/smehrete/Juncbase2nd_try/092925.wtdmsovswtcsc.drimformat.mysamples.filtered.tsv --formula /private/groups/brookslab/smehrete/Juncbase2nd_try/wtdmso_wtcsc_manifest_drim.tsv --outDir /drim_output/ --prefix 100725.wtdmsovswtcsc.drimresults

Rscript /private/groups/brookslab/cafelton/git-flair/flair/src/flair/diffExp_drimseq.R --threads 12 --group1 wtdmso --group2 mtdmso --batch b1 --matrix /private/groups/brookslab/smehrete/Juncbase2nd_try/092925.wtdmsovsmtdmso.drimformat.mysamples.filtered.tsv --formula /private/groups/brookslab/smehrete/Juncbase2nd_try/wtdmso_mtdmso_manifest_drim.tsv --outDir /drim_output/ --prefix 100725.wtdmsovsmtdmso.drimresults

script /private/groups/brookslab/cafelton/git-flair/flair/src/flair/diffExp_drimseq.R --threads 12 --group1 wtdmso --group2 mtcsc --batch b1 --matrix /private/groups/brookslab/smehrete/Juncbase2nd_try/092925.wtdmsovsmtcsc.drimformat.mysamples.filtered.tsv --formula /private/groups/brookslab/smehrete/Juncbase2nd_try/wtdmso_mtcsc_manifest_drim.tsv --outDir /drim_output/ --prefix 100725.wtdmsovsmtcsc.drimresults
