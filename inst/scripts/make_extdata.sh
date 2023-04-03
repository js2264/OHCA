cooler dump -t chroms /home/rsg/.cache/R/fourDNData/4d434d8538a0_4DNFI9FVHJZQ.mcool::/resolutions/250000 | grep -P 'chr17\t' > inst/extdata/chromsizes
cooler dump -t pixels -r chr17 --join /home/rsg/.cache/R/fourDNData/4d434d8538a0_4DNFI9FVHJZQ.mcool::/resolutions/5000 > inst/extdata/chr17_5000.bg2
cooler cload pairs inst/extdata/chromsizes:5000 inst/extdata/chr17_5000.bg2 -c1 1 -p1 2 -c2 4 -p2 5 -0 inst/extdata/chr17.cool
cooler zoomify --resolutions 5000,250000 --balance inst/extdata/chr17.cool -o inst/extdata/chr17.mcool
