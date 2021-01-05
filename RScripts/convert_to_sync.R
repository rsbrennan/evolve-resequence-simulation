#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

option_list = list(
    make_option(c("-n", "--nqtl"), type="character", default=NA, 
              help="number of QTL simulated"),
    make_option(c("-s", "--sln"), type="character", default=NA, 
              help="strength of selection, for example 10, 20, 30 (equal to 0.1, 0.2, 0.3)")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.na(opt$nqtl)){
  cat("ERROR: No argument provided for nqtl\n\n")
  print_help(opt_parser)
  quit(status=1)
}

if (is.na(opt$sln)){
  cat("ERROR: No argument provided for sln\n\n")
  print_help(opt_parser)
  quit(status=1)
}


# read in mutation files

SimReps <- 100

for(i in 1:SimReps){

    setwd(paste0("~/evolve-resequence-simulation/simulations/NQTL", opt$nqtl, "/SimRep", i,"/", opt$sln, "_percent"))
    #setwd("~/evolve-resequence-simulation/simulations/NQTL20/SimRep1/10_percent")
    # read in sim data to a llist
    temp = list.files(pattern="*mutations.txt")
    dfl = lapply(temp, read.csv, header=F, sep=" ")

    # set colnames
    colnm <- c("mut_id", "mut_id2", "mut_type", "position", "S", "dominance", "subpop", "generation", "count")
    dfl <- lapply(dfl, setNames, colnm)

    #get freqs of the mutations
    for(j in 1:length(dfl)){
        dfl[[j]]$frequency <- dfl[[j]]$count/20000
    }

    # generate base counts for sync, 100x coverage

    df2 <- lapply(dfl, function(x) { x$ct1 <- round(100*x$frequency);
                          x$ct2 <- round(100*(1-x$frequency));
                          return(x)})

    df2 <- lapply(df2, function(x) { x[order(x$pos),]})

    ## make sync file:
    # sync file
    ## column 1: reference contig
    ## column 2: position in the reference contig
    ## column 3: refernce character
    ## column >3: allele frequencies for all populations in the form    A-count:T-count:C-count:G-count:N-count:deletion-count

    out <- data.frame( chr = rep(c("chr1"),nrow(df2[[1]])),
                pos = df2[[1]]$position,
                reference = rep(c("A"), nrow(df2[[1]])),
                ctr1 = paste(df2[[1]]$ct1, df2[[1]]$ct2, "0", "0", "0", "0", sep=":"),
                ctr2 = paste(df2[[2]]$ct1, df2[[2]]$ct2, "0", "0", "0", "0", sep=":"),
                ctr3 = paste(df2[[3]]$ct1, df2[[3]]$ct2, "0", "0", "0", "0", sep=":"),
                ctr4 = paste(df2[[4]]$ct1, df2[[4]]$ct2, "0", "0", "0", "0", sep=":"),
                sel1 = paste(df2[[5]]$ct1, df2[[5]]$ct2, "0", "0", "0", "0", sep=":"),
                sel2 = paste(df2[[6]]$ct1, df2[[6]]$ct2, "0", "0", "0", "0", sep=":"),
                sel3 = paste(df2[[7]]$ct1, df2[[7]]$ct2, "0", "0", "0", "0", sep=":"),
                sel4 = paste(df2[[8]]$ct1, df2[[8]]$ct2, "0", "0", "0", "0", sep=":")
                )

    write.table(file=paste0("~/evolve-resequence-simulation/simulations/NQTL", opt$nqtl, "/SimRep", i, "/", opt$sln, "_percent","/SimRep", i, ".sync"), out,
                row.names=F, col.names=F, quote=F, sep="\t")

        print(paste0("Converted file:", "~/evolve-resequence-simulation/simulations/NQTL", opt$nqtl, "/SimRep", i, "/", opt$sln,"_percent/SimRep", i, ".sync"))

}
