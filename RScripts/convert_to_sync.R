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
    #setwd("~/evolve-resequence-simulation/simulations/NQTL20/SimRep1/80_percent")
    # read in sim data to a llist
    temp <-  list.files(pattern="*mutations.txt")
    temp2 <- temp[grep("all_indiv_mutations.txt",temp, invert=T)]
    dfl = lapply(temp2, read.csv, header=F, sep=" ")
    all_muts <- read.csv("all_indiv_mutations.txt", header=F, sep=" ")
    # set colnames
    colnm <- c("mut_id", "mut_id_perm", "mut_type", "position", "S", "dominance", "subpop", "generation", "count")
    colnames(all_muts) <- colnm
    dfl <- lapply(dfl, setNames, colnm)

    #read in genome counts
    genome_counts <- read.csv("genome_number.txt", skip=1, header=F) # skipping number of all indivs

    # get number of mutations that should be present:
    mut_count <- nrow(all_muts)

    #get freqs of the mutations
    for(j in 1:length(dfl)){
        # first check that all the mutations are present:
        ## mut_id_perm is the permanent mutation id that matches across all samples
        ## if it doesn't match, add in the missing ones with a frequency of 0
        if(nrow(dfl[[j]]) != mut_count){
            missing <- subset(all_muts, !(all_muts$mut_id_perm %in% dfl[[j]]$mut_id_perm))
            missing$count <- 0
            dfl[[j]] <- rbind(dfl[[j]], missing)
            dfl[[j]]$frequency <- dfl[[j]]$count/genome_counts$V1[j]
        }else{
            dfl[[j]]$frequency <- dfl[[j]]$count/genome_counts$V1[j]
            }


    }

    # generate base counts for sync, 200x coverage

    df2 <- lapply(dfl, function(x) { x$ct1 <- round(200*x$frequency);
                          x$ct2 <- round(200*(1-x$frequency));
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
