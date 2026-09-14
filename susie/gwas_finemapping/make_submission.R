egfr_list <- readRDS('egfr_list_1m.rds')
loci <- names(egfr_list)

dir.create('outputs')

for (i in loci) {
    sample_folder <- paste0('outputs/', i)
    dir.create(sample_folder)

    # write submission.sh
    submission_file <- paste0(sample_folder, "/submission.sh")
    write("#! /bin/bash", submission_file,append=FALSE)
    write("#SBATCH --partition=standard", submission_file,append=TRUE)
    write("#SBATCH --ntasks=1", submission_file,append=TRUE)
    write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
    write("#SBATCH --mem=64G", submission_file,append=TRUE)
    write(paste0("#SBATCH --job-name=",i), submission_file,append=TRUE)
    write(paste0("#SBATCH --output=",i,".%j.out"), submission_file,append=TRUE)
    write(paste0("#SBATCH --error=",i,".%j.err"), submission_file,append=TRUE)
    write("source ~/.bashrc", submission_file, append=TRUE)
    write(sprintf("cd /path/to/susie/gwas_finemapping"), submission_file, append=TRUE) # dir containing run_susie.R

    cmd1 <- sprintf("Rscript run_susie.R %s", i)    
    write(cmd1, submission_file, append=TRUE)
    system(paste0('sbatch ', submission_file)) # submit the job
}

#####################
# summarize results #
#####################
loci <- list.files('outputs/')

cols <- rownames(read.csv('outputs/rs10055349/meta.csv', row.names=1))
meta <- do.call(rbind, lapply(loci, function(x) {
    tmp <- read.csv(sprintf('outputs/%s/meta.csv', x), row.names=1)[,1]
    return(tmp)
}))
rownames(meta) <- loci
colnames(meta) <- cols
meta <- data.frame(meta)
write.csv(meta, 'susie_meta.csv')