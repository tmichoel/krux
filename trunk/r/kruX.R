CalculateKruskalWallisWithMatrix <- function(geno.para, mrna.para, 
                                             p.value.threshold = 1e-3,
                                             num.transcript = 1000,
                                             shuffled.mrna = F){
  # Compute Kruskal-Wallis test p-values for pairs of transcripts and SNPs with
  # matrix operations.
  # Tested on 2.14.1
  # Args:
  #           geno.para: Matrix for SNPs where each column are SNPs of a particular
  #                      sample and each row are the values of a particular SNP 
  #		         cross all samples. The column and row names of the matrix
  #		         represent sample names and SNP ID, respectively.
  #           mrna.para: Matrix for expressiion values of transcrtips. The format 
  #                      is defined same as genotype matrix.
  #   p.value.threshold: Only pairs of SNP-transcript  with p-values less than the threshold 
  #                      are included in the calculation
  #      num.transcript: The number of transcripts to be calculated in one round.
  #			 The parameter is determined by available  memory size  
  #       shuffled.mrna: If T, mrna is shuffled, so the result is based on 
  #                      perturbation data. Default value is F
  # Returns: 
  #   Ordered transcript-SNP pairs by their p-values
  # Date
  #   Jun 4, 2012
  #   V1.0 Modify script to handle missing value in expression values
  #        and genotype values
  #   Oct 2, 2012
  #   V1.1 Fixed bug about index variable
  #        Add parameter to support shuffle mna data

  PreCalculate <-function(SNPs.id, gene.start, gene.end){
    #This function calculates freedome for chi-square test, numbers of SNPs 
    #for no-missing  transcripts, and threshold for the corresponding p-value
    #Note that the number of SNPs with a particular value
    #should be corrected for missing values in transcripts

    genes <- row.names(mrna)[gene.start:gene.end] #Genes precessed 
    mrna.missing <- matrix(data= 0, nrow=length(genes), ncol=dim(mrna)[2]) 
    mrna.missing[is.na(mrna[genes, ])] <- 1

    #Calculate the number of missing transcripts for each SNP with a 
    #value from 0,1,2
    snps.count.correct.0  <- geno.rank.0[SNPs.id, ] %*% t(mrna.missing)
    snps.count.correct.1  <- geno.rank.1[SNPs.id, ] %*% t(mrna.missing)
    snps.count.correct.2  <- geno.rank.2[SNPs.id, ] %*% t(mrna.missing)

    #The degree of freedom of chi-square test between a pair 
    # of SNP-transcript is the number of different 
    #values of the SNP in all samples with non-missing value 
    #of the transcript
    chi.square.fd <- (((snps[SNPs.id, "0"] - snps.count.correct.0) > 0)
                     + ((snps[SNPs.id, "1"] - snps.count.correct.1) > 0) 
                     + ((snps[SNPs.id, "2"] - snps.count.correct.2) > 0)
                     - 1)
    
    #Calculate the total number of effctive samples in a test between 
    #a pair of SNP-transcript
    N <- snps[SNPs.id, "0"] - snps.count.correct.0 +
         snps[SNPs.id, "1"] - snps.count.correct.1 +
         snps[SNPs.id, "2"] - snps.count.correct.2
    ks.test.threshold<-(qchisq(p.value.threshold, chi.square.fd, lower.tail=F) +
                        3 * (N + 1)) * (N * (N + 1)) / 12
    ks.test.threshold[which(chi.square.fd < 1)] <- (2^.Machine$integer.max)

    result <-list(snps.count.correct.0 = snps.count.correct.0, 
                  snps.count.correct.1 = snps.count.correct.1,
                  snps.count.correct.2 = snps.count.correct.2,
                  chi.square.fd = chi.square.fd, N = N, 
                  ks.test.threshold = ks.test.threshold, genes = genes)
    #browser()
    return(result)
  }
  
  CalculateNoMissing <- function(SNP.id, pre.result){
    #This function handles SNPs without missing value. 

    #calculate ranks for selected transcripts
    #NA's rank is replaced by 0
    expression.rank <- apply(mrna[pre.result[["genes"]],,drop=FALSE],
                             1, rank, na.last = 'keep')
    expression.rank[which(is.na(expression.rank))] <- 0

    #calulate rank statistics for different genotype values (0,1,2) 
    #Bug fixed: add drop=FALSE, April 10, 2013
    ks.test <- matrix(0, nrow=dim(geno[SNP.id,,drop=FALSE])[1], 
                      ncol=length(pre.result[["genes"]]))
    rs.0 <- (((geno.rank.0[SNP.id, ] %*% expression.rank) ^ 2) / 
            (snps[SNP.id, "0"] - pre.result[["snps.count.correct.0"]]))
    rs.0[which(is.na(rs.0))] <- 0 #remove NaN by 0/0
    rs.1 <- (((geno.rank.1[SNP.id, ] %*% expression.rank) ^ 2) /
            (snps[SNP.id, "1"] - pre.result[["snps.count.correct.1"]]))
    rs.1[which(is.na(rs.1))] <- 0
    rs.2 <- (((geno.rank.2[SNP.id, ] %*% expression.rank) ^ 2) /
            (snps[SNP.id, "2"] - pre.result[["snps.count.correct.2"]]))
    rs.2[which(is.na(rs.2))] <- 0
    ks.test <- rs.0 + rs.1 + rs.2
    pre.result[["ks.test"]] <- ks.test
    #browser()
    return(pre.result)
  }
  
  CalculateMissing <- function(SNP.id, pre.result){
    #This function handles SNPs with missing value. 
    #It process one gene in each call. 

    #Order SNPs according transcript order (NA is last) 
    #There is only one gene.
    #browser()
    gene <- pre.result[["genes"]][1]
    geno.missing.order <- geno[SNP.id, order(mrna[gene, ]), drop=F]
    geno.rank.0.order <- geno.rank.0[SNP.id, order(mrna[gene, ])]
    geno.rank.1.order <- geno.rank.1[SNP.id, order(mrna[gene, ])]
    geno.rank.2.order <- geno.rank.2[SNP.id ,order(mrna[gene, ])]

    #Calculate correction due to missing value in SNPs and transcripts
    tmp<-matrix(0,nrow=dim(geno.missing.order)[1], 
                  ncol=dim(geno.missing.order)[2])
    tmp[which(is.na(geno.missing.order), arr.ind = T)] <- 1
    cul.sum.missing <- t(apply(tmp, 1, cumsum))
    order.seq <- seq(1,dim(mrna)[2])

    #Since missing value of in the transcript are represented by 
    #zeor in order.seq. The corresponding part in cul.sum.missing
    #is replaced by zeor, too. 
    #Note that NA is in the tail of the output of order function.
    if (T %in% is.na(mrna[gene, ])){
      cul.sum.missing[, (sum(!(is.na(mrna[gene, ]))) + 1):
                        dim(cul.sum.missing)[2]] <- 0
      order.seq[(sum(!(is.na(mrna[gene, ]))) + 1):length(order.seq)] <- 0
    } 

    #calulate rank statistics for different genotype values (0,1,2) 
    #with correction for missing value in SNPs and transcripts.
    rs.0 <- (geno.rank.0.order %*% order.seq - 
             rowSums(geno.rank.0.order * cul.sum.missing))^2 / 
             (snps[SNP.id, "0"] - pre.result[["snps.count.correct.0"]])
    rs.0[which(is.na(rs.0))]<-0 #remove NaN by 0/0
    rs.1 <- (geno.rank.1.order %*% order.seq - 
             rowSums(geno.rank.1.order * cul.sum.missing))^2 /
             (snps[SNP.id, "1"] - pre.result[["snps.count.correct.1"]])
    rs.1[which(is.na(rs.1))] <- 0
    rs.2 <- (geno.rank.2.order %*% order.seq -
             rowSums(geno.rank.2.order * cul.sum.missing))^2 /
             (snps[SNP.id, "2"] - pre.result[["snps.count.correct.2"]])
    rs.2[which(is.na(rs.2))] <- 0
    ks.test <- rs.0 + rs.1 + rs.2
    pre.result[["ks.test"]] <- ks.test

    #browser()
    return(pre.result)
  }

  PostCalculate <- function(SNP.id, pre.result){
    #This function caculates P-vauels for pairs with statistics more than 
    #the threshold and save the result to output list
    snps.selected <- which(pre.result[["ks.test"]] >= 
                           pre.result[["ks.test.threshold"]], 
                           arr.ind = T)

    #Fix bug to add drop=F April 10, 2013
    output.round <- data.frame(snps = row.names(geno[SNP.id,,drop=F])[snps.selected[, 1]],
                               gene = pre.result[["genes"]][snps.selected[, 2]],
                               fd = pre.result[["chi.square.fd"]][snps.selected],
                               chi = pre.result[["ks.test"]][snps.selected] * 
                                     12 / (pre.result[["N"]][snps.selected] * 
                                     (pre.result[["N"]][snps.selected] + 1)) - 
                                     3 * (pre.result[["N"]][snps.selected] + 1),
                               stringsAsFactors = FALSE)
    output.round$pvalue <- pchisq(output.round$chi, output.round$fd, 
                                  lower.tail=F)
    #browser()
    return(output.round)
  }
  
  #Start main program
  #Output list
  t.start = proc.time()
  output <- list()  

  #Make sure expression matrix and SNP matrix with valuse for same samples 
  geno <- geno.para[,sort(intersect(colnames(geno.para), colnames(mrna.para)))]
  mrna <- mrna.para[,sort(intersect(colnames(geno.para), colnames(mrna.para)))]

  #shuffle mrna if shuffled.mrna for FDR calculation
  if (shuffled.mrna == T){
    cat("Expression data is shuffled for FDR calculation\n")
    mrna <- mrna[, sample(seq(1, dim(mrna)[2]))]
    colnames(mrna) <- colnames(geno)
  }

  #seperate SNPs with missing values from the others
  id.missing <- apply(geno, 1, function(x){return(sum(is.na(x)))})


  #Genotype index matrix
  geno.rank.0 <- matrix(0, nrow=dim(geno)[1], ncol=dim(geno)[2])
  geno.rank.0[which(geno==0)] <- 1
  geno.rank.1 <- matrix(0, nrow=dim(geno)[1], ncol=dim(geno)[2])
  geno.rank.1[which(geno==1)] <- 1
  geno.rank.2 <- matrix(0, nrow=dim(geno)[1], ncol=dim(geno)[2])
  geno.rank.2[which(geno==2)] <- 1

  #The number of SNPs with a particular value
  snps <- matrix(data=NA, nrow=dim(geno)[1], ncol=3, 
                        dimnames=list(row.names(geno),c("0","1","2")))
  snps[, "0"] <- rowSums(geno == 0, na.rm = TRUE)
  snps[, "1"] <- rowSums(geno == 1, na.rm = TRUE)
  snps[, "2"] <- rowSums(geno == 2, na.rm = TRUE)

  index <- 0
  #Start to process SNPs without missing value
  if ( TRUE %in% (id.missing == 0)){
    t.1<-proc.time()
    id <- 0 #fix bug at version 1.1
    #In each round, only num.transcript of transcripts are processed. 
    while((id * num.transcript + 1) <= dim(mrna)[1]){
     ind.start <- id * num.transcript + 1
     ind.end <- min((id + 1) * num.transcript, dim(mrna)[1])
     cat("Processing from transcipt ID from ", ind.start, " to ", 
         ind.end," ",format(Sys.time(), "%a %b %d %X %Y"),"\n")
     pre.result <- PreCalculate(id.missing == 0, ind.start, ind.end) 
     pre.result <- CalculateNoMissing(id.missing == 0, pre.result)
     result <- PostCalculate(id.missing == 0, pre.result)
     #Chech if result is empty before adding it to output
     if (dim(result)[1] > 0){
      output[[index + 1]] <- result
      index <- index + 1
     }
     id <- id +1
    }
    t.2<-proc.time()
    #cat("Total running time for SNPS without missing values ", t.2-t.1,"\n")
  }

  #Start to process SNPs with missing value
  #Have process SNPs with missing value against each transcript
  if ( TRUE %in% (id.missing > 0)){
    t.1<-proc.time()
    for (i in seq(1, dim(mrna)[1])){
      #cat("Processing from transcipt ID ", 
      #    i ," ",format(Sys.time(), "%a %b %d %X %Y"),"\n")
      pre.result <- PreCalculate(id.missing > 0, i, i) 
      pre.result <- CalculateMissing(id.missing > 0, pre.result)
      result <- PostCalculate(id.missing > 0, pre.result)
      #Chech if result is empty before adding it to output
      if (dim(result)[1] > 0){
        output[[index + 1]] <- result
        index <- index + 1
      }
    }
    t.2<-proc.time()
    #cat("Total running time for SNPS with missing values ", t.2-t.1,"\n")
  }
  
   #cat("Total running time for R matrix-based test: ", proc.time()-t.start,"\n") 
  return(output)
}  


brutal.force.kruskal <-function(geno.matrix, mrna.matrix, p.value.threshold=1e-3 ){
  #This function calls Kruskal-Wallis. Use it to test matrix-based function
  #Parameters
  # 1. SNPs matrix
  # 2. mrna matrix
  #output
  # 1. data frame for ks test values
  t.start = proc.time()
  test.value <- data.frame()
  for (snp in row.names(geno.matrix)){
    for (gene in row.names(mrna.matrix)){
      missing <- is.na(geno.matrix[snp, ]) | is.na(mrna.matrix[gene, ])
      if (length(table(t(geno.matrix[snp, !missing]))) < 2){
        next();
      }
      k.test <- kruskal.test(t(mrna.matrix[gene, !missing]), 
                           t(geno.matrix[snp, !missing]))
      if (k.test$p.value < p.value.threshold){
      	test.value <-rbind(test.value,
                       data.frame(snps = snp, gene = gene,
                       fd = k.test$parameter, chi = k.test$statistic, 
                       pvalue = k.test$p.value, stringsAsFactors=F))
     }
    }
  }
  test.value <- test.value[order(test.value$snps, test.value$gene), ]
  cat("Total running time for R built-in test: ", proc.time()-t.start,"\n") 
  return(test.value)
}
