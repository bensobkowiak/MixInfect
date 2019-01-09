############################################################################################################
################ ESTIMATION OF MIXTURES AND MIXTURE PROPORTIONS FROM VARIANT CALLING #######################
############################################################################################################

MixInfect<-function(InputVCF,ExcludedRegions.present=FALSE,ExcludedRegions,output,Qual=20,LowCov=10,BICoutput=TRUE){
  
  options(stringsAsFactors = F)
  #### Load required packages
  if (!require(stringr)){
    install.packages("stringr")
    library(stringr)
  }
  if (!require(seqinr)){
    install.packages("seqinr")
    library(seqinr)
  }
  if (!require(Rsamtools)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rsamtools")
    a
    y
    library("Rsamtools")
  }
  if (!require(mclust)){
    install.packages("mclust")
    library("mclust", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
  }
  
  ############################ DATA INPUTS AND QC ######################################
  
  #### Read in vcf file and isolate header
  input<-read.table(file = InputVCF)
  header_input<-as.matrix(read.table(file = InputVCF,comment.char=" ",sep="\n"))
  end_head<-which(grepl("#CHROM",header_input)==TRUE)
  names<-unlist(strsplit(header_input[end_head],"\t"))
  format<-which(names=='FORMAT')
  names<-names[10:length(names)]
  rm(header_input,end_head)
  
  #### Determine fields of GT and DP (or DP4)
  s<-unlist(str_split(input[1,format], ":"))
  GT<-which(s=='GT')
  if (length(which(s=='DP'))==0){
    DP<-which(s=='DP4')
  } else {
    DP<-which(s=='DP')
  }
  rm(format,s)
  
  #### Make new matrices of separared GT:PL:DP:GQ
  if (ncol(input)>10){
  geno<-input[,10:ncol(input)]
  } else {
    geno<-as.data.frame(input[,10:ncol(input)])
  }
  s<-unlist(str_split(geno[1,1], ":"))
  no_of_columns<-length(s)     
  GT1<-matrix(0,nrow(geno),ncol(geno)*no_of_columns)
  l<-seq(1,ncol(GT1),no_of_columns)
  for (i in 1:nrow(geno)){
    for (j in 1:ncol(geno)){
      s<-unlist(str_split(geno[i,j], ":"))
      for (k in 1:no_of_columns){ 
        GT1[i,l[j]+k-1]<-s[k]
      } 
    }
  }
  sep_geno<-cbind(input[,1:9],GT1) 
  rm(geno,i,j,k,l,s,input,GT1)
  
  ##### Remove low qual, >biallelic SNPs and all Indels
  
  #### Remove low quality variants (<q20)
  lowqual<-Qual
  qual<-as.integer(sep_geno[,6])
  b<-qual < lowqual
  c<-which(b == 'FALSE')
  sep_geno<-sep_geno[c,]
  rm(b,c,qual,lowqual)
  
  #### Remove Indels
  desc<-sep_geno[,8]
  new_desc<-grepl('INDEL',desc)
  f<-which(new_desc=='FALSE')
  sep_geno<-sep_geno[f,]
  rm(desc,f,new_desc)
  
  #### Remove > biallelic SNPs
  alt<-sep_geno[,5]
  a<-which(grepl(",",alt)==FALSE)
  sep_geno<-sep_geno[a,]
  rm(a,alt)
  
  #### Mark low coverage (less than 10) as '?'
  chromosome<-as.character(sep_geno[,1])
  position<-sep_geno[,2]
  ref<-as.character(sep_geno[,4])
  alt<-as.character(sep_geno[,5])
  geno<-as.matrix(sep_geno[,seq(GT+9,ncol(sep_geno),no_of_columns)])
  read<-data.matrix(as.matrix(sep_geno[,seq(DP+9,ncol(sep_geno),no_of_columns)]))
  class(read)<-"numeric"
  snprd<-read<LowCov
  lowread<-which(snprd == 'TRUE') 
  geno[lowread]<-'?'
  genotype<-as.matrix(cbind(chromosome,position,ref,alt,geno))
  rm(read,snprd,lowread,sep_geno,alt,chromosome,position,ref,GT,DP,no_of_columns)
  row.names(geno)<-NULL
  row.names(genotype)<-NULL
  
  #### Removal of PE/PPE and known resistance genes (optional)
  if (ExcludedRegions.present==TRUE){
    exclude<-read.csv(file=ExcludedRegions,header=T)
    excluded_sites<-numeric()
    for (i in 1:nrow(exclude)){
      excluded_sites<-c(excluded_sites,exclude[i,1]:exclude[i,2])
    }
    b<-which(is.element(as.numeric(genotype[,2]),excluded_sites)==T)
    if (length(b)!=0){
      geno<-as.data.frame(geno[-b,])
      genotype<-genotype[-b,]
    }
    rm(b,exclude,excluded_sites,i)
  }
  
  
  ######################### ESTIMATION OF MIXED SAMPLES ############################
  
  ### Create output file
  
  outputfile<-as.data.frame(matrix(NA,nrow=length(names),ncol=11))
  colnames(outputfile)<-c("Sample name","Mix or Non-mix","Mixed SNPs","Total SNPs","Proportion het/total SNPs","Number of strains in mix","Major strain proportion","SD Major strain proportion","SEM Major strain proportion","CI low Major strain proportion","CI high Major strain proportion")
  outputfile[,1]<-names
  
  #### No of het SNPs and total and proportions
  mixes<-matrix(0,ncol=ncol(geno),nrow=4)
  for (i in 1:ncol(geno)){
    mixes[1,i]<-length(which(geno[,i]=="0/1"))
  }
  for (i in 1:ncol(geno)){
    mixes[2,i]<-length(which(geno[,i]=="1/1"))
  }
  for (i in 1:ncol(geno)){
    mixes[3,i]<-mixes[1,i]+mixes[2,i]
  }
  for (i in 1:ncol(geno)){
    mixes[4,i]<-(mixes[1,i]/mixes[3,i])*100
  }
  
  outputfile[,3]<-mixes[1,]
  outputfile[,4]<-mixes[3,]
  outputfile[,5]<-mixes[4,]
  outputfile[,2]<-'Non-mix'
  outputfile[,6]<-1
  
  
  #################### ESTIMATE PROPORTIONS OF MIXED SAMPLES (up to 3 mixes in sample) #######################
  
  mixgeno<-as.data.frame(geno[,which(outputfile[,5]>0.5 & outputfile[,3]>10)])
  mixgenotype<-cbind(genotype[,1:4],genotype[,which(outputfile[,5]>0.5 & outputfile[,3]>10)+4])
  mixnames<-names[which(outputfile[,5]>0.5 & outputfile[,3]>10)]
  
  #### Remove any SNPs that do not have a hetero call at any sample
  new22<-matrix(0,nrow(mixgeno),ncol(mixgeno))
  for (i in 1:nrow(mixgeno)){
    for (j in 1:ncol(mixgeno)){
      if (mixgeno[i,j]=="0/1"){
        new22[i,j]='FALSE'}
      else  {new22[i,j]='TRUE'
      }
    }
  }
  d<- lapply(1:nrow(new22), function(i){
    all(as.logical(new22[i,]))
  }
  )
  g<-do.call(rbind,d)
  hetero_present<-which(g=='FALSE')
  genotype<-mixgenotype[hetero_present,]
  colnames(genotype)<-c("chromosome","position","ref","alt",mixnames)
  rm(g,geno,new22,d,hetero_present,i,j,mixgeno,mixgenotype,mixes,names)
  
  
  ##### Make matrix with proportion of ref to alt (0 is ref, 1 alt)
  # First mark all ref and alt with 0 or 1 - mark hetero as '-'
  
  output_prop=matrix(0,nrow(genotype),ncol(genotype))
  rownames(output_prop)=rownames(genotype)
  colnames(output_prop)=colnames(genotype)
  
  for (i in 1:nrow(genotype)){
    for (j in 5:ncol(genotype)){
      if (genotype[i,j]=='0/0'){
        output_prop[i,j]=0}
      else  {if (genotype[i,j]=='1/1'){
        output_prop[i,j]=1
      }
        else  {if (genotype[i,j]=='0/1'){
          output_prop[i,j]='-'
        }
          else {
            output_prop[i,j]='?'}
        }
      }
    }
  }
  
  output_prop[,1:4]<-genotype[,1:4]
  colnames(output_prop)<-c("Chromosome","Position","Ref","Alt",mixnames)
  rm(i,j)
  row.names(output_prop)<-NULL
  
  #### Replace all '-' with proportion based on BAM files
  
  what<-c("pos", "seq","qual")
  bamnames<-as.character(colnames(output_prop))
  poorqual<-c('!','“','#','$','%','&','‘','(',')','*','+',',','–','.')
  
  #loop across all samples
  for (i in 5:ncol(output_prop)){
    bamname<-bamnames[i]
    a<-which(output_prop[,i]=='-')
    
    #loop across all heteroSNPs
    for (j in 1:length(a)){
      which<-RangesList(Chromosome=IRanges(as.numeric(output_prop[a[j],2]),as.numeric(output_prop[a[j],2])))
      param<-ScanBamParam(what=what,which=which)
      reads<-scanBam(bamname,param = param)
      SNP<-as.numeric(output_prop[a[j],2])
      
      sequence<-reads[[1]]$seq
      position<-reads[[1]]$pos
      quality<-reads[[1]]$qual
      
      bases<-NULL
      
      # loop across all reads
      for (m in (1:length(position))){
        z<-unlist(sequence[m])
        qual<-unlist(quality[m])
        if (is.na(position[m])!=TRUE){
          if (((SNP+1)-(position[m])<=length(z))==TRUE){
            if (is.element(as.character(qual[(SNP+1)-(position[m])]), poorqual)==FALSE){
              pos<-as.character(z[(SNP+1)-(position[m])])
              if (pos==output_prop[a[j],3]){
                base<-0
              } else if (pos==output_prop[a[j],4]){
                base<-1
              } else {
                base<-'N'
              }
              if ((base==0 | base==1)==TRUE){
                bases<-c(bases,base)
              }
            }
          }
        }
      }
      output_prop[a[j],i]<-length(which(bases==1))/length(bases)
    }
  }
  rm(a,bamname,base,bases,i,j,m,param,pos,position,reads,sequence,SNP,which,z,qual,quality,poorqual,what)
  
  
  #### Proportions of hetero SNP calls
  proportions_alt<-as.data.frame(output_prop[,5:ncol(output_prop)])
  proportions_alt_pos<-output_prop
  
  ########## Do ref allele proportion
  
  output_prop_ref=matrix(0,nrow(genotype),ncol(genotype))
  rownames(output_prop_ref)=rownames(genotype)
  colnames(output_prop_ref)=colnames(genotype)
  
  for (i in 1:nrow(genotype)){
    for (j in 5:ncol(genotype)){
      if (genotype[i,j]=='0/0'){
        output_prop_ref[i,j]=1}
      else  {if (genotype[i,j]=='1/1'){
        output_prop_ref[i,j]=0
      }
        else  {if (genotype[i,j]=='0/1'){
          output_prop_ref[i,j]='-'
        }
          else {
            output_prop_ref[i,j]='?'}
        }
      }
    }
  }
  
  output_prop_ref[,1:4]<-genotype[,1:4]
  colnames(output_prop_ref)<-c("Chromosome","Position","Ref","Alt",mixnames)
  rm(i,j)
  row.names(output_prop_ref)<-NULL
  
  #### Replace all '-' with proportion based on BAM files
  
  what<-c("pos", "seq","qual")
  bamnames<-as.character(colnames(output_prop_ref))
  poorqual<-c('!','“','#','$','%','&','‘','(',')','*','+',',','–','.')
  
  #loop across all samples
  for (i in 5:ncol(output_prop_ref)){
    bamname<-bamnames[i]
    a<-which(output_prop_ref[,i]=='-')
    
    #loop across all heteroSNPs
    for (j in 1:length(a)){
      which<-RangesList(Chromosome=IRanges(as.numeric(output_prop_ref[a[j],2]),as.numeric(output_prop_ref[a[j],2])))
      param<-ScanBamParam(what=what,which=which)
      reads<-scanBam(bamname,param = param)
      SNP<-as.numeric(output_prop_ref[a[j],2])
      
      sequence<-reads[[1]]$seq
      position<-reads[[1]]$pos
      quality<-reads[[1]]$qual
      
      bases<-NULL
      
      # loop across all reads
      for (m in (1:length(position))){
        z<-unlist(sequence[m])
        qual<-unlist(quality[m])
        if (is.na(position[m])!=TRUE){
          if (((SNP+1)-(position[m])<=length(z))==TRUE){
            if (is.element(as.character(qual[(SNP+1)-(position[m])]), poorqual)==FALSE){
              pos<-as.character(z[(SNP+1)-(position[m])])
              if (pos==output_prop_ref[a[j],3]){
                base<-1
              } else if (pos==output_prop_ref[a[j],4]){
                base<-0
              } else {
                base<-'N'
              }
              if ((base==0 | base==1)==TRUE){
                bases<-c(bases,base)
              }
            }
          }
        }
      }
      output_prop_ref[a[j],i]<-length(which(bases==1))/length(bases)
    }
  }
  rm(a,bamname,base,bases,i,j,m,param,pos,position,reads,sequence,SNP,which,z,qual,quality,poorqual,what)
  
  proportions_ref<-as.data.frame(output_prop_ref[,5:ncol(output_prop_ref)])
  proportions_ref_pos<-output_prop_ref
  
  write.csv(proportions_alt,file=paste(output,"_proportions_alt.csv",sep = ""),row.names = FALSE)
  write.csv(proportions_ref,file=paste(output,"_proportions_ref.csv",sep = ""),row.names = FALSE)
  
  ####### Run Guassian Mclust with 2:6 groups
  
  BICvalues<-as.data.frame(matrix(0,ncol=4,nrow = ncol(proportions_alt)))
  BICvalues[,1]<-mixnames
  colnames(BICvalues)<-c("Sample identifier","G=2","G=4","G=6")
  
  for (i in 1:nrow(BICvalues)){
    b<-c(as.numeric(as.character(proportions_ref[which(proportions_ref[,i]!=0 & proportions_ref[,i]!=1 & proportions_ref[,i]!='?'),i])),as.numeric(as.character(proportions_alt[which(proportions_alt[,i]!=0 & proportions_alt[,i]!=1 & proportions_alt[,i]!='?'),i])))
    if (length(which(is.na(b)))>0){
      d<-b[-which(is.na(b))]
    } else {d<-b}
    a<-mclustBIC(d,G=c(2,4,6))[,2]
    if (length(a)==3){
      BICvalues[i,2:4]<-a
    } else { BICvalues[i,2:4]<-c(a,rep(NA,3-length(a)))
    }
  }

  
  ## if over BiC 20 for G=2 then accept, otherwise if under 20 is G=4 > 10, etc.... 
  
  G_number<-as.data.frame(matrix(0,ncol=2,nrow = nrow(BICvalues)))
  G_number[,1]<-BICvalues[,1]
  colnames(G_number)<-c("Sample identifier","G_number")
  
  for (i in 1:nrow(BICvalues)){
    if (BICvalues[i,2]>=20 & is.na(BICvalues[i,2])==FALSE){
      G_number[i,2]<-2
    } else if (BICvalues[i,3]>=20 & is.na(BICvalues[i,3])==FALSE){
      G_number[i,2]<-4
    } else if (BICvalues[i,4]>=20 & is.na(BICvalues[i,4])==FALSE){
      G_number[i,2]<-6
    } else {
      G_number[i,2]<-1
    }
  }
  
  #####
  
  ## THIS IS WHERE THE COUNTING OF CLUSTER POINTS GOES!!!!
  
  ### G = 2 samples
  
  proportions_alt_G2<-as.data.frame(proportions_alt[,which(G_number[,2]==2)])
  proportions_ref_G2<-as.data.frame(proportions_ref[,which(G_number[,2]==2)])
  colnames(proportions_alt_G2)<-as.character(mixnames[which(G_number[,2]==2)])
  colnames(proportions_ref_G2)<-as.character(mixnames[which(G_number[,2]==2)])
  G2_majority<-as.data.frame(matrix(0,ncol=6,nrow = ncol(proportions_alt_G2)))
  G2_majority[,1]<-colnames(proportions_alt_G2)
  
  if (nrow(G2_majority)>0){
    for (i in 1:ncol(proportions_alt_G2)){
      a<-c(as.numeric(as.character(proportions_ref_G2[which(proportions_ref_G2[,i]!=0 & proportions_ref_G2[,i]!=1 & proportions_ref_G2[,i]!='?'),i])),as.numeric(as.character(proportions_alt_G2[which(proportions_alt_G2[,i]!=0 & proportions_alt_G2[,i]!=1 & proportions_alt_G2[,i]!='?'),i])))
      if (length(which(is.na(a)))>0){
        d<-a[-which(is.na(a))]
      } else {d<-a}
      b<-Mclust(d,G=2)
      e<-b$parameters$mean[order(b$parameters$mean,decreasing = T)]
      high_mean_group<-as.integer(which(b$parameters$mean==e[1]))
      high_mean_group_points<-d[which(b$classification==high_mean_group)]
      G2_majority[i,3]<-sd(high_mean_group_points)
      sem<-sd(high_mean_group_points)/sqrt(length(high_mean_group_points))
      CI<-c(mean(high_mean_group_points)-2*sem,mean(high_mean_group_points)+2*sem)
      G2_majority[i,4]<-sem
      G2_majority[i,5]<-CI[1]
      G2_majority[i,6]<-CI[2]
      G2_majority[i,2]<-e[1]
    }
    
    for (i in 1:nrow(G2_majority)){
      a<-which(outputfile[,1]==G2_majority[i,1])
      outputfile[a,2]<-"Mix"
      outputfile[a,6]<-2
      outputfile[a,7]<-G2_majority[i,2]
      outputfile[a,8]<-G2_majority[i,3]
      outputfile[a,9]<-G2_majority[i,4]
      outputfile[a,10]<-G2_majority[i,5]
      outputfile[a,11]<-G2_majority[i,6]
    }
  }
  
  ## G = 4 samples
  
  proportions_alt_G4<-as.data.frame(proportions_alt[,which(G_number[,2]==4)])
  proportions_ref_G4<-as.data.frame(proportions_ref[,which(G_number[,2]==4)])
  colnames(proportions_alt_G4)<-as.character(mixnames[which(G_number[,2]==4)])
  colnames(proportions_ref_G4)<-as.character(mixnames[which(G_number[,2]==4)])
  G4_majority<-as.data.frame(matrix(0,ncol=6,nrow = ncol(proportions_alt_G4)))
  G4_majority[,1]<-colnames(proportions_alt_G4)
  
  if (nrow(G4_majority)>0){
    for (i in 1:ncol(proportions_alt_G4)){
      a<-c(as.numeric(as.character(proportions_ref_G4[which(proportions_ref_G4[,i]!=0 & proportions_ref_G4[,i]!=1 & proportions_ref_G4[,i]!='?'),i])),as.numeric(as.character(proportions_alt_G4[which(proportions_alt_G4[,i]!=0 & proportions_alt_G4[,i]!=1 & proportions_alt_G4[,i]!='?'),i])))
      if (length(which(is.na(a)))>0){
        d<-a[-which(is.na(a))]
      } else {d<-a}
      b<-Mclust(d,G=4)
      e<-b$parameters$mean[order(b$parameters$mean,decreasing = T)]
      high_mean_group<-as.integer(which(b$parameters$mean==e[1]))
      high_mean_group_points<-d[which(b$classification==high_mean_group)]
      G4_majority[i,3]<-sd(high_mean_group_points)
      sem<-sd(high_mean_group_points)/sqrt(length(high_mean_group_points))
      CI<-c(mean(high_mean_group_points)-2*sem,mean(high_mean_group_points)+2*sem)
      G4_majority[i,4]<-sem
      G4_majority[i,5]<-CI[1]
      G4_majority[i,6]<-CI[2]
      G4_majority[i,2]<-e[1]
    }
    for (i in 1:nrow(G4_majority)){
      a<-which(outputfile[,1]==G4_majority[i,1])
      outputfile[a,2]<-"Mix"
      outputfile[a,6]<-3
      outputfile[a,7]<-G4_majority[i,2]
      outputfile[a,8]<-G4_majority[i,3]
      outputfile[a,9]<-G4_majority[i,4]
      outputfile[a,10]<-G4_majority[i,5]
      outputfile[a,11]<-G4_majority[i,6]
    }
  }
  
  # G = 6 samples
  
  proportions_alt_G6<-as.data.frame(proportions_alt[,which(G_number[,2]==6)])
  proportions_ref_G6<-as.data.frame(proportions_ref[,which(G_number[,2]==6)])
  colnames(proportions_alt_G6)<-as.character(mixnames[which(G_number[,2]==6)])
  colnames(proportions_ref_G6)<-as.character(mixnames[which(G_number[,2]==6)])
  G6_majority<-as.data.frame(matrix(0,ncol=6,nrow = ncol(proportions_alt_G6)))
  G6_majority[,1]<-colnames(proportions_alt_G6)
  
  if (nrow(G6_majority)>0){
    for (i in 1:ncol(proportions_alt_G6)){
      a<-c(as.numeric(as.character(proportions_ref_G6[which(proportions_ref_G6[,i]!=0 & proportions_ref_G6[,i]!=1 & proportions_ref_G6[,i]!='?'),i])),as.numeric(as.character(proportions_alt_G6[which(proportions_alt_G6[,i]!=0 & proportions_alt_G6[,i]!=1 & proportions_alt_G6[,i]!='?'),i])))
      if (length(which(is.na(a)))>0){
        d<-a[-which(is.na(a))]
      } else {d<-a}
      b<-Mclust(d,G=6)
      e<-b$parameters$mean[order(b$parameters$mean,decreasing = T)]
      high_mean_group<-as.integer(which(b$parameters$mean==e[1]))
      high_mean_group_points<-d[which(b$classification==high_mean_group)]
      G6_majority[i,3]<-sd(high_mean_group_points)
      sem<-sd(high_mean_group_points)/sqrt(length(high_mean_group_points))
      CI<-c(mean(high_mean_group_points)-2*sem,mean(high_mean_group_points)+2*sem)
      G6_majority[i,4]<-sem
      G6_majority[i,5]<-CI[1]
      G6_majority[i,6]<-CI[2]
      G6_majority[i,2]<-e[1]
    }
    
    for (i in 1:nrow(G6_majority)){
      a<-which(outputfile[,1]==G6_majority[i,1])
      outputfile[a,2]<-"Mix"
      outputfile[a,6]<-4
      outputfile[a,7]<-G6_majority[i,2]
      outputfile[a,8]<-G6_majority[i,3]
      outputfile[a,9]<-G6_majority[i,4]
      outputfile[a,10]<-G6_majority[i,5]
      outputfile[a,11]<-G6_majority[i,6]
    }
  }
  
  ### Number of genes containing het SNPs
  
  gene_coordinates<-read.csv(file="H37Rv_annotation_per_base.csv.gz",header=T)
  SNPs_genes<-as.data.frame(matrix(NA,nrow=nrow(proportions_alt_pos),ncol=3))
  colnames(SNPs_genes)<-c("SNP_coordinate","Gene","Overlap_gene")
  SNPs_genes[,1]<-proportions_alt_pos[,2]
  
  for (i in 1:nrow(SNPs_genes)){
    SNPs_genes[i,2]<-as.character(gene_coordinates[as.numeric(as.character(SNPs_genes[i,1])),2])
    SNPs_genes[i,3]<-as.character(gene_coordinates[as.numeric(as.character(SNPs_genes[i,1])),3])
    
  }
  
  No_genes<-as.data.frame(matrix(NA,nrow=nrow(outputfile),ncol=3))
  No_genes[,1]<-outputfile[,1]
  columnnames<-colnames(proportions_alt)
  for (i in 1:nrow(No_genes)){
    a<-which(columnnames==No_genes[i,1])
    if (length(a)==0){
      No_genes[i,2]<-NA
    } else {
      b<-which(proportions_alt[,a]!=0 & proportions_alt[,a]!=1 & proportions_alt[,a]!="?")
      genenames<-SNPs_genes[b,2]
      c<-which(genenames=="No_gene")
      No_of_noncoding<-as.numeric(length(c))
      if (length(c)!=0){
        genenames<-genenames[-c]
      }
      if (length(genenames)==0){
        No_genes[i,2]<-No_of_noncoding
      } else {
        number_genes<-as.numeric(length(unique(genenames)))
        No_genes[i,2]<-number_genes+No_of_noncoding
      }
    }
    No_genes[i,3]<-No_genes[i,2]/outputfile[i,3]
  }
  
  columnnames<-colnames(outputfile)
  outputfile<-cbind(outputfile,No_genes[,2:3])
  outputfile[which(outputfile[,5]<0.5),9]<-NA
  colnames(outputfile)<-c(columnnames,"No different genes","Proportion hetSNPs to genes")
  
  #### Write output file
  write.csv(outputfile,file=paste(output,"_mixes_separate.csv",sep = ""),row.names = FALSE)
  
  ######### combine only het SNP props if in same gene
  
  ### combine all SNPs in genes
  
  nongenic<-which(SNPs_genes[,2]=="No_gene")
  genicSNPs<-SNPs_genes[-nongenic,]
  uniqGenes<-genicSNPs[!duplicated(genicSNPs[,2]), ] 
  tocombine<-list()
  for (i in 1:nrow(uniqGenes)){
    tocombine[i]<-list(which(as.character(SNPs_genes[,2])==as.character(uniqGenes[i,2])))
  }
  allSNPs<-c(tocombine,as.list(nongenic))
  
  ## Make all alt >=0.5 and all ref <0.5
  
  proportion_alt_high<-as.matrix(proportions_alt)
  proportion_ref_low<-as.matrix(proportions_ref)
  
  for (i in 1:nrow(proportion_alt_high)){
    for (j in 1:ncol(proportion_alt_high)){
      if (proportion_alt_high[i,j] < 0.5){
        proportion_alt_high[i,j]<-1-as.numeric(proportion_alt_high[i,j])
      }
    }
  }
  
  for (i in 1:nrow(proportion_ref_low)){
    for (j in 1:ncol(proportion_ref_low)){
      if (proportion_ref_low[i,j] >= 0.5){
        proportion_ref_low[i,j]<-1-as.numeric(proportion_alt_high[i,j])
      }
    }
  }
  
  
  
  
  proportions_alt_het_combined<-proportion_alt_high
  proportions_ref_het_combined<-proportion_ref_low
  
  for (i in 1:length(allSNPs)){
    a<-unlist(allSNPs[i])
    for (j in 1:ncol(proportions_alt_het_combined)){
      hets<-which(as.numeric(proportion_alt_high[a,j])!=0 & as.numeric(proportion_alt_high[a,j])!=1 & as.numeric(proportion_alt_high[a,j])!='?')
      if (length(hets)>1){
        proportions_alt_het_combined[a[hets[1]],j]<-sum(as.numeric(proportion_alt_high[a[hets],j]))/length(hets)
        for (k in 2:length(hets)){
          proportions_alt_het_combined[a[hets[k]],j]<-"?"
        }
      }
    }
  }
  proportions_alt_het_combined[is.na(proportions_alt_het_combined)] <- "?"
  
  for (i in 1:length(allSNPs)){
    a<-unlist(allSNPs[i])
    for (j in 1:ncol(proportions_ref_het_combined)){
      hets<-which(as.numeric(proportion_ref_low[a,j])!=0 & as.numeric(proportion_ref_low[a,j])!=1 & as.numeric(proportion_ref_low[a,j])!='?')
      if (length(hets)>1){
        proportions_ref_het_combined[a[hets[1]],j]<-sum(as.numeric(proportion_ref_low[a[hets],j]))/length(hets)
        for (k in 2:length(hets)){
          proportions_ref_het_combined[a[hets[k]],j]<-"?"
        }
      }
    }
  }
  proportions_ref_het_combined[is.na(proportions_ref_het_combined)] <- "?"
  
  ##### Run Gaussian on het combined 
  
  BICvalues_c<-as.data.frame(matrix(0,ncol=4,nrow = ncol(proportions_alt_het_combined)))
  BICvalues_c[,1]<-mixnames
  colnames(BICvalues_c)<-c("Sample identifier","G=2","G=4","G=6")
  
  for (i in 1:nrow(BICvalues_c)){
    b<-c(as.numeric(as.character(proportions_ref_het_combined[which(proportions_ref_het_combined[,i]!=0 & proportions_ref_het_combined[,i]!=1 & proportions_ref_het_combined[,i]!='?'),i])),as.numeric(as.character(proportions_alt_het_combined[which(proportions_alt_het_combined[,i]!=0 & proportions_alt_het_combined[,i]!=1 & proportions_alt_het_combined[,i]!='?'),i])))
    if (length(which(is.na(b)))>0){
      d<-b[-which(is.na(b))]
    } else {d<-b}
    a<-mclustBIC(d,G=c(2,4,6))[,2]
    if (length(a)==3){
      BICvalues_c[i,2:4]<-a
    } else { BICvalues_c[i,2:4]<-c(a,rep(NA,3-length(a)))
    }
  }
  
  ## if over BiC 20 for G=2 then accept, otherwise if under 20 is G=4 > 10, etc.... 
  
  G_number_c<-as.data.frame(matrix(0,ncol=2,nrow = nrow(BICvalues_c)))
  G_number_c[,1]<-BICvalues_c[,1]
  colnames(G_number_c)<-c("Sample identifier","G_number")
  
  for (i in 1:nrow(BICvalues_c)){
    if (BICvalues_c[i,2]>=20 & is.na(BICvalues_c[i,2])==FALSE){
      G_number_c[i,2]<-2
    } else if (BICvalues_c[i,3]>=20 & is.na(BICvalues_c[i,3])==FALSE){
      G_number_c[i,2]<-4
    } else if (BICvalues_c[i,4]>=20 & is.na(BICvalues_c[i,4])==FALSE){
      G_number_c[i,2]<-6
    } else {
      G_number_c[i,2]<-1
    }
  }
  
  outputfile_c<-outputfile
  outputfile_c[,2]<-"Non-mix"
  outputfile_c[,6]<-1
  
  ### G = 2 samples
  
  proportions_alt_G2_c<-as.data.frame(proportions_alt_het_combined[,which(G_number_c[,2]==2)])
  proportions_ref_G2_c<-as.data.frame(proportions_ref_het_combined[,which(G_number_c[,2]==2)])
  colnames(proportions_alt_G2_c)<-as.character(mixnames[which(G_number_c[,2]==2)])
  colnames(proportions_ref_G2_c)<-as.character(mixnames[which(G_number_c[,2]==2)])
  G2_majority_c<-as.data.frame(matrix(0,ncol=2,nrow = ncol(proportions_alt_G2_c)))
  G2_majority_c[,1]<-colnames(proportions_alt_G2_c)
  
  if (nrow(G2_majority_c)>0){
    for (i in 1:ncol(proportions_alt_G2_c)){
      a<-c(as.numeric(as.character(proportions_ref_G2_c[which(proportions_ref_G2_c[,i]!=0 & proportions_ref_G2_c[,i]!=1 & proportions_ref_G2_c[,i]!='?'),i])),as.numeric(as.character(proportions_alt_G2_c[which(proportions_alt_G2_c[,i]!=0 & proportions_alt_G2_c[,i]!=1 & proportions_alt_G2_c[,i]!='?'),i])))
      if (length(which(is.na(a)))>0){
        d<-a[-which(is.na(a))]
      } else {d<-a}
      b<-Mclust(a,G=2)
      e<-b$parameters$mean[order(b$parameters$mean,decreasing = T)]
      G2_majority_c[i,2]<-e[1]
    }
    
    for (i in 1:nrow(G2_majority_c)){
      a<-which(outputfile_c[,1]==G2_majority_c[i,1])
      outputfile_c[a,2]<-"Mix"
      outputfile_c[a,6]<-2
      outputfile_c[a,7]<-G2_majority_c[i,2]
    }
  }
  
  ## G = 4 samples
  
  proportions_alt_G4_c<-as.data.frame(proportions_alt_het_combined[,which(G_number_c[,2]==4)])
  proportions_ref_G4_c<-as.data.frame(proportions_ref_het_combined[,which(G_number_c[,2]==4)])
  colnames(proportions_alt_G4_c)<-as.character(mixnames[which(G_number_c[,2]==4)])
  colnames(proportions_ref_G4_c)<-as.character(mixnames[which(G_number_c[,2]==4)])
  G4_majority_c<-as.data.frame(matrix(0,ncol=2,nrow = ncol(proportions_alt_G4_c)))
  G4_majority_c[,1]<-colnames(proportions_alt_G4_c)
  
  if (nrow(G4_majority_c)>0){
    for (i in 1:ncol(proportions_alt_G4_c)){
      a<-c(as.numeric(as.character(proportions_ref_G4_c[which(proportions_ref_G4_c[,i]!=0 & proportions_ref_G4_c[,i]!=1 & proportions_ref_G4_c[,i]!='?'),i])),as.numeric(as.character(proportions_alt_G4_c[which(proportions_alt_G4_c[,i]!=0 & proportions_alt_G4_c[,i]!=1 & proportions_alt_G4_c[,i]!='?'),i])))
      if (length(which(is.na(a)))>0){
        d<-a[-which(is.na(a))]
      } else {d<-a}
      b<-Mclust(a,G=4)
      e<-b$parameters$mean[order(b$parameters$mean,decreasing = T)]
      G4_majority_c[i,2]<-e[1]
    }
    for (i in 1:nrow(G4_majority_c)){
      a<-which(outputfile_c[,1]==G4_majority_c[i,1])
      outputfile_c[a,2]<-"Mix"
      outputfile_c[a,6]<-3
      outputfile_c[a,7]<-G4_majority_c[i,2]
    }
  }
  
  # G = 6 samples
  
  proportions_alt_G6_c<-as.data.frame(proportions_alt_het_combined[,which(G_number_c[,2]==6)])
  proportions_ref_G6_c<-as.data.frame(proportions_ref_het_combined[,which(G_number_c[,2]==6)])
  colnames(proportions_alt_G6_c)<-as.character(mixnames[which(G_number_c[,2]==6)])
  colnames(proportions_ref_G6_c)<-as.character(mixnames[which(G_number_c[,2]==6)])
  G6_majority_c<-as.data.frame(matrix(0,ncol=2,nrow = ncol(proportions_alt_G6_c)))
  G6_majority_c[,1]<-colnames(proportions_alt_G6_c)
  
  if (nrow(G6_majority_c)>0){
    for (i in 1:ncol(proportions_alt_G6_c)){
      a<-c(as.numeric(as.character(proportions_ref_G6_c[which(proportions_ref_G6_c[,i]!=0 & proportions_ref_G6_c[,i]!=1 & proportions_ref_G6_c[,i]!='?'),i])),as.numeric(as.character(proportions_alt_G6_c[which(proportions_alt_G6_c[,i]!=0 & proportions_alt_G6_c[,i]!=1 & proportions_alt_G6_c[,i]!='?'),i])))
      if (length(which(is.na(a)))>0){
        d<-a[-which(is.na(a))]
      } else {d<-a}
      b<-Mclust(a,G=6)
      e<-b$parameters$mean[order(b$parameters$mean,decreasing = T)]
      G6_majority_c[i,2]<-e[1]
    }
    
    for (i in 1:nrow(G6_majority_c)){
      a<-which(outputfile_c[,1]==G6_majority_c[i,1])
      outputfile_c[a,2]<-"Mix"
      outputfile_c[a,6]<-4
      outputfile_c[a,7]<-G6_majority_c[i,2]
    }
  }
  
  ## Remove false positives through gene combine
  for (i in 1:nrow(outputfile)){
    if (outputfile[i,2]=="Non-mix"){
      outputfile_c[i,2]<-"Non-mix"
      outputfile_c[i,6]<-1
      outputfile_c[i,7]<-NA
    }
  }
  
  ### BIC only for Mixtures
  sepmixes<-as.character(outputfile[which(outputfile[,2]=="Mix"),1])
  combinedmixes<-as.character(outputfile_c[which(outputfile_c[,2]=="Mix"),1])
  a<-which(is.element(BICvalues[,1],sepmixes))
  BICvalues<-BICvalues[a,]
  a<-which(is.element(BICvalues_c[,1],combinedmixes))
  BICvalues_c<-BICvalues_c[a,]
  
  #### Write output files
  if (BICoutput==T){
    write.csv(BICvalues_c,file=paste(output,"_BICvalues_het_combined.csv",sep = ""),row.names = FALSE)
    write.csv(BICvalues,file=paste(output,"_BICvalues_separate.csv",sep = ""),row.names = FALSE)
  }
  write.csv(outputfile_c,file=paste(output,"_output_genes_combined.csv",sep = ""),row.names = FALSE)
}
