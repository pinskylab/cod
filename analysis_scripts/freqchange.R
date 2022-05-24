library(tidyverse)

Can40_afreqs=read.table("data/Can40Loci_afreqs/plink2.Can40.afreq",header=T)
Can13_afreqs=read.table("data/Can40Loci_afreqs/plink2.Can13.afreq",header=T)

Can40_afreqs$ContFreq=Can13_afreqs$ALT_FREQS
Can40_afreqs$FreqChange=Can40_afreqs$ContFreq-Can40_afreqs$ALT_FREQS

plot(Can40_afreqs$ALT_FREQS,Can40_afreqs$FreqChange,pch=16,col="gray",cex=0.3,xlab="Canada 1940 Frequency",ylab="Relative Change in Frequency")

plot(Can40_afreqs$ALT_FREQS,abs(Can40_afreqs$FreqChange),pch=16,col="gray",cex=0.3,xlab="Canada 1940 Frequency",ylab="Absolute Change in Frequency")

Can40_afreqs$AbsFreqChange=abs(Can40_afreqs$ContFreq-Can40_afreqs$ALT_FREQS)

### change vs frequency

tags <- c("[0-0.1)","[0.1-0.2)", "[0.2-0.3)", "[0.3-0.4)", "[0.4-0.5)", "[0.5-0.6)","[0.6-0.7)", "[0.7-0.8)","[0.8-0.9)", "[0.9-1.0)")

v <- Can40_afreqs %>% select(ALT_FREQS,AbsFreqChange) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    ALT_FREQS < 0.1 ~ tags[1],
    ALT_FREQS >= 0.1 & ALT_FREQS < 0.2 ~ tags[2],
    ALT_FREQS >= 0.2 & ALT_FREQS < 0.3 ~ tags[3],
    ALT_FREQS >= 0.3 & ALT_FREQS < 0.4 ~ tags[4],
    ALT_FREQS >= 0.4 & ALT_FREQS < 0.5 ~ tags[5],
    ALT_FREQS >= 0.5 & ALT_FREQS < 0.6 ~ tags[6],
    ALT_FREQS >= 0.6 & ALT_FREQS < 0.7 ~ tags[7],
    ALT_FREQS >= 0.7 & ALT_FREQS < 0.8 ~ tags[8],
    ALT_FREQS >= 0.8 & ALT_FREQS < 0.9 ~ tags[9],
    ALT_FREQS >= 0.9 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=AbsFreqChange)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="bisque",color="black",alpha=0.3) + 
  labs(x='Canada 1940 Frequency',y='Absolute Frequency Change') +
  guides(color=FALSE) +
  theme_minimal()

#### change vs major frequency

Can40_afreqs$RefFreq=1-Can40_afreqs$ALT_FREQS

Can40_afreqs %>% rowwise() %>% mutate(MajorAF=max(RefFreq,ALT_FREQS)) -> Can40MAF

tags <- c("[0.5-0.55)","[0.55-0.6)", "[0.6-0.65)", "[0.65-0.7)", "[0.7-0.75)", "[0.75-0.8)","[0.8-0.85)", "[0.85-0.9)","[0.9-0.95)", "[0.95-1.0)")

v <- Can40MAF %>% select(MajorAF,AbsFreqChange) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    MajorAF < 0.55 ~ tags[1],
    MajorAF >= 0.55 & MajorAF < 0.6 ~ tags[2],
    MajorAF >= 0.6 & MajorAF < 0.65 ~ tags[3],
    MajorAF >= 0.65 & MajorAF < 0.7 ~ tags[4],
    MajorAF >= 0.7 & MajorAF < 0.75 ~ tags[5],
    MajorAF >= 0.75 & MajorAF < 0.8 ~ tags[6],
    MajorAF >= 0.8 & MajorAF < 0.85 ~ tags[7],
    MajorAF >= 0.85 & MajorAF < 0.9 ~ tags[8],
    MajorAF >= 0.9 & MajorAF < 0.95 ~ tags[9],
    MajorAF >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=AbsFreqChange)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="bisque",color="black",alpha=0.3) + 
  labs(x='Canada 1940 Major Allele Frequency',y='Absolute Frequency Change') +
  guides(color=FALSE) +
  theme_minimal()

#######

Lof07_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof07.afreq",header=T)
Lof11_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof11.afreq",header=T)
Lof14_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof14.afreq",header=T)


Lof07_afreqs$Freq11=Lof11_afreqs$ALT_FREQS
Lof07_afreqs$Freq14=Lof14_afreqs$ALT_FREQS

Lof07_afreqs$FreqChange11=Lof07_afreqs$Freq11-Lof07_afreqs$ALT_FREQS
Lof07_afreqs$FreqChange14=Lof07_afreqs$Freq14-Lof07_afreqs$ALT_FREQS


plot(Lof07_afreqs$ALT_FREQS,Lof07_afreqs$FreqChange11,pch=16,col="gray",cex=0.3,xlab="Norway 1907 Frequency",ylab="Relative Change in Frequency, 2011")
plot(Lof07_afreqs$ALT_FREQS,abs(Lof07_afreqs$FreqChange11),pch=16,col="gray",cex=0.3,xlab="Norway 1907 Frequency",ylab="Absolute Change in Frequency, 2011")

plot(Lof07_afreqs$ALT_FREQS,Lof07_afreqs$FreqChange14,pch=16,col="gray",cex=0.3,xlab="Norway 1907 Frequency",ylab="Relative Change in Frequency, 2014")
plot(Lof07_afreqs$ALT_FREQS,abs(Lof07_afreqs$FreqChange14),pch=16,col="gray",cex=0.3,xlab="Norway 1907 Frequency",ylab="Absolute Change in Frequency, 2014")


#### change vs major frequency

Lof07_afreqs$RefFreq=1-Lof07_afreqs$ALT_FREQS

Lof07_afreqs$AbsFreqChange11=abs(Lof07_afreqs$FreqChange11)

Lof07_afreqs$AbsFreqChange14=abs(Lof07_afreqs$FreqChange14)

Lof07_afreqs %>% rowwise() %>% mutate(MajorAF=max(RefFreq,ALT_FREQS)) -> Lof07MAF

tags <- c("[0.5-0.55)","[0.55-0.6)", "[0.6-0.65)", "[0.65-0.7)", "[0.7-0.75)", "[0.75-0.8)","[0.8-0.85)", "[0.85-0.9)","[0.9-0.95)", "[0.95-1.0)")

v <- Lof07MAF %>% select(MajorAF,AbsFreqChange11) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    MajorAF < 0.55 ~ tags[1],
    MajorAF >= 0.55 & MajorAF < 0.6 ~ tags[2],
    MajorAF >= 0.6 & MajorAF < 0.65 ~ tags[3],
    MajorAF >= 0.65 & MajorAF < 0.7 ~ tags[4],
    MajorAF >= 0.7 & MajorAF < 0.75 ~ tags[5],
    MajorAF >= 0.75 & MajorAF < 0.8 ~ tags[6],
    MajorAF >= 0.8 & MajorAF < 0.85 ~ tags[7],
    MajorAF >= 0.85 & MajorAF < 0.9 ~ tags[8],
    MajorAF >= 0.9 & MajorAF < 0.95 ~ tags[9],
    MajorAF >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=AbsFreqChange11)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="bisque",color="black",alpha=0.3) + 
  labs(x='Norway 1907 Major Allele Frequency',y='Absolute Frequency Change - 2011') +
  guides(color="none") +
  theme_minimal()




v <- Lof07MAF %>% select(MajorAF,AbsFreqChange14) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    MajorAF < 0.55 ~ tags[1],
    MajorAF >= 0.55 & MajorAF < 0.6 ~ tags[2],
    MajorAF >= 0.6 & MajorAF < 0.65 ~ tags[3],
    MajorAF >= 0.65 & MajorAF < 0.7 ~ tags[4],
    MajorAF >= 0.7 & MajorAF < 0.75 ~ tags[5],
    MajorAF >= 0.75 & MajorAF < 0.8 ~ tags[6],
    MajorAF >= 0.8 & MajorAF < 0.85 ~ tags[7],
    MajorAF >= 0.85 & MajorAF < 0.9 ~ tags[8],
    MajorAF >= 0.9 & MajorAF < 0.95 ~ tags[9],
    MajorAF >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=AbsFreqChange14)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="bisque",color="black",alpha=0.3) + 
  labs(x='Norway 1907 Major Allele Frequency',y='Absolute Frequency Change - 2014') +
  guides(color="none") +
  theme_minimal()


#### change vs major frequency - Norway contemps

Lof07_afreqs$FreqChangeContemp=Lof07_afreqs$Freq14-Lof07_afreqs$Freq11

Lof07_afreqs$RefFreq11=1-Lof07_afreqs$Freq11

Lof07_afreqs$AbsFreqChangeContemp=abs(Lof07_afreqs$FreqChangeContemp)

Lof07_afreqs %>% rowwise() %>% mutate(MajorAF11=max(RefFreq11,Freq11)) -> Lof11MAF

tags <- c("[0.5-0.55)","[0.55-0.6)", "[0.6-0.65)", "[0.65-0.7)", "[0.7-0.75)", "[0.75-0.8)","[0.8-0.85)", "[0.85-0.9)","[0.9-0.95)", "[0.95-1.0)")

v <- Lof11MAF %>% select(MajorAF11,AbsFreqChangeContemp) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    MajorAF11 < 0.55 ~ tags[1],
    MajorAF11 >= 0.55 & MajorAF11 < 0.6 ~ tags[2],
    MajorAF11 >= 0.6 & MajorAF11 < 0.65 ~ tags[3],
    MajorAF11 >= 0.65 & MajorAF11 < 0.7 ~ tags[4],
    MajorAF11 >= 0.7 & MajorAF11 < 0.75 ~ tags[5],
    MajorAF11 >= 0.75 & MajorAF11 < 0.8 ~ tags[6],
    MajorAF11 >= 0.8 & MajorAF11 < 0.85 ~ tags[7],
    MajorAF11 >= 0.85 & MajorAF11 < 0.9 ~ tags[8],
    MajorAF11 >= 0.9 & MajorAF11 < 0.95 ~ tags[9],
    MajorAF11 >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=AbsFreqChangeContemp)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="bisque",color="black",alpha=0.3) + 
  labs(x='Norway 2011 Major Allele Frequency',y='Absolute Frequency Change - 2011 to 2014') +
  guides(color="none") +
  theme_minimal()

### change covariance - Canada vs Norways!

Can40_afreqs$FreqChange
Lof07_afreqs$FreqChange11
Lof07_afreqs$FreqChange14

covdf=data.frame(Can40_afreqs$CHROM,Can40_afreqs$FreqChange,Lof07_afreqs$FreqChange11,Lof07_afreqs$FreqChange14)

covdf=na.omit(covdf)

cor_Can40_Lof11=cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)/sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))
cor_Can40_Lof14=cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14)/sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14))

cor_all=(cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)+cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14))/(sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))+sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14)))

LGs=unique(covdf$Can40_afreqs.CHROM)

cordf=data.frame(matrix(ncol=3,nrow=23))
colnames(cordf)=c("LG","cor11","cor14")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  cordf$LG[i]=c
  cordf$cor11[i]=cov(cor_4011_LG$Can40_afreqs.FreqChange,cor_4011_LG$Lof07_afreqs.FreqChange11)/sqrt(var(cor_4011_LG$Can40_afreqs.FreqChange)*var(cor_4011_LG$Lof07_afreqs.FreqChange11))
  cordf$cor14[i]=cov(cor_4011_LG$Can40_afreqs.FreqChange,cor_4011_LG$Lof07_afreqs.FreqChange14)/sqrt(var(cor_4011_LG$Can40_afreqs.FreqChange)*var(cor_4011_LG$Lof07_afreqs.FreqChange14))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.15,0.15),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordf$cor11,col="blue",pch=19)
points(x=seq(1,23)+0.1,y=cordf$cor14,col="green",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Canada-Norway 2011","Canada-Norway 2014"),col=c("blue","green"),pch=c(19,19))

### change covariance - Norway07-11 vs Norway01-14!

cor_Lof11_Lof14=cov(covdf$Lof07_afreqs.FreqChange11,covdf$Lof07_afreqs.FreqChange14)/sqrt(var(covdf$Lof07_afreqs.FreqChange11)*var(covdf$Lof07_afreqs.FreqChange14))

LGs=unique(covdf$Can40_afreqs.CHROM)

cordf1114=data.frame(matrix(ncol=3,nrow=23))
colnames(cordf1114)=c("LG","cor1114")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  cordf1114$LG[i]=c
  cordf1114$cor1114[i]=cov(cor_4011_LG$Lof07_afreqs.FreqChange11,cor_4011_LG$Lof07_afreqs.FreqChange14)/sqrt(var(cor_4011_LG$Lof07_afreqs.FreqChange11)*var(cor_4011_LG$Lof07_afreqs.FreqChange14))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.65,0.65),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordf1114$cor1114,col="purple",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Norway 2011 - Norway 2014"),col=c("purple"),pch=c(19,19))

### change covariance - Norway11-14 vs Canada!

Lof07_afreqs$FreqChange1114=Lof07_afreqs$Freq14-Lof07_afreqs$Freq11

covdf1114=data.frame(Can40_afreqs$CHROM,Can40_afreqs$FreqChange,Lof07_afreqs$FreqChange1114)

covdf1114=na.omit(covdf1114)

cor_Lof1114_Can=cov(covdf1114$Lof07_afreqs.FreqChange1114,covdf1114$Can40_afreqs.FreqChange)/sqrt(var(covdf1114$Lof07_afreqs.FreqChange1114)*var(covdf1114$Can40_afreqs.FreqChange))

LGs=unique(covdf1114$Can40_afreqs.CHROM)

cordfCan1114=data.frame(matrix(ncol=3,nrow=23))
colnames(cordfCan1114)=c("LG","cor1114")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf1114[which(covdf1114$Can40_afreqs.CHROM==c),]
  cordfCan1114$LG[i]=c
  cordfCan1114$cor1114[i]=cov(cor_4011_LG$Lof07_afreqs.FreqChange1114,cor_4011_LG$Can40_afreqs.FreqChange)/sqrt(var(cor_4011_LG$Lof07_afreqs.FreqChange1114)*var(cor_4011_LG$Can40_afreqs.FreqChange))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.65,0.65),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordfCan1114$cor1114,col="purple",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Norway 2011 - Norway 2014"),col=c("purple"),pch=c(19,19))

#### compare ####

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.65,0.65),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23),y=cordf1114$cor1114,col="purple",pch=19)
points(x=seq(1,23)-0.1,y=cordf$cor11,col="blue",pch=19)
points(x=seq(1,23)-0.1,y=cordfCan1114$cor1114,col="orange",pch=19)
points(x=seq(1,23)+0.1,y=cordf$cor14,col="green",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=4,y=-0.2,legend=c("Canada 1940-2013 vs Norway 1907-2011","Canada 1940-2013 vs Norway 1907-2014","Norway 1907-2011 vs Norway 1907-2014","Canada 1940-2013 vs Norway 2011-2014"),col=c("blue","green","purple","orange"),pch=c(19,19,19,19))

##### CV for allele freq change - Canada

Can40_CVcalc=data.frame(matrix(nrow=length(unique(Can40_afreqs$ALT_FREQS)),ncol=5))
names=c("p","invp","var_p","CV","n")

colnames(Can40_CVcalc)=names

Can40_CVcalc$p=sort(unique(Can40_afreqs$ALT_FREQS))

for(i in 1:length(Can40_CVcalc$p)) {
  pframe=na.omit(Can40_afreqs[which(Can40_afreqs$ALT_FREQS==Can40_CVcalc$p[i]),])
  Can40_CVcalc$invp[i]=1-Can40_CVcalc$p[i]
  Can40_CVcalc$var_p[i]=var(pframe$FreqChange)
  Can40_CVcalc$CV[i]=Can40_CVcalc$var_p[i]/(Can40_CVcalc$p[i]*Can40_CVcalc$invp[i])
  Can40_CVcalc$n[i]=nrow(pframe)
}

plot(Can40_CVcalc$p,Can40_CVcalc$CV)

Can40_CVcalc <- transform(Can40_CVcalc,majaf=pmax(Can40_CVcalc$p,Can40_CVcalc$invp))

plot(Can40_CVcalc$majaf,Can40_CVcalc$CV,cex=log10(Can40_CVcalc$n)/2,ylab="Variance in 1940-2013 Frequency Change",xlab="Canada 1940 Major Allele Frequency",col="orange")

##### binned variance vs major frequency - Canada

tags <- c("[0.5-0.55)","[0.55-0.6)", "[0.6-0.65)", "[0.65-0.7)", "[0.7-0.75)", "[0.75-0.8)","[0.8-0.85)", "[0.85-0.9)","[0.9-0.95)", "[0.95-1.0)")

v <- Can40_CVcalc %>% select(majaf,CV) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    majaf < 0.55 ~ tags[1],
    majaf >= 0.55 & majaf < 0.6 ~ tags[2],
    majaf >= 0.6 & majaf < 0.65 ~ tags[3],
    majaf >= 0.65 & majaf < 0.7 ~ tags[4],
    majaf >= 0.7 & majaf < 0.75 ~ tags[5],
    majaf >= 0.75 & majaf < 0.8 ~ tags[6],
    majaf >= 0.8 & majaf < 0.85 ~ tags[7],
    majaf >= 0.85 & majaf < 0.9 ~ tags[8],
    majaf >= 0.9 & majaf < 0.95 ~ tags[9],
    majaf >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=CV)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="orange",color="black",alpha=0.3) + 
  labs(x='Canada 1940 Major Allele Frequency',y='CV of Allele Frequency Change') +
  guides(color=FALSE) +
  theme_minimal()


##### CV for allele freq change - Norway

Lof07_CVcalc=data.frame(matrix(nrow=length(na.omit(unique(Lof07_afreqs$ALT_FREQS))),ncol=8))
names=c("p","invp","var_p11","var_p14","CV11","CV14","n11","n14")

colnames(Lof07_CVcalc)=names

Lof07_CVcalc$p=sort(unique(Lof07_afreqs$ALT_FREQS))

for(i in 1:length(Lof07_CVcalc$p)) {
  pframe=na.omit(Lof07_afreqs[which(Lof07_afreqs$ALT_FREQS==Lof07_CVcalc$p[i]),])
  Lof07_CVcalc$invp[i]=1-Lof07_CVcalc$p[i]
  Lof07_CVcalc$var_p11[i]=var(pframe$FreqChange11)
  Lof07_CVcalc$CV11[i]=Lof07_CVcalc$var_p11[i]/(Lof07_CVcalc$p[i]*Lof07_CVcalc$invp[i])
  Lof07_CVcalc$n11[i]=nrow(pframe)
  Lof07_CVcalc$var_p14[i]=var(pframe$FreqChange14)
  Lof07_CVcalc$CV14[i]=Lof07_CVcalc$var_p14[i]/(Lof07_CVcalc$p[i]*Lof07_CVcalc$invp[i])
  Lof07_CVcalc$n14[i]=nrow(pframe)
}

plot(Lof07_CVcalc$p,Lof07_CVcalc$CV11)

plot(Lof07_CVcalc$p,Lof07_CVcalc$CV14)

Lof07_CVcalc <- transform(Lof07_CVcalc,majaf=pmax(Lof07_CVcalc$p,Lof07_CVcalc$invp))

plot(Lof07_CVcalc$majaf,Lof07_CVcalc$CV11,cex=log10(Lof07_CVcalc$n11)/2,ylab="Variance in 1907-2011 Frequency Change",xlab="Norway 1907 Major Allele Frequency",col="blue")
plot(Lof07_CVcalc$majaf,Lof07_CVcalc$CV14,cex=log10(Lof07_CVcalc$n14)/2,ylab="Variance in 1907-2014 Frequency Change",xlab="Norway 1907 Major Allele Frequency",col="green")

##### binned CV vs major frequency - Norway

tags <- c("[0.5-0.55)","[0.55-0.6)", "[0.6-0.65)", "[0.65-0.7)", "[0.7-0.75)", "[0.75-0.8)","[0.8-0.85)", "[0.85-0.9)","[0.9-0.95)", "[0.95-1.0)")

v <- Lof07_CVcalc %>% select(majaf,CV11) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    majaf < 0.55 ~ tags[1],
    majaf >= 0.55 & majaf < 0.6 ~ tags[2],
    majaf >= 0.6 & majaf < 0.65 ~ tags[3],
    majaf >= 0.65 & majaf < 0.7 ~ tags[4],
    majaf >= 0.7 & majaf < 0.75 ~ tags[5],
    majaf >= 0.75 & majaf < 0.8 ~ tags[6],
    majaf >= 0.8 & majaf < 0.85 ~ tags[7],
    majaf >= 0.85 & majaf < 0.9 ~ tags[8],
    majaf >= 0.9 & majaf < 0.95 ~ tags[9],
    majaf >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=CV11)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="blue",color="black",alpha=0.3) + 
  labs(x='Norway 1907 Major Allele Frequency',y='CV of 1907-2011 Allele Frequency Change') +
  guides(color=FALSE) +
  theme_minimal()

v <- Lof07_CVcalc %>% select(majaf,CV14) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    majaf < 0.55 ~ tags[1],
    majaf >= 0.55 & majaf < 0.6 ~ tags[2],
    majaf >= 0.6 & majaf < 0.65 ~ tags[3],
    majaf >= 0.65 & majaf < 0.7 ~ tags[4],
    majaf >= 0.7 & majaf < 0.75 ~ tags[5],
    majaf >= 0.75 & majaf < 0.8 ~ tags[6],
    majaf >= 0.8 & majaf < 0.85 ~ tags[7],
    majaf >= 0.85 & majaf < 0.9 ~ tags[8],
    majaf >= 0.9 & majaf < 0.95 ~ tags[9],
    majaf >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=CV14)) + 
  geom_jitter(aes(color='blue'),alpha=0.2) +
  geom_boxplot(fill="green",color="black",alpha=0.3) + 
  labs(x='Norway 1907 Major Allele Frequency',y='CV of 1907-2014 Allele Frequency Change') +
  guides(color=FALSE) +
  theme_minimal()


##### CV for allele freq change - Norway Contemporary!

LofContemp_CVcalc=data.frame(matrix(nrow=length(na.omit(unique(Lof07_afreqs$Freq11))),ncol=5))
names=c("p","invp","var","CV","n11")

colnames(LofContemp_CVcalc)=names

LofContemp_CVcalc$p=sort(unique(Lof07_afreqs$Freq11))

for(i in 1:length(LofContemp_CVcalc$p)) {
  pframe=na.omit(Lof07_afreqs[which(Lof07_afreqs$Freq11==LofContemp_CVcalc$p[i]),])
  LofContemp_CVcalc$invp[i]=1-LofContemp_CVcalc$p[i]
  LofContemp_CVcalc$var[i]=var(pframe$FreqChange11)
  LofContemp_CVcalc$CV[i]=LofContemp_CVcalc$var[i]/(LofContemp_CVcalc$p[i]*LofContemp_CVcalc$invp[i])
  LofContemp_CVcalc$n[i]=nrow(pframe)
}

plot(LofContemp_CVcalc$p,LofContemp_CVcalc$CV)

LofContemp_CVcalc <- transform(LofContemp_CVcalc,majaf=pmax(LofContemp_CVcalc$p,LofContemp_CVcalc$invp))

plot(LofContemp_CVcalc$majaf,LofContemp_CVcalc$CV,cex=log10(LofContemp_CVcalc$n)/2,ylab="Variance in 2011-2014 Frequency Change",xlab="Norway 2011 Major Allele Frequency",col="purple")

##### binned CV vs major frequency - Norway

tags <- c("[0.5-0.55)","[0.55-0.6)", "[0.6-0.65)", "[0.65-0.7)", "[0.7-0.75)", "[0.75-0.8)","[0.8-0.85)", "[0.85-0.9)","[0.9-0.95)", "[0.95-1.0)")

v <- LofContemp_CVcalc %>% select(majaf,CV) #pick the variable 
vgroup <- as_tibble(v) %>% 
  mutate(tag = case_when(
    majaf < 0.55 ~ tags[1],
    majaf >= 0.55 & majaf < 0.6 ~ tags[2],
    majaf >= 0.6 & majaf < 0.65 ~ tags[3],
    majaf >= 0.65 & majaf < 0.7 ~ tags[4],
    majaf >= 0.7 & majaf < 0.75 ~ tags[5],
    majaf >= 0.75 & majaf < 0.8 ~ tags[6],
    majaf >= 0.8 & majaf < 0.85 ~ tags[7],
    majaf >= 0.85 & majaf < 0.9 ~ tags[8],
    majaf >= 0.9 & majaf < 0.95 ~ tags[9],
    majaf >= 0.95 ~ tags[10]
  ))
summary(vgroup)

vgroup$tag <- factor(vgroup$tag,
                     levels = tags,
                     ordered = FALSE)
summary(vgroup$tag)


ggplot(data = vgroup, mapping = aes(x=tag,y=CV)) + 
  geom_jitter(aes(color='purple'),alpha=0.2) +
  geom_boxplot(fill="purple",color="black",alpha=0.3) + 
  labs(x='Norway 2011 Major Allele Frequency',y='CV of 2011-2014 Allele Frequency Change') +
  guides(color=FALSE) +
  theme_minimal()


###### plotting all CVs
par(mfrow=c(1,3))

plot(x=NULL,y=NULL,xlim=c(0.49,1.01),ylim=c(0,1.5),xlab="",ylab="CV in Frequency Change, Norway 2011-2014")
points(LofContemp_CVcalc$majaf,LofContemp_CVcalc$CV,cex=log10(LofContemp_CVcalc$n)/2,col="purple")

plot(x=NULL,y=NULL,xlim=c(0.49,1.01),ylim=c(0,1.5),xlab="MajorAF",ylab="CV in Frequency Change, Norway 1907-2011")
points(Lof07_CVcalc$majaf,Lof07_CVcalc$CV11,cex=log10(Lof07_CVcalc$n11)/2,col="blue")

plot(x=NULL,y=NULL,xlim=c(0.49,1.01),ylim=c(0,1.5),xlab="",ylab="CV in Frequency Change, Norway 1907-2014")
points(Lof07_CVcalc$majaf,Lof07_CVcalc$CV14,cex=log10(Lof07_CVcalc$n14)/2,col="green")

