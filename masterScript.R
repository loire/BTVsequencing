require(tidyverse)
###### First, import depth
cov = read.csv("Tables/All_depth.tsv",sep="\t",h=F)
colnames(cov) = c("Sample","Chrom","pos","cov")


##### import metadata

meta = read.csv("metadata.csv",sep="\t")
meta %>% select(Sample,AT....)
meta
### Ceiling at 2500
cov$cov = ifelse(cov$cov>2500,2500,cov$cov)
samples = unique(cov$Sample)
pdf("Samples_coverage.pdf",height = 12,width=15)
for(i in samples){
  print(i)
  cov50 = cov %>% filter(Sample==i) %>% 
  group_by(Chrom) %>% mutate(cov50= ifelse(cov>50,TRUE,FALSE)) %>% 
  summarize(cov50=100*round(sum(cov50)/n(),3),meancov = round(mean(cov),2))
  
  p=cov %>% filter(Sample ==i) %>% ggplot() + 
  geom_point(aes(x=pos,y=cov,color= ifelse(cov>50,"black","red")),size=0.05,stat="identity",show.legend = F) + 
  geom_line(aes(x=pos,y=cov),color="grey",size=0.2) + 
  facet_wrap(. ~ Chrom,scales="free_x",nrow=2) + theme_bw() + 
  ggtitle(paste("Coverage analysis of sample ",i)) + xlab("Position") + ylab("Coverage") + 
  geom_text(data = cov50,
            aes(x=500,y=2000,
            label=paste(paste("mean=",meancov),
                        paste("%cov>50=",cov50),sep="\n"))) +
  ylim(0,2500) + 
  scale_color_manual(values=c("black","red"))
  print(p)
  }
dev.off()

###### Detect low coverage samples 
require(ggpubr)
cov %>% group_by(Sample) %>% 
  summarize(meancov = round(mean(cov),2),fractionN50 = round(sum(cov>50)/n(),2)) %>% 
  left_join(meta %>%  select(Sample,species,island), by = c("Sample")) %>% 
  arrange(species,desc(fractionN50)) %>% ggtexttable(rows=NULL)
ggsave("Allsamples_coverage_table.pdf",height = 15,width=10)
cov %>% group_by(Sample) %>% 
  summarize(meancov = round(mean(cov),2),fractionN50 = round(sum(cov>50)/n(),2)) %>% 
  left_join(meta %>%  select(Sample,species,island), by = c("Sample")) %>% 
  arrange(species,desc(fractionN50)) %>% 
  filter(fractionN50>=0.95) %>% 
  ggtexttable(rows=NULL)
ggsave("GoodSamples_coverage_table.pdf",height = 15,width=10)

goodsample = cov %>% group_by(Sample) %>% 
  summarize(meancov = round(mean(cov),2),fractionN50 = round(sum(cov>50)/n(),2)) %>% 
  left_join(meta %>%  select(Sample,species,island), by = c("Sample")) %>% 
  arrange(species,desc(fractionN50)) %>% 
  filter(fractionN50>=0.95) 
goodsample$Sample = droplevels(goodsample$Sample)
goodsample$Sample
####### 

##### Then, import filtered vcf ##### Wait for updated 
fvcf = read.csv("Tables/merged_vcf.tsv",sep="\t",h=T)

datavcf = fvcf %>% mutate(AO = as.character(AO)) %>% 
  separate_rows(sep=",",AO,convert =T) %>% 
  filter(SAMPLE %in% goodsample$Sample) %>% 
  mutate(freq = AO/DP) %>% 
  mutate(id = paste(X.CHROM,POS,sep=":")) %>% 
  left_join(meta %>%  select(Sample,species,island), by = c("SAMPLE"="Sample"))



samplecount = dim(goodsample)[1]

# Search for fixed difference to refs in all samples  
fixed = datavcf %>% filter(freq>0.90) %>%  
  group_by(X.CHROM,POS) %>% 
  summarize(count = n())  %>% 
  filter(count>=samplecount) %>% mutate(id = paste(X.CHROM,POS,sep=":"))

fixed$id %>% length  ### 28 fixed position relative to reference

datavcf$id %>% unique %>% length #### 122 variable among samples positions

datavcf %>% filter(!id %in% fixed$id) %>% filter(DP>100) %>% 
  ggplot() + geom_point(aes(x=POS,y=freq,color=species)) + 
  facet_wrap(X.CHROM ~ .,drop=T, scales = "free_x",nrow=2) + theme_minimal() + ylim(0,1) + 
  theme(panel.border = element_rect(fill=NA))


ggsave("SNP_CHROM_by_species.pdf",height=7,width=10)

ggplot()


dev.off()

chrom = cov %>% group_by(Chrom) %>% 
  summarize(size=max(pos))
heat1 = datavcf %>% filter(!id %in% fixed$id) %>% 
  filter(DP>100) %>% 
  left_join(chrom, by=c("X.CHROM"="Chrom")) %>% 
  group_by(SAMPLE,X.CHROM) %>% summarize(species,SNPcount = as.character(n()),density =round(n()/size*100,2))

heat2 = datavcf %>% filter(!id %in% fixed$id) %>% 
  filter(DP>100) %>% 
  left_join(chrom, by=c("X.CHROM"="Chrom")) %>% 
  group_by(SAMPLE) %>% summarize(X.CHROM="All",species,SNPcount = as.character(n()),density =round(n()/sum(chrom$size)*100,2)) %>% unique



rbind(heat2,heat1) %>% 
  mutate(type = ifelse(X.CHROM=="All","All","Chromosomes")) %>% 
  ggplot() + geom_tile(aes(x=X.CHROM,y=SAMPLE,fill=density)) + 
  geom_text(aes(x=X.CHROM,y=SAMPLE,label=SNPcount)) +
  scale_fill_viridis_c(direction = -1,option="C") + 
  theme_minimal() + theme(panel.grid=element_blank()) + 
  facet_grid(rows=vars(species),cols = vars(type),scales = "free",space="free") + 
  xlab(NULL)
ggsave("HeatMapSNP.pdf",width=8,height=6)


datavcf %>% filter(!id %in% fixed$id) %>% filter(DP>100) %>%
  #mutate(species = fct_recode(species,mammals="mouton",mammals = "vache")) %>% 
  filter(freq> 0.9) %>% 
  group_by(id,species) %>% summarize(REF,ALT,fixedSpecies = n()) -> dataFixedSNP

datavcf %>% select(SAMPLE,species) %>% 
  unique %>% group_by(species) %>% summarize(nb=n()) -> NbSamples

require(stringr)
left_join(dataFixedSNP,NbSamples,by="species") %>% 
  mutate(freq = fixedSpecies/nb) %>% 
  mutate(chrom = as.integer(str_split(id,":")[[1]][1])) %>% 
  mutate(pos = as.integer(str_split(id,":")[[1]][2])) -> dataFreqMut
  
dataFreqMut %>%   
ggplot() + geom_tile(aes(x=species, y=reorder(as.factor(pos), desc(as.factor(pos))),fill=freq)) + 
  geom_text(aes(x=species,y=reorder(as.factor(pos), desc(as.factor(pos))),label=paste(fixedSpecies,nb,sep="/"))) + 
  scale_fill_viridis_c(name = "proportion of samples\n with mutation",option = "D",direction = -1,begin = 0.25)  + 
  facet_grid(rows=vars(chrom),cols = vars(species),scales = "free",space="free") + 
  theme_minimal() + 
  theme(panel.background = element_rect(fill="grey90"), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank(),text = element_text(family="Times",size=14)
  ) + scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + ylab("positions of mutations") + 
  ggtitle("Distribution of fixed mutations among samples")
ggsave("FixedMutationsInSamples.pdf",height = 6,width =7)


######################### Now polymorphics


datavcf %>% filter(!id %in% fixed$id) %>% filter(DP>100) %>%
  #mutate(species = fct_recode(species,mammals="mouton",mammals = "vache")) %>% 
  filter(freq< 0.9) %>% 
  group_by(id,species) %>% summarize(REF,ALT,fixedSpecies = n()) -> dataPolySNP

datavcf %>% select(SAMPLE,species) %>% 
  unique %>% group_by(species) %>% summarize(nb=n()) -> NbSamples

require(stringr)
left_join(dataPolySNP,NbSamples,by="species") %>% 
  mutate(freq = fixedSpecies/nb) %>% 
  mutate(chrom = as.integer(str_split(id,":")[[1]][1])) %>% 
  mutate(pos = as.integer(str_split(id,":")[[1]][2])) -> dataFreqPoly

dataFreqPoly %>%   
  ggplot() + geom_tile(aes(x=species, y=reorder(as.factor(pos), desc(as.factor(pos))),fill=freq)) + 
  geom_text(aes(x=species,y=reorder(as.factor(pos), desc(as.factor(pos))),label=paste(fixedSpecies,nb,sep="/"))) + 
  scale_fill_viridis_c(name = "proportion of samples\n with mutation",option = "D",direction = -1,begin = 0.25)  + 
  facet_grid(rows=vars(chrom),cols = vars(species),scales = "free",space="free") + 
  theme_minimal() + 
  theme(panel.background = element_rect(fill="grey90"), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank(),text = element_text(family="Times",size=14)
  ) + scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + ylab("positions of mutations") + 
  ggtitle("Distribution of polymorphic mutations among samples")
ggsave("PolymorphicMutationsInSamples.pdf",height = 6,width =7)



#################### All mutations ###############
datavcf %>% filter(!id %in% fixed$id) %>% filter(DP>100) %>%
  #mutate(species = fct_recode(species,mammals="mouton",mammals = "vache")) %>% 
  group_by(id,species) %>% summarize(REF,ALT,fixedSpecies = n()) -> dataPolySNP

datavcf %>% select(SAMPLE,species) %>% 
  unique %>% group_by(species) %>% summarize(nb=n()) -> NbSamples

require(stringr)
left_join(dataPolySNP,NbSamples,by="species") %>% 
  mutate(freq = fixedSpecies/nb) %>% 
  mutate(chrom = as.integer(str_split(id,":")[[1]][1])) %>% 
  mutate(pos = as.integer(str_split(id,":")[[1]][2])) -> dataFreqPoly

dataFreqPoly %>%   
  ggplot() + geom_tile(aes(x=species, y=reorder(as.factor(pos), desc(as.factor(pos))),fill=freq)) + 
  geom_text(aes(x=species,y=reorder(as.factor(pos), desc(as.factor(pos))),label=paste(fixedSpecies,nb,sep="/"))) + 
  scale_fill_viridis_c(name = "proportion of samples\n with mutation",option = "D",direction = -1,begin = 0.25)  + 
  facet_grid(rows=vars(chrom),cols = vars(species),scales = "free",space="free") + 
  theme_minimal() + 
  theme(panel.background = element_rect(fill="grey90"), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.grid = element_blank(),text = element_text(family="Times",size=14)
  ) + scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) + ylab("positions of mutations") + 
  ggtitle("Distribution of mutations among samples")
ggsave("AllcMutationsInSamples.pdf",height = 8,width =7)




