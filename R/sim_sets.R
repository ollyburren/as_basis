

## function to compute relatedness of basis diseases
#this will be computed from 1K genome chromosome 1 blocks removing blocks
#with less than 10 SNPs.

tf<-function(b,i,gs){
    ns<-2^(i-1)
    ss<-ceiling(gs/2^i)
    return (head(split(sample(b), ceiling(seq_along(b)/ss)),n=ns))
}

##sv sample space (integer vector)
##ng number of groups
##gs group size

selectGroups<-function(sv,ng,gs){
    branches<-ceiling(log(ng)/log(2))
    bag<-sv
    po<-list()
    ## create po a list of sets for each branch point that have the correct characteristics
    for (i in 1:(branches+1)){
        # select 2^i set of size casual/2^i 
        l<-tf(bag,i,gs)
        bag<-bag[!bag %in% unlist(l)]
        po[[i]]<-l
    }
    ## create indices for selecting from po to get correct tree structure    
    bs<-function(n) sort(rep(1:(ng/(2^n)),2^n))
    indices<-cbind(1,do.call('cbind',lapply((branches-1):0,bs)))
    ## create samples
    fs<-list()
    for(i in 1:nrow(indices)){
        r<-lapply(seq_along(indices[i,]),function(y){
            unlist(po[[y]][indices[i,y]])
        })
        fs[[i]]<-do.call('c',r)
    }
    ## add random selection to make up to n.causal.blocks
    lo<-sample(bag[!bag %in% unique(unlist(fs))])
    for(i in 1:length(fs)){
        ll<-gs-length(fs[[i]])
        a<-lo[1:ll]
        fs[[i]]<-sort(as.numeric(c(fs[[i]],a)))
        lo<-lo[-(1:ll)]
    }
    return(fs)
}



#n.blocks<-1:length(chr1.r2s)
## select groups of LD blocks such that there are 8 groups (representing a disease) 
#fs<-selectGroups(n.blocks,8,50)

