### GenStatM.R

GenStatM <- function(dnaM.m,pheno.v,chiptype="450k"){

    if (chiptype == "450k"){
        data("probe450kfemanno")
        probefemanno <- probe450kfemanno
    }
    else if (chiptype == "EPIC" ){
        data("probeEPICfemanno")
        probefemanno <- probeEPICfemanno
    }
    else{
        print("ERROR: Please indicate correct data type!")
        break
    }

    extractFn <- function(tmp.v, ext.idx) {
        return(tmp.v[ext.idx])
    }
    map.idx <- match(rownames(dnaM.m), probefemanno$probeID);
    probeInfo.lv <- lapply(probefemanno, extractFn, map.idx)
    beta.lm <- list()
    for (g in 1:6) {
        group.idx <- which(probeInfo.lv[[3]] == g)
        tmp.m <- dnaM.m[group.idx, ]
        rownames(tmp.m) <- probeInfo.lv$eid[group.idx];
        sel.idx <- which(is.na(rownames(tmp.m)) == FALSE);
        tmp.m <- tmp.m[sel.idx,];
        nL <- length(factor(rownames(tmp.m)));
        nspg.v <- summary(factor(rownames(tmp.m)),maxsum=nL);
        beta.lm[[g]] <- rowsum(tmp.m,group=rownames(tmp.m))/nspg.v;
        print(paste("Done for regional gene group ", g, sep = ""))
    }
    unqEID.v <- unique(c(rownames(beta.lm[[2]]), rownames(beta.lm[[4]]), 
        rownames(beta.lm[[1]])))
    avbeta.m <- matrix(nrow = length(unqEID.v), ncol = ncol(dnaM.m))
    colnames(avbeta.m) <- colnames(dnaM.m)
    rownames(avbeta.m) <- unqEID.v
    for (gr in c(1, 4, 2)) {
        avbeta.m[match(rownames(beta.lm[[gr]]), rownames(avbeta.m)), 
            ] <- beta.lm[[gr]]
    }
    data.m <- avbeta.m
    sampletype.f <- as.factor(pheno.v)
    design.sample <- model.matrix(~0 + sampletype.f)
    colnames(design.sample) <- levels(sampletype.f)
    sampletypes.v <- levels(sampletype.f)
    lmf.o <- lmFit(data.m, design.sample)
    ntypes <- length(levels(sampletype.f))
    ncomp <- 0.5 * ntypes * (ntypes - 1)
    cont.m <- matrix(0, nrow = ncol(design.sample), ncol = ncomp)
    tmp.v <- vector()
    c <- 1
    for (i1 in 1:(ntypes - 1)) {
        for (i2 in (i1 + 1):ntypes) {
            cont.m[i1, c] <- -1
            cont.m[i2, c] <- 1
            tmp.v[c] <- paste(sampletypes.v[i2], "--", sampletypes.v[i1], 
                sep = "")
            c <- c + 1
        }
    }
    rownames(cont.m) <- sampletypes.v
    colnames(cont.m) <- tmp.v
    lmf2.o <- contrasts.fit(lmf.o, cont.m)
    bay.o <- eBayes(lmf2.o)
    top.lm <- list()
    for (c in 1:ncol(cont.m)) {
        top.lm[[c]] <- topTable(bay.o, coef = c, adjust.method = "fdr", 
            number = nrow(data.m))
    }
    
    return(list(top = top.lm, cont = cont.m, avbeta = avbeta.m))

}
