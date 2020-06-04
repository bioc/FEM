### GenStatMsp.R

GenStatMsp <- function(dnaM.m,pheno.v,chiptype="450k"){

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
    probeEID.lv <- list()
    for (g in 1:6) {
        group.idx <- which(probeInfo.lv[[3]] == g) 
        tmp.m <- dnaM.m[group.idx, ]
        probes.v <- rownames(tmp.m)
        rownames(tmp.m) <- probeInfo.lv$eid[group.idx]
        sel.idx <- which(is.na(rownames(tmp.m)) == FALSE)
        tmp.m <- tmp.m[sel.idx,]
        beta.lm[[g]] <- tmp.m
        probes.v <- probes.v[sel.idx]
        names(probes.v) <- rownames(tmp.m)
        probeEID.lv[[g]] <- probes.v
        print(paste("Done for regional gene group ", g, sep = ""))
    }

    betaProm.m <- rbind(beta.lm[[2]], beta.lm[[4]])
    betaTS1500.m <- beta.lm[[1]]

    probeProm.v <- c(probeEID.lv[[2]], probeEID.lv[[4]])
    probeTS1500.v <- probeEID.lv[[1]]

    unqEID.v <- sort(unique(c(rownames(betaProm.m), rownames(betaTS1500.m))))

    sampletype.f <- as.factor(pheno.v)
    design.sample <- model.matrix(~0 + sampletype.f)
    colnames(design.sample) <- levels(sampletype.f)
    sampletypes.v <- levels(sampletype.f)

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
	
    lmfProm.o <- lmFit(betaProm.m, design.sample)
    lmfProm2.o <- contrasts.fit(lmfProm.o, cont.m)
    bayProm.o <- eBayes(lmfProm2.o)
	
    lmfTS1500.o <- lmFit(betaTS1500.m, design.sample)
    lmfTS1500_2.o <- contrasts.fit(lmfTS1500.o, cont.m)
    bayTS1500.o <- eBayes(lmfTS1500_2.o)
	
    top.lm <- list()
    probeID.lv <- list()
	
    for (c in 1:ncol(cont.m)) {
        tmp.m <- topTable(bayProm.o, coef = c, adjust.method = "fdr",
            number = nrow(betaProm.m))
        sp.lm <- split(tmp.m, tmp.m[,1])
        tmp1.m  <- matrix(as.numeric(unlist(lapply(sp.lm, FUN=function(x) x[which(abs(x$t) == max(abs(x$t)) ), ] ))), ncol= 7, byrow=T)
        tmp1.idx <- as.numeric(unlist(lapply(sp.lm, FUN=function(x) rownames(x[which(abs(x$t) == max(abs(x$t)) ), ]) )))
        topProbeProm.v <- probeProm.v[tmp1.idx]
        rownames(tmp1.m) <- tmp1.m[,1]
        colnames(tmp1.m) <- colnames(tmp.m)
		
        tmp.m <- topTable(bayTS1500.o, coef = c, adjust.method = "fdr",
            number = nrow(betaTS1500.m))
        sp.lm <- split(tmp.m, tmp.m[,1])
        tmp2.m  <- matrix(as.numeric(unlist(lapply(sp.lm, FUN=function(x) x[which(abs(x$t) == max(abs(x$t)) ), ] ))), ncol= 7, byrow=T)
        tmp2.idx <- as.numeric(unlist(lapply(sp.lm, FUN=function(x) rownames(x[which(abs(x$t) == max(abs(x$t)) ), ]) )))
        topProbeTS1500.v <- probeTS1500.v[tmp2.idx]
        rownames(tmp2.m) <- tmp2.m[,1]
        colnames(tmp2.m) <- colnames(tmp.m)
		
        top.m <- matrix(nrow = length(unqEID.v), ncol = ncol(tmp.m))
        rownames(top.m) <- unqEID.v
        colnames(top.m) <- colnames(tmp.m)
        top.m[match(rownames(tmp2.m), rownames(top.m)),] <- tmp2.m
        top.m[match(rownames(tmp1.m), rownames(top.m)),] <- tmp1.m
        top.m <- top.m[,-1]

        probeID.v <- rep(NA, length(unqEID.v))
        tmp.idx <- match(names(topProbeTS1500.v), unqEID.v)
        probeID.v[tmp.idx] <- topProbeTS1500.v
        names(probeID.v)[tmp.idx] <- names(topProbeTS1500.v)
        tmp.idx <- match(names(topProbeProm.v), unqEID.v)
        probeID.v[tmp.idx] <- topProbeProm.v
        names(probeID.v)[tmp.idx] <- names(topProbeProm.v)

        top.lm[[c]] <- top.m
        probeID.lv[[c]] <- probeID.v
		
    }
	
    return(list(top = top.lm, cont = cont.m, probeID = probeID.lv))

}


