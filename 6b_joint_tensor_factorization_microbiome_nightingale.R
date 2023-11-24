library(np)

tempted_all <- function(featuretable_one, featuretable_two, timepoint, subjectID,
                        threshold=1, pseudo=1, transform="clr",
                        r = 3, smooth=1e-6,
                        interval = NULL, resolution = 51,
                        maxiter=100, epsilon=1e-6, r_svd=1,
                        do_ratio=TRUE, pct_ratio=0.05, absolute=FALSE,
                        pct_aggregate=1, contrast=NULL){

  datlist <- format_tempted(featuretable=featuretable_one, timepoint=timepoint, subjectID=subjectID,
                            threshold=threshold, pseudo=pseudo, transform=transform)
  datalist_raw <- format_tempted(featuretable=featuretable_one, timepoint=timepoint, subjectID=subjectID,
                                 threshold=threshold, pseudo=pseudo, transform="none")
  mean_svd <- svd_centralize(datlist, r_svd)

  datlist2 <- format_tempted(featuretable=featuretable_two, timepoint=timepoint, subjectID=subjectID,
                             threshold=threshold, pseudo=pseudo, transform=transform)
  datalist_raw2 <- format_tempted(featuretable=featuretable_two, timepoint=timepoint, subjectID=subjectID,
                                  threshold=threshold, pseudo=pseudo, transform="none")
  mean_svd2 <- svd_centralize(datlist2, r_svd)

  res_tempted <- ftsvd_2_v0(datlist=mean_svd$datlist,
                            datlist_2=mean_svd2$datlist,
                            r = r, smooth=smooth,
                            interval = interval, resolution = resolution,
                            maxiter=maxiter, epsilon=epsilon)
  return(res_tempted)
}

svd_centralize <- function(datlist, r = 1){
  n <- length(datlist)
  p <- nrow(datlist[[1]])-1
  mean_hat <- matrix(0,n,p)
  for (i in 1:length(datlist)){
    mean_hat[i,] <- apply(datlist[[i]][-1,], 1, mean)
  }
  mean_hat_svd <- svd(mean_hat, nu=r, nv=r)
  mean_hat_svd1 <- mean_hat_svd$u %*% t(mean_hat_svd$v * mean_hat_svd$d[1:r])
  mf.new <- datlist
  for (i in 1:length(datlist)){
    mf.new[[i]][-1,] <- datlist[[i]][-1,] - mean_hat_svd1[i,]
  }
  results <- list("datlist" = mf.new, "A_tilde" = mean_hat_svd$u,
                  "B_tilde" = mean_hat_svd$v, "lambda_tilde" = mean_hat_svd$d[1:r])
  return(results)
}

##########################
# RKHS functional regression
#########################
bernoulli_kernel <- function(x, y){
  k1.x = x-0.5
  k1.y = y-0.5
  k2.x = 0.5*(k1.x^2-1/12)
  k2.y = 0.5*(k1.y^2-1/12)
  xy = abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4.xy = 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern.xy = k1.x %*% t(k1.y) + k2.x %*% t(k2.y) - k4.xy + 1
  return(kern.xy)
}

freg_rkhs <- function(Lt, Ly, z, interval, resolution = 101, smooth=1e-8){
  tm = NULL
  for (i in 1:length(Lt)) {
    tm = c(tm, Lt[[i]])
  }
  A = matrix(0, length(tm), length(tm))
  A2 = NULL
  c = rep(0, length(tm))
  for(i in 1:length(Lt)){
    Ki = bernoulli_kernel(Lt[[i]], tm)
    A = A + z[i]^2 * t(Ki)%*%Ki
    A2 = rbind(A2, Ki)
    c = c + z[i] * t(Ki) %*% Ly[[i]]
  }
  A.temp = A + smooth * A2
  A.temp.eig = eigen(A.temp, symmetric = TRUE)
  A.d = A.temp.eig$value
  A.d[A.d<1e-10] = 1e-10
  
  beta = solve( (A.temp.eig$vector)%*%(t(A.temp.eig$vector)*A.d) ) %*% c
  
  phi.est = bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm) %*% beta
  #plot(phi.est)
  return(phi.est)
}

format_tempted <- function(featuretable, timepoint, subjectID,
                           threshold=1, pseudo=NULL, transform="clr"){
  ntm <- which(table(subjectID)==1)
  if(length(ntm)>0)
    stop(paste('Please remove these subjects with only one time point:',
               paste(names(ntm), collapse=', ')))
  if (length(subjectID)!=nrow(featuretable))
    stop('length of subjectID does not match featuretable!')
  if (length(timepoint)!=nrow(featuretable))
    stop('length of timepoint does not match featuretable!')
  # get pseudo count
  if (is.null(pseudo) & (transform %in% c("clr", "logcomp", "logit"))){
    pseudo <- apply(featuretable, 1, function(x){
      min(x[x!=0])/2
    })
  }
  # keep taxon that has non-zeros in >1-threshold samples
  featuretable <- featuretable[,colMeans(featuretable==0)<=threshold]
  if(transform=='logcomp'){
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable)))
  }else if(transform=='comp'){
    featuretable <- featuretable
    featuretable <- t(featuretable/rowSums(featuretable))
  }else if(transform=='ast'){
    featuretable <- featuretable
    featuretable <- t(asin(sqrt(featuretable/rowSums(featuretable))))
  }else if(transform=='clr'){
    featuretable <- featuretable+pseudo
    featuretable <- log(featuretable/rowSums(featuretable))
    featuretable <- t(featuretable-rowMeans(featuretable))
  }else if(transform=='logit'){
    featuretable <- featuretable+pseudo
    featuretable <- t(featuretable/rowSums(featuretable))
    featuretable <- log(featuretable/(1-featuretable))
  }else if(transform=='none'){
    featuretable <- t(featuretable)
  }else{
    print('Input transformation method is wrong! logcomp is applied instead')
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable)))
  }
  featuretable <- rbind(timepoint, featuretable)
  rownames(featuretable)[1] <- 'timepoint'
  subID <- unique(subjectID)
  nsub <- length(subID)

  # construct list of data matrices, each element representing one subject
  datlist <- vector("list", length = nsub)
  names(datlist) <- subID

  # Each slice represents an individual (unequal sized matrix).
  for (i in 1:nsub){
    # print(i)
    datlist[[i]] <- featuretable[, subjectID==subID[i]]
    datlist[[i]] <- datlist[[i]][,order(datlist[[i]][1,])]
    datlist[[i]] <- datlist[[i]][,!duplicated(datlist[[i]][1,])]
  }
  return(datlist)
}

ftsvd_2_v0 <- function(datlist, datlist_2, interval = NULL, r = 3, resolution = 251, smooth=1e-8,
                    maxiter=20, epsilon=1e-4){
  
  n = length(datlist) #42
  p = nrow(datlist[[1]])-1 # number of features in data 1
  p_2 = nrow(datlist_2[[1]])-1#ADDTION # number of features in data 2
  
  Lambda = rep(0, r) # Scalar
  Lambda_2 = rep(0, r)#ADDTION # Scalar
  A = matrix(0, n, r) # Subject singular vectors across datasets
  B = matrix(0, p, r) # Feature singular vectors data 1
  B_2 = matrix(0, p_2, r)#ADDTION # Feature singular vectors data 2
  
  Phi = matrix(0, resolution, r) # Time singular functions data 1
  Phi_2 = matrix(0, resolution, r)#ADDTION # Time singular functions data 2
  
  PCname <- paste('Component', 1:r)
  colnames(A) = PCname
  colnames(B) = PCname
  colnames(B_2) = PCname#ADDTION
  colnames(Phi) = PCname
  colnames(Phi_2) = PCname#ADDTION
  
  rownames(A) = names(datlist)
  rownames(B) = rownames(datlist[[1]])[-1]
  rownames(B_2) = rownames(datlist_2[[1]])[-1]#ADDTION
  
  
  # Calculate range.
  timestamps.all = NULL
  for (i in 1:n){
    timestamps.all = c(timestamps.all, datlist[[i]][1,])
  }
  
  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }
  
  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
    datlist_2[[i]][1,] = (datlist_2[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])#ADDTION
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  
  res = NULL
  Lambda = rep(0, r)
  X = NULL
  y = NULL
  
  #Second modality
  Lambda_2 = rep(0, r)#ADDTION
  X_2 = NULL#ADDTION
  y_2 = NULL#ADDTION
  
  ti <- vector(mode='list', length=n)
  for (i in 1:n){
    temp = 1 + round((resolution-1) * (datlist[[i]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    ti[[i]] <- temp
  }
  
  tipos <- vector(mode='list', length=n)
  for (i in 1:n){
    keep <- ti[[i]]>0
    y <- c(y, as.vector(t(datlist[[i]][-1,keep])))
    y_2 <- c(y_2, as.vector(t(datlist_2[[i]][-1,keep])))#ADDTION
    tipos[[i]] <- keep
  }
  
  for (s in 1:r){ 
    # calculate rank-1 component sequentially.
    # Step 1: initialization.
    print(sprintf("Calculate the %dth Component", s))
    
    # b.hat = b.initials[,s]
    # obtain intialization of b
    data.unfold = NULL
    data.unfold_2 = NULL#ADDTION
    for (i in 1:n){
      data.unfold = cbind(data.unfold, datlist[[i]][2:(p+1),])
      data.unfold_2 = cbind(data.unfold_2, datlist_2[[i]][2:(p_2+1),])#ADDTION
    }
    b.initials <- svd(data.unfold, nu=r, nv=r)$u
    b.hat = b.initials[,1]
    
    # obtain intialization of b2
    b_2.initials <- svd(data.unfold_2, nu=r, nv=r)$u#ADDTION
    b_2.hat = b_2.initials[,1]#ADDTION
    
    # compress the data and apply function PCA for Phi.
    Lt = list()
    Ly = list()
    for (i in 1:n){
      Lt = c(Lt, list(datlist[[i]][1,]))
      Ly = c(Ly, list(as.numeric(b.hat%*%datlist[[i]][2:(p+1),])))
    }
    phi.hat = freg_rkhs(Lt,Ly,rep(1,length(Lt))/sqrt(length(Lt)),
                        interval=interval, 
                        resolution = resolution, smooth=smooth)
    phi.hat = phi.hat / sqrt(sum(phi.hat^2))
    
    
    # compress the data and apply function PCA for Phi_2.#ADDTION
    Lt_2 = list() # change 1
    Ly_2 = list() # change 2
    for (i in 1:n){
      Lt_2 = c(Lt_2, list(datlist_2[[i]][1,])) # change 3
      Ly_2 = c(Ly_2, list(as.numeric(b_2.hat%*%datlist_2[[i]][2:(p_2+1),]))) # change 3, 4, 5
    }
    phi_2.hat = freg_rkhs(Lt_2,Ly_2,rep(1,length(Lt_2))/sqrt(length(Lt_2)),
                          interval=interval, 
                          resolution = resolution, smooth=smooth) # change 6, 7
    phi_2.hat = phi_2.hat / sqrt(sum(phi_2.hat^2)) 

    ###Liat - to check why Lambda.hat and Lambda_2.hat are initialized this way
    Lambda.hat = 4000
    Lambda_2.hat = 4000
    
    # iteratively update a,b,b2,phi, phi2
    a.hat <- rep(1,n)
    t <- 0
    dif <- 1
    while(t<=maxiter & dif>epsilon){
      # update a:
      a.tilde <- rep(0,n)
      z_1 = rep(0,n)
      z_2 = rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        
        z_1[i] <- b.hat %*% datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]] %*% Lambda.hat
        
        z_2[i] <- b_2.hat %*% datlist_2[[i]][2:(p_2+1),t.temp] %*% phi_2.hat[ti[[i]][t.temp]] %*% Lambda_2.hat
        
        a.tilde[i] <-  (z_1[i] +  z_2[i])/(sum((phi.hat[ti[[i]][t.temp]])^2) %*% Lambda.hat^2 + sum((phi_2.hat[ti[[i]][t.temp]])^2) %*% Lambda_2.hat^2)

      }
            
      a.new <- a.tilde / sqrt(sum(a.tilde^2)) #not needed here but keeping anyway
      dif <- sum((a.hat - a.new)^2) #let's verify this is zero
      a.hat <- a.new
      
      num_tmp = c()
      denom_tmp = c()
      
      num_tmp_2 = c()
      denom_tmp_2 = c()
      
      for(i in 1:n){
        t.temp <- tipos[[i]]
        
        num_tmp[i] =  b.hat %*% datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]] 
        denom_tmp[i] = sum((phi.hat[ti[[i]][t.temp]])^2) 
        
        num_tmp_2[i] =  b_2.hat %*% datlist_2[[i]][2:(p_2+1),t.temp] %*% phi_2.hat[ti[[i]][t.temp]] 
        denom_tmp_2[i] = sum((phi_2.hat[ti[[i]][t.temp]])^2) 
      }
      
      Lambda.hat <- as.numeric(num_tmp%*%(as.numeric(a.hat))) / as.numeric(denom_tmp%*%(as.numeric(a.hat)^2))
      Lambda_2.hat <- as.numeric(num_tmp_2%*%(as.numeric(a.hat))) / as.numeric(denom_tmp_2%*%(as.numeric(a.hat)^2))
      # update b
      temp.num <- matrix(0,p,n)
      temp.denom <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        temp.num[,i] <- datlist[[i]][2:(p+1),t.temp] %*% phi.hat[ti[[i]][t.temp]] 
        temp.denom[i] <-sum((phi.hat[ti[[i]][t.temp]])^2)
      }
      
      ##Original
      b.tilde <- as.numeric(temp.num%*%(as.numeric(a.hat))) / as.numeric(temp.denom%*%(as.numeric(a.hat)^2))
      b.new <- b.tilde / sqrt(sum(b.tilde^2))
      dif <- max(dif, sum((b.hat - b.new)^2))
      b.hat <- b.new
      # update b2
      temp_2.num <- matrix(0,p_2,n)
      temp_2.denom <- rep(0,n)
      for (i in 1:n){
        t.temp <- tipos[[i]]
        temp_2.num[,i] <- datlist_2[[i]][2:(p_2+1),t.temp] %*% phi_2.hat[ti[[i]][t.temp]]
        temp_2.denom[i] <- sum((phi_2.hat[ti[[i]][t.temp]])^2)
      }
      
      ##Original
      b_2.tilde <- as.numeric(temp_2.num%*%(as.numeric(a.hat))) / as.numeric(temp.denom%*%(as.numeric(a.hat)^2))
      b_2.new <- b_2.tilde / sqrt(sum(b_2.tilde^2))
      dif <- max(dif, sum((b_2.hat - b_2.new)^2))
      b_2.hat <- b_2.new
      # updata phi:
      Lt = list()
      Ly = list()
      for (i in 1:n){
        Lt = c(Lt, list(datlist[[i]][1,]))
        Ly = c(Ly, list(as.numeric(b.hat%*%datlist[[i]][2:(p+1),])))
      }
      
      # updata phi:
      ##Original
      phi.hat = freg_rkhs(Lt,Ly,a.hat,interval=interval,
                          resolution = resolution, smooth=smooth)
      phi.hat = phi.hat / sqrt(sum(phi.hat^2))
      # updata phi2:
      phi_2.hat = freg_rkhs(Lt_2,Ly_2,a.hat,interval=interval,
                            resolution = resolution, smooth=smooth)
      phi_2.hat = phi_2.hat / sqrt(sum(phi_2.hat^2))
      
      t <- t+1
    }

    A[,s] = a.hat
    #P.A = P.A - a.hat %*% t(a.hat)
    B[,s] = b.hat
    B_2[,s] = b_2.hat
    Phi[,s] = t(phi.hat)
    Phi_2[,s] = t(phi_2.hat)
    Lambda[s] = Lambda.hat
    Lambda_2[s] = Lambda_2.hat
    
    # update datlist
    for (i in 1:n){
      temp <- tipos[[i]]
      datlist[[i]][2:(p+1),which(temp)] = datlist[[i]][2:(p+1),which(temp)] - 
        Lambda[s] * A[i,s] * (B[,s] %*% t(Phi[ti[[i]][temp],s])) 
      
      datlist_2[[i]][2:(p_2+1),which(temp)] = datlist_2[[i]][2:(p_2+1),which(temp)] - 
        Lambda_2[s] * A[i,s] * (B_2[,s] %*% t(Phi_2[ti[[i]][temp],s]))#ADDITION 
    }
    print(paste0("Convergence reached at dif=", dif, ', iter=', t))
  }

  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]
  results = list("A.hat" = A, "B.hat" = B, "B_2.hat" = B_2,
                 "Phi.hat" = Phi,"Phi_2.hat" = Phi_2, "time" = time.return,
                 "Lambda" = Lambda,"Lambda_2" = Lambda_2)
  return(results)
}

tdenoise_2 <- function(res_ftsvd, mean_svd=NULL){
  n <- nrow(res_ftsvd$A.hat)
  p <- nrow(res_ftsvd$B_2.hat)
  resol <- nrow(res_ftsvd$Phi.hat)
  tensor.est <- array(0,dim=c(n,p,resol))
  if (!is.null(mean_svd))
    tensor.est <- (mean_svd$A.tilde %*% t(mean_svd$B.tilde * mean_svd$lambda.tilde)) %o%
    rep(1, resol)
  for (i in 1:ncol(res_ftsvd$A.hat)){
    tensor.est <- tensor.est+res_ftsvd$A.hat[,i]%o%res_ftsvd$B_2.hat[,i]%o%res_ftsvd$Phi_2.hat[,i]*res_ftsvd$Lambda_2[i]
  }
  dimnames(tensor.est)[[3]] <- res_ftsvd$time
  return(tensor.est)
}

# load and match tables
count_table_one <- read.csv("../data/matched-data/without-csf/table-nightingale-serum-joint-tensors.tsv", sep = '\t', row.names = 1,  header = TRUE)
count_table_one <- t(count_table_one)
count_table_two <- read.csv("../data/matched-data/without-csf/table-microbiome-fecal-joint-tensors.tsv", sep = '\t', row.names = 1, header = TRUE)
count_table_two <- t(count_table_two)
rownames(count_table_one) <- sub("^X", "", rownames(count_table_one))
rownames(count_table_two) <- sub("^X", "", rownames(count_table_two))
meta_table <- read.csv("../data/matched-data/without-csf/metadata-joint-tensors.tsv", sep = '\t', row.names = 1, header = TRUE)
samples_order <- rownames(meta_table)
count_table_one <- count_table_one[samples_order, ]
count_table_two <- count_table_two[samples_order, ]
meta_table <- meta_table[samples_order, ]
# run joint-tempted
tempted_res <- tempted_all(count_table_one,
                           count_table_two,
                           meta_table$timepoint_encoded,
                           meta_table$SubjectID,
                           r=15,
                           resolution=8,
                           smooth=1e-5,
                           pseudo=1,
                           transform="clr")
# get results to save
write.csv(tempted_res$A.hat, "../results/joint-tensor-factor-microbiome-nightingale/a_hat.csv")
write.csv(tempted_res$B.hat, "../results/joint-tensor-factor-microbiome-nightingale/b_hat.csv")
write.csv(tempted_res$B_2.hat, "../results/joint-tensor-factor-microbiome-nightingale/b_two_hat.csv")
write.csv(tempted_res$Phi.hat,  "../results/joint-tensor-factor-microbiome-nightingale/phi_hat.csv")
write.csv(tempted_res$Phi_2.hat,  "../results/joint-tensor-factor-microbiome-nightingale/phi_two_hat.csv")
write.csv(tempted_res$time,  "../results/joint-tensor-factor-microbiome-nightingale/time.csv")
write.csv(tempted_res$Lambda,  "../results/joint-tensor-factor-microbiome-nightingale/lambda.csv")
write.csv(tempted_res$Lambda_2,  "../results/joint-tensor-factor-microbiome-nightingale/lambda_two.csv" )
