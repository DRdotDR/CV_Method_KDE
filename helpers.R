library(compiler)

genunifp <- function(n, M) {
  p <- matrix(stats::rexp(n * M), nrow = n,dimnames=list(1:n,LETTERS[1:M]))
  p <- t(apply(p, 1, function(x) x/sum(x)))
}
genunifp <- cmpfun(genunifp)

genormixt <- function(p, mean, sd) {
  nc <- 1:ncol(p)
  obj <- apply(p, 1, function(x) sample(nc, 1, prob = x))
  smp <- stats::rnorm(nrow(p), mean[obj], sd[obj])
}
genormixt <- cmpfun(genormixt)

gengamixt <- function(p, shape, rate) {
  nc <- 1:ncol(p)
  obj <- apply(p, 1, function(x) sample(nc, 1, prob = x))
  smp <- stats::rgamma(nrow(p), shape[obj], rate[obj])
}
gengamixt <- cmpfun(gengamixt)

lsweight <- function(p) {
  p <- as.matrix(p)
  a <- as.data.frame(p %*% MASS::ginv(t(p) %*% p))
  names(a)<-colnames(p)
  a
}
lsweight <- cmpfun(lsweight)

wtsamp <- function(x, cumm = NULL, indiv = NULL) {
  o <- order(x)
  xo <- c(-Inf, x[o])
  if (!is.null(cumm))
    cumm <- rbind(rep(0, ncol(cumm)), cumm[o, ])
  if (!is.null(indiv))
    indiv <- indiv[o, ]
  structure(list(xo = xo, cumm = cumm, indiv = indiv), class = "wtsamp")
}
wtsamp <- cmpfun(wtsamp)

wtcens <- function(x, delta, cumm = NULL, indiv = NULL) {
  o <- order(x)
  xo <- c(-Inf, x[o])
  deltao <- delta[o]
  if (!is.null(cumm))
    cumm <- rbind(rep(0, ncol(cumm)), cumm[o, ])
  if (!is.null(indiv))
    indiv <- indiv[o, ]
  structure(list(xo = xo, deltao = deltao, cumm = cumm, indiv = indiv), class = "wtcens")
}
wtcens <- cmpfun(wtcens)

indiv2cumm <- function(xs) {
  xs$cumm <- apply(xs$indiv, 2, cumsum)
  xs$cumm <- rbind(rep(0, ncol(xs$cumm)), xs$cumm)
  xs
}
indiv2cumm <- cmpfun(indiv2cumm)

cumm2indiv <- function(xs) {
  xs$indiv <- apply(xs$cumm, 2, diff)
  xs
}
cumm2indiv <- cmpfun(cumm2indiv)

edfgen <- function(xs, m) {
  if (is.null(xs$cumm))
    xs <- indiv2cumm(xs)
  a <- xs$cumm[, m]
  x <- xs$xo
  f <- function(t) {
    a[findInterval(t, x)]
  }
}
edfgen <- cmpfun(edfgen)

corsup <- function(xs) {
  if (is.null(xs$cumm))
    xs <- indiv2cumm(xs)
  xs$cumm <- pmin(apply(xs$cumm, 2, cummax), 1)
  xs$indiv <- NULL
  xs
}
corsup <- cmpfun(corsup)

corinf <- function(xs) {
  if (is.null(xs$cumm))
    xs <- indiv2cumm(xs)
  n <- length(xs$x) - 1
  xs$cumm <- xs$cumm[-1, ]
  xs$cumm <- rbind(rep(0, ncol(xs$cumm)), pmax(apply(xs$cumm[n:1, ], 2, cummin), 0)[n:1,
  ])
  xs$indiv <- NULL
  xs
}
corinf <- cmpfun(corinf)

cormid <- function(xs) {
  if (is.null(xs$cumm))
    xs <- indiv2cumm(xs)
  cummin <- pmin(apply(xs$cumm, 2, cummax), 1)
  n <- length(xs$x) - 1
  xs$cumm <- xs$cumm[-1, ]
  xs$cumm <- (rbind(rep(0, ncol(xs$cumm)), pmax(apply(xs$cumm[n:1, ], 2, cummin), 0)[n:1,
  ]) + cummin)/2
  xs$indiv <- NULL
  xs
}
cormid <- cmpfun(cormid)

gensampmixt <- function(p, xs) {
  randwt <- randwtgen(xs)
  crange <- 1:ncol(p)
  components <- apply(p, 1, function(pr) sample(crange, 1, replace = TRUE, prob = pr))
  sapply(components, randwt, n = 1)
}
gensampmixt <- cmpfun(gensampmixt)

randwtgen <- function(xs) {
  x <- xs$x[-1]
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  prob <- xs$indiv
  delta0 = min(diff(x))/2
  randwt <- function(m, n, delta = delta0) {
    r <- sample(x, n, prob = prob[, m], replace = TRUE)
    if (delta > 0)
      r <- r + stats::runif(n, -delta, delta)
    r
  }
}
randwtgen <- cmpfun(randwtgen)

meanw <- function(xs) {
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  mx<-as.vector(xs$xo[-1] %*% as.matrix(xs$indiv))
  names(mx)<-colnames(xs$indiv)
  return(mx)
}
meanw <- cmpfun(meanw)

varw <- function(xs){
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  sx<-as.vector((xs$xo[-1])^2 %*% as.matrix(xs$indiv)) -
    (meanw(xs))^2
  names(sx)<-colnames(xs$indiv)
  return(sx)
}
varw <- cmpfun(varw)

sdw <- function(xs,corr=cormid){
  if(is.null(corr))
  {if (is.null(xs$indiv)) xs <- cumm2indiv(xs)}
  else
  {xs <- cumm2indiv(corr(xs))}
  sx<-sqrt(as.vector((xs$xo[-1])^2 %*% as.matrix(xs$indiv)) -
             (meanw(xs))^2)
  names(sx)<-colnames(xs$indiv)
  return(sx)
}
sdw <- cmpfun(sdw)

quantilew <- function(xs,prob){
  if (is.null(xs$cumm))
    xs <- indiv2cumm(xs)
  n <- nrow(xs$cumm)
  M <- ncol(xs$cumm)
  q <- numeric(M)
  for( m in 1:M ){
    j <- 1
    while(j<n & xs$cumm[j,m]<prob) j <- j+1
    q_left <- xs$xo[j]
    j <- n
    while( j>1 & xs$cumm[j,m]>prob ) j <- j-1
    #    q_right <- xs$xo[j]
    q_right <- xs$xo[min(j+1,n)]
    q[m] <- ( q_left + q_right )/2
  }
  names(q)<-colnames(xs$cumm)
  return(q)
}
quantilew <- cmpfun(quantilew)

medianw <- function(xs)quantilew(xs,0.5)
medianw <- cmpfun(medianw)

IQRw <- function(xs){
  quantilew(xs,0.75) - quantilew(xs,0.25)
}
IQRw <- cmpfun(IQRw)

Epanechn<-function(x)ifelse(abs(x)<1,0.75*(1-x^2),0)
Epanechn <- cmpfun(Epanechn)

densgen<-function(xs,m,Kern=Epanechn){
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  a <- xs$indiv[, m]
  x <- xs$xo[-1]
  f <- Vectorize(
    function(t,h) {
      sum(a*Kern((t-x)/h))/h
    },
    vectorize.args ="t"
  )
}
densgen <- cmpfun(densgen)

silvbw<-function(xs,m,delta=1.7188){
  if (is.null(xs$indiv))
    xs <- cumm2indiv(xs)
  const<-1.011354 # (8*sqrt(pi)/3)^(1/5)/(stats::qnorm(0.75)-stats::qnorm(0.25))
  const*delta*(sum((xs$indiv[,m])^2))^(1/5)*min(IQRw(xs)[m],sdw(xs)[m])
}
silvbw <- cmpfun(silvbw)

sdMean<-function(x,p,comp=1:ncol(p),
                 means=FALSE,CI=FALSE,alpha=0.05){
  M<-ncol(p)
  D<-rep(NA,M)
  if(is.vector(x)&is.numeric(x))
    sx<-wtsamp(x,indiv=lsweight(p))
  else{
    if(class(x)=="wtsamp")
      sx<-x
    else{
      warning("x must be vector or wtsamp")
      return(NA)
    }
  }
  m<-meanw(sx)
  m2<-varw(sx)+m^2
  for(k in comp){
    app<-matrix(ncol=M,nrow=M)
    for(i in 1:M){
      for(l in 1:i){
        app[i,l]<-sum(sx$indiv[,k]^2*p[,i]*p[,l])
        app[l,i]<-app[i,l]
      }
    }
    ap<-apply(app,1,sum)
    D[k]=sum(ap*m2)-m%*%app%*%m
    if(D[k]<0){
      warning("Negative estimate of variance is obtained",D[k])
      D[k]=NA
    }
  }
  D=sqrt(D)
  if(!(means|CI)){
    names(D)<-colnames(p)
    return(D)
  }else{
    R<-data.frame(sd=D)
    if(means)R$means<-m
    if(CI){
      lambda=stats::qnorm(1-alpha/2)
      R$lower<-m-lambda*D
      R$upper<-m+lambda*D
    }
    row.names(R)<-colnames(p)
    return(R)
  }
}
sdMean <- cmpfun(sdMean)

sdMedian<-function(x,p,comp=1:ncol(p),
                   medians=FALSE,CI=FALSE,alpha=0.05){
  M<-ncol(p)
  D<-rep(NA,M)
  if(is.vector(x)&is.numeric(x)){
    sx<-wtsamp(x,indiv=lsweight(p))
    sx <- indiv2cumm(sx)
  }
  else{
    if(class(x)=="wtsamp")
      sx<-x
    if(is.null(sx$cumm)) sx <- indiv2cumm(sx)
    else{
      warning("x must be vector or wtsamp")
      return(NA)
    }
  }
  med<-medianw(sx)
  for(k in comp){
    Fm<-sx$cumm[findInterval(med[k], sx$xo),]
    app<-matrix(ncol=M,nrow=M)
    for(i in 1:M){
      for(l in 1:i){
        app[i,l]<-sum(sx$indiv[,k]^2*p[,i]*p[,l])
        app[l,i]<-app[i,l]
      }
    }
    ap<-apply(app,1,sum)
    f<-densgen(sx,k)
    zz<-sum(ap*Fm)-Fm%*%app%*%Fm
    if(zz<0){
      warning("Negative estimate of variance is obtained",zz)
      zz=NA
    }
    else
      D[k]<-sqrt(zz)/f(med[k],silvbw(sx,k))
  }
  if(!(medians|CI)){
    names(D)<-colnames(p)
    return(D)
  }else{
    R<-data.frame(sd=D)
    if(medians)R$medians<-med
    if(CI){
      lambda=stats::qnorm(1-alpha/2)
      R$lower<-med-lambda*D
      R$upper<-med+lambda*D
    }
    row.names(R)<-colnames(p)
    return(R)
  }
}
sdMedian <- cmpfun(sdMedian)

gencensg <- function(p, shape, rate, shapec, ratec) {
  nc <- 1:ncol(p)
  obj <- apply(p, 1, function(x) sample(nc, 1, prob = x))
  x <- stats::rgamma(nrow(p), shape = shape[obj], rate = rate[obj])
  c <- stats::rgamma(nrow(p), shape = shapec[obj], rate = ratec[obj])
  delta <- x < c
  x <- pmin(x, c)
  cs <- data.frame(x, delta)
}
gencensg <- cmpfun(gencensg)

KMcdf <- function(cs) {
  if (is.null(cs$cumm))
    cs <- indiv2cumm(cs)
  if (is.null(cs$indiv))
    cs <- cumm2indiv(cs)
  wt <- 1 - apply(1 - cs$indiv * cs$deltao/(1 - utils::head(cs$cumm, n = -1)), 2, cumprod)
  wt[is.nan(wt)] <- 1
  wtsamp(cs$xo[-1], cumm = wt, indiv = NULL)
}
KMcdf <- cmpfun(KMcdf)

get_r_density_function <- function(type, param1, param2) {
  switch(
    type,
    "normal" = function()
      rnorm(1, mean = param1, sd = param2),
    "uniform" = function()
      runif(
        1,
        min = min(param1, param2),
        max = max(param1, param2)
      ),
    "cauchy" = function()
      rcauchy(1, location = param1, scale = param2),
    "multi_normal" = function(){
      if (runif(1, 0, 1) > 0.5) {
        rnorm(1, -param1, param2)
      } else{
        rnorm(1, param1, param2)
      }
    }
  )
}

get_density_function <- function(type, param1, param2) {
  switch(
    type,
    "normal" = function(x)
      dnorm(x, mean = param1, sd = param2),
    "uniform" = function(x)
      dunif(
        x,
        min = min(param1, param2),
        max = max(param1, param2)
      ),
    "cauchy" = function(x)
      dcauchy(x, location = param1, scale = param2),
    "multi_normal" = function(x){
      dnorm(x, mean = -param1, sd = param2)/2 +
        dnorm(x, mean = param1, sd = param2)/2
    }
  )
}

  
mixture_density <- function(x,
                            type1,
                            param1_1,
                            param1_2,
                            type2,
                            param2_1,
                            param2_2,
                            w1 = 0.5) {
  dens1 <- get_density_function(type1, param1_1, param1_2)
  dens2 <- get_density_function(type2, param2_1, param2_2)
  w1 * dens1(x) + (1 - w1) * dens2(x)
}

get_x_range <- function(type, param1, param2, margin = 1) {
  if (type == "uniform") {
    p1 <- min(param1, param2)
    p2 <- max(param1, param2)
    c(p1 - margin, p2 + margin)
  } else {
    # For normal and cauchy, use param1 as center and param2 as scale
    c(param1 - 3 * param2 - margin, param1 + 3 * param2 + margin)
  }
}

silvrule <- function(X, a, delta = 0.7764) {
  xs <- wtsamp(X, indiv = a)
  
  const <- 1.364296
  S <- sdw(xs)
  I <- IQRw(xs)
  
  h <- numeric(ncol(xs$indiv))
  
  for (m in 1:ncol(xs$indiv)) {
    h[m] <- const * delta *
      (sum((xs$indiv[, m])^2))^(1/5) *
      min(0.7413011 * I[m], S[m])
  }
  
  h
}

K <- function(x) {
  dnorm(x, mean = 0, sd = 1)
}

gen_plot_1 <- function(res, second=FALSE) {
  L <- res$L
  
  set.seed(42)
  Z_plot <- rep(0, L)
  W1 <- (1:L) / L
  W2 <- 1 - W1
  
  for (j in 1:L)
  {
    if (runif(1, 0, 1) > W1[j]) {
      Z_plot[j] <- get_r_density_function(res$p2$type, res$p2$param1, res$p2$param2)()
    } else{
      Z_plot[j] <- get_r_density_function(res$p1$type,res$p1$param1,  res$p1$param2)()
    }
  }

  V <- matrix(c(mean(W1^2), mean(W1 * W2), mean(W1 * W2), mean(W2^2)), nrow = 2)
  V1 <- ginv(V)
  A1_plot <- V1[1,1] * W1 + V1[1,2] * W2
  A2_plot <- V1[2,1] * W1 + V1[2,2] * W2

  Phat1_plot <- function(x, hN) {
    S <- 0
    for (j in 1:L) {
      S <- S + A1_plot[j] / L / hN * K((x - Z_plot[j]) / hN)
    }
    return(S)
  }
  
  Phat2_plot <- function(x, hN) {
    S <- 0
    for (j in 1:L) {
      S <- S + A2_plot[j] / L / hN * K((x - Z_plot[j]) / hN)
    }
    return(S)
  }
  
  if (second) {
    dens <- get_density_function(res$p2$type, res$p2$param1, res$p2$param2)
    Min <- mean(sapply(res$results, `[[`, "Min2"))
    Sil <- mean(sapply(res$results, `[[`, "Sil2"))
    Phat_cv <- function(x) {
      Phat2_plot(x, Min)
    }
    Phat_sil <- function(x) {
      Phat2_plot(x, Sil)
    }
  } else {
    dens <- get_density_function(res$p1$type, res$p1$param1, res$p1$param2)
    Min <- mean(sapply(res$results, `[[`, "Min1"))
    Sil <- mean(sapply(res$results, `[[`, "Sil1"))
    Phat_cv <- function(x) {
      Phat1_plot(x, Min)
    }
    Phat_sil <- function(x) {
      Phat1_plot(x, Sil)
    }
  }

  ggplot(data = data.frame(x = c(-5, 5)), aes(x = x)) +
    stat_function(
      aes(linetype = "Справжня щільність"),
      fun = dens,
      linewidth = 1.3,
      color = "black"
    ) +
    stat_function(
      aes(linetype = "Оцінка методом крос-валідації"),
      fun = Phat_cv,
      linewidth = 1.2,
      color = "red"
    ) +
    stat_function(
      aes(linetype = "Оцінка методом Сільвермана"),
      fun = Phat_sil,
      linewidth = 1.2,
      color = "blue"
    ) +
    scale_linetype_manual(
      values = c(
        "Справжня щільність" = "solid",
        "Оцінка методом крос-валідації" = "dashed",
        "Оцінка методом Сільвермана" = "dotdash"
      )
    ) +
    labs(
      title = "Компонента: графіки щільності та її оцінок",
      x = "x",
      y = "Щільність",
      linetype = "Метод"
    ) +
    theme_excel() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
}

run_iteration <- function(iter, type1, param1_1, param1_2, type2, param2_1, param2_2, L, T, sigma) {
  library(Rcpp)
  library(functions)
  
  Z<-rep(0,L)
  
  ######Concentrations of the components W
  W1 <- (1:L) / L
  W2 <- 1 - W1
  
  V <- matrix(0, nrow = 2, ncol = 2)
  
  V[1,1] <- mean(W1^2)
  V[2,2] <- mean(W2^2)
  V[1,2] <- mean(W1 * W2)
  V[2,1] <- V[1,2]
  
  V1 <- ginv(V)
  
  A1 <- V1[1,1] * W1 + V1[1,2] * W2
  A2 <- V1[2,1] * W1 + V1[2,2] * W2
  
  dens1 <- get_density_function(type1, param1_1, param1_2)
  dens2 <- get_density_function(type2, param2_1, param2_2)
  
  #Create mixture of 2 components
  for (j in 1:L)
  {
    if (runif(1, 0, 1) > W1[j]) {
      Z[j] <- get_r_density_function(type2, param2_1, param2_2)()
    } else{
      Z[j] <- get_r_density_function(type1, param1_1, param1_2)()
    }
  }
  
  ################ Distance matrix ###########################################
  # Замість подвійного цикла
  DZ <- outer(Z, Z, "-")
  DZ2 <- DZ^2
  
  ################ CV #########################################################
  CV1 <- function(hN) {
    CV1_cpp(hN, DZ, DZ2, W1, W2, V, A1, A2, sigma)
  }
  
  CV2 <- function(hN) {
    CV2_cpp(hN, DZ, DZ2, W1, W2, V, A1, A2, sigma)
  }
  
  ################ Optimization ##############################################
  O1 <- optimize(CV1, c(0, 15))
  O2 <- optimize(CV2, c(0, 15))
  
  Min1 <- O1$minimum
  Min2 <- O2$minimum
  
  ################ KDE estimators ############################################
  Phat1 <- function(x, hN) {
    S <- 0
    for (j in 1:L) {
      S <- S + A1[j] / L / hN * K((x - Z[j]) / hN)
    }
    return(S)
  }
  Phat2 <- function(x, hN) {
    S <- 0
    for (j in 1:L) {
      S <- S + A2[j] / L / hN * K((x - Z[j]) / hN)
    }
    return(S)
  }
  
  ################ ISE ########################################################
  fu1 <- function(x, hN) {
    (dens1(x) - Phat1(x, hN))^2
  }
  
  fu2 <- function(x, hN) {
    (dens2(x) - Phat2(x, hN))^2
  }
  
  ISE1Min <- tryCatch(
    integrate(function(x) fu1(x, Min1), -40, 40)$value,
    error = function(e) NA_real_
  )
  
  ISE2Min <- tryCatch(
    integrate(function(x) fu2(x, Min2), -50, 50)$value,
    error = function(e) NA_real_
  )
  
  ################ Silverman ##################################################
  A <- as.matrix(cbind(A1, A2) / L)
  
  S <- silvrule(Z, A)
  
  Sil1 <- S[1]
  Sil2 <- S[2]
  
  ISE1Sil <- tryCatch(
    integrate(function(x) fu1(x, Sil1), -40, 40)$value,
    error = function(e) NA_real_
  )
  
  ISE2Sil <- tryCatch(
    integrate(function(x) fu2(x, Sil2), -50, 50)$value,
    error = function(e) NA_real_
  )
  
  if (any(is.na(c(ISE1Min, ISE2Min, ISE1Sil, ISE2Sil)))) {
    return(NULL)
  }
  
  #h_seq <- seq(0.1, 5, length.out = 100) 
  #cv1_values <- sapply(h_seq, CV1)
  #cv2_values <- sapply(h_seq, CV2)
  
  
  ################ Return #####################################################
  list(
    Min1 = Min1,
    Min2 = Min2,
    Sil1 = Sil1,
    Sil2 = Sil2,
    ISE1Min = ISE1Min,
    ISE2Min = ISE2Min,
    ISE1Sil = ISE1Sil,
    ISE2Sil = ISE2Sil
    #h_grid = h_seq,
    #cv1_plot = cv1_values,
    #cv2_plot = cv2_values
  )
}
