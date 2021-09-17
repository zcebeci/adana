# Kod 1.1: Tek tepeli fonksiyon grafiğinde tepeler ve vadileri bulma
findoptima = function(x, type="max", pflag=TRUE){
  if(type=="max"){
    if(pflag)
       results = which(diff(c(TRUE, diff(x)>=0,FALSE))<0)
    else
      results = which(diff(diff(x)>=0)<0)+1
  }else{
    if(pflag)
      results = which(diff(c(FALSE, diff(x)> 0, TRUE))>0)
    else
      results = which(diff(diff(x)>0)>0)+1
  }
  return(results)
}

# Kod 2.1: Onlu tamsayıdan ikili sayıya dönüştürme
int2bin = function(int, m){
  int = round(int)
  mc = floor(log(int, base=2)+1)
  if(missing(m)) m = mc
  bin = rep(0, m)
  i = 0
  while(int >= 1){
    i = i + 1
    bin[i] = int %% 2
    int = int %/% 2
  }
  return(rev(bin))
}

# Kod 2.2: İkili sayıdan onlu tamsayıya dönüştürme 
bin2int = function(bin) 
  sum(2^(which(rev(unlist(
  strsplit(as.character(bin), "")) == 1))-1))

# Kod 2.3: İkili sayıdan gri kodlu sayıya dönüştürme
bin2gray = function(bin){
  gray = rep(NA, length(bin))
  gray[1] = bin[1]
  for(i in 2:length(bin))
    gray[i]= xor(bin[i], bin[i-1])
  return(gray)
}

# Kod 2.4: Gri kodlu sayıdan ikili sayıya dönüştürme 1
gray2bin = function(gray){
  bin = rep(NA, length(bin))
  bin[1] = gray[1] 
  for(i in 2:length(bin))
    bin[i]= xor(gray[i], bin[i-1])
  return(bin)
}

# Kod 2.5: Gri kodlu sayıdan ikili sayıya dönüştürme 2
gray2bin2 = function(gray){
  bin = rep(NA, length(bin))
  bin[1] = bitvalue = gray[1] 
  for(i in 2:length(bin)){
    if(gray[i]==1) bitvalue=ifelse(bitvalue==1, 0, 1)
    else bitvalue=gray[i-1]
    bin[i]= bitvalue
  }
  return(bin)
}

# Kod 2.6: Gerçel sayıdan ikili sayıya dönüştürme
# Bağımlılık: Kod 2.1
encode = function(real, lb, ub, m){
  num = (real-lb)*(2^m-1)/(ub-lb)
  bin = int2bin(num, m)
  return(bin)
}

# Kod 2.7: İkili sayıdan gerçel sayıya dönüştürme
decode = function(bin, lb, ub, m){
  if(missing(lb) | missing(ub)) stop("lb, ub eksik")
  if(missing(m)) m=length(bin)
  real = lb + (ub-lb)/(2^m-1) * sum(c(bin)*c(2^seq(m-1,0)))
  return(real)
}

# Kod 2.8: Tamsayı vektörleri ikili vektörlere dönüştürme
# Bağımlılık: Kod 2.1
encode4int = function(x, M, ...){
  nx = length(x)
  xbin = c()
  for(i in 1:nx)
    xbin = c(xbin, int2bin(x[i], m=M[i]))
  return(xbin)
}

# Kod 2.9: Değişken uzunluklarını saptama
calcM = function(ub, ...){
  nx = length(lb)
  M = c()
  for(i in 1:nx)
    M = c(M, length(int2bin(ub[i])))
  return(M)
}

# Kod 2.10: İkili vektörleri tamsayı vektörlere dönüştürme
# Bağımlılık: Kod 2.2
decode4int = function(x, M, ...){
  nx = length(M)
  xint = integer(nx)
  j = 1 
  for(i in 1:nx){
    k = j+M[i]-1
    xint[i] = bin2int(x[j:k])
    j = k+1
  }
  return(xint)
}

# Kod 2.11: Gerçel sayı matrisini ikili kodlama
# Bağımlılık: Kod 2.6
encodepop = function(x, lb, ub, eps, ...){
  n = nrow(x)
  if(missing(lb)) stop("lb, alt sınır verilmeli")
  if(missing(ub)) stop("ub, üst sınır verilmeli")
  if(missing(eps)) eps=0.1
  v = length(lb)
  if(length(eps)==1) eps = rep(eps, v)
  m = floor(log((ub-lb)/eps,2)+1)
  M = sum(m)
  binmat = matrix(NA, nrow=n, ncol=M)
  for(i in 1:n){
    binvec = c()
    for(j in 1:v){
      vbin = encode(x[i,j], lb[j], ub[j], m[j])
      binvec = c(binvec, vbin)
    }
    binmat[i,] = binvec
  }
  return(list(binmat=binmat, m=m))
}

# Kod 2.12: İkili sayı matrisini gerçel sayıya dönüştürme
# Bağımlılık: Kod 2.7
decodepop = function(x, lb, ub, m,  ...){
  n = nrow(x)
  if(missing(lb)) stop("lb, alt sınır verilmeli")
  if(missing(ub)) stop("ub, üst sınır verilmeli")
  if(missing(m)) stop("krom. uz. (m) eksik")
  v = length(lb)
  xreal = matrix(NA, nrow=n, ncol=v)
  for(i in 1:n){
    realvec = c()
    for(j in 1:v){
      if(j==1){
         idx1 = 1
         idx2 = m[1]
      }else{
        idx1 = sum(m[1: (j-1)])+1
        idx2 = idx1+m[j]-1
      }
      vreal = decode(x[i,idx1:idx2], lb[j], ub[j], m[j])
      realvec = c(realvec, vreal)
    }
    xreal[i,] = realvec
  }
  return(xreal)
}

# Kod 2.13: İkili kodlama ile başlatma
initbin = function(n, m, prevpop, type, ...){
  if(missing(m)) m = 8
  if(missing(prevpop)) prevpop = NULL
  if(is.null(prevpop)){
     nprev = 0
  }else{
     prevpop = as.matrix(unname(prevpop))
     nprev = nrow(prevpop)
     m = ncol(prevpop)
  }
  if(missing(n))
    if(is.null(prevpop)) n = 4*m else n=nprev
  if(missing(type)) type = 1
  initpop = matrix(NA, ncol=m+2, nrow=n) # Başlangıç popülasyonu
  if(nprev==0){ #Rastlantısal başlatma
    for(i in 1:n)  
      initpop[i, 1:m] = sample(0:1, m, replace=TRUE) 
  }else if(nprev<n){ #Hibrit başlatma
    initpop[1:nprev,1:m] = prevpop[1:nprev,]
    for(i in (nprev+1): n) 
      initpop[i, 1:m] = sample(0:1, m, replace=TRUE) 
  }else if(nprev>n){ 
    initpop[1:n,1:m] = prevpop[1:n,]
  }else{ #Sezgisel başlatma
    initpop[1:nprev,1:m] = prevpop[1:nprev,]
  }
  initpop[,m+1] = 0
  if(nprev==0){
    rnames = paste0("T0", 1:n)
  }else if(nprev<n){
    rnames = c(paste0("Pre", 1:nprev), paste0("T0.", 1: (n-nprev)))
  }else{
    rnames = paste0("Pre", 1:n)
  }
  rownames(initpop) = rnames
  colnames(initpop) = c(paste0("gen", 1:m), "t", "fitval")
  if(type==2)
    initpop = initpop[,1:m]
  return(population=initpop)
}

# Kod 2.14: Değer kodlamalı başlatma 
initval = function(n, m, prevpop, lb, ub, nmode="real", type, ...){
  if(missing(m)){
    if(!missing(lb))
      m = length(lb)
    else
      m = 8 # Varsayılan gen sayısı
  }
  if(missing(prevpop)) prevpop = NULL
  if(is.null(prevpop)){
     nprev = 0
  }else{
     prevpop = as.matrix(unname(prevpop))
     nprev = nrow(prevpop)
     m = ncol(prevpop)
  }
  if(missing(n))
    if(is.null(prevpop)) n = 4*m else n=nprev
  if(missing(lb)) lb = rep(0, m)
  if(missing(ub)) ub = rep(1, m)
  if(length(lb) != length(ub)) 
     stop("ub, lb uzunluğu eşit olmalı!")
  if(missing(type)) type = 1
  initpop = matrix(NA, ncol=m+2, nrow=n) # Başlangıç popülasyonu
  if(nprev==0){
    for(i in 1:n) 
      for(j in 1:m)
        initpop[i,j] = runif(1, lb[j], ub[j]) 
  }else if(nprev<n){
    initpop[1:nprev,1:m] = prevpop[1:nprev,]
    for(i in (nprev+1): n) 
      for(j in 1:m)
        initpop[i,j] = runif(1, lb[j], ub[j]) 
  }else if(nprev>n){
    initpop[1:n,1:m] = prevpop[1:n,]
  }else{
    initpop[1:nprev,1:m] = prevpop[1:nprev,]
  }
  if(nmode=="integer") 
    initpop[,1:m] = apply(initpop[,1:m], 2, round)
  initpop[,m+1] = 0
  if(nprev==0){
    rnames = paste0("T0.", 1:n)
  }else if(nprev<n){
    rnames = c(paste0("Pre.", 1:nprev), 
      paste0("T0.", 1:(n-nprev)))
  }else{
    rnames = paste0("Pre.", 1:n)
  }
  rownames(initpop) = rnames
  colnames(initpop) = c(paste0("gen", 1:m), "t", "fitval")
  if(type==2)
    initpop = initpop[,1:m]
  return(population=initpop)
}

# Kod 2.15: Permütasyon kodlamalı başlatma
initperm = function(n, permset, prevpop, type, ...){
  if(missing(permset)) stperm=sample(0:9, 10, replace=FALSE)
  m = length(permset) # Gen sayısı
  if(missing(prevpop)) prevpop = NULL
  if(is.null(prevpop)){
     nprev = 0
  }else{
     prevpop = as.matrix(unname(prevpop))
     nprev = nrow(prevpop)
     m = ncol(prevpop)
  }
  if(missing(n))
    if(is.null(prevpop)) n = 4*m else n=nprev
  if(missing(type)) type = 1
  initpop = matrix(NA, ncol=m+2, nrow=n) #Başlangıç pop.
  if(nprev==0){
    for(i in 1:n) 
      initpop[i,1:m] = sample(permset, m, replace=FALSE)
  }else if(nprev<n){
    initpop[1:nprev,1:m] = prevpop[1:nprev,]
    for(i in (nprev+1): n) 
      initpop[i,1:m] = sample(permset, m, replace=FALSE)
  }else if(nprev>n){
    initpop[1:n,1:m] = prevpop[1:n,]
  }else{
    initpop[1:nprev,1:m] = prevpop[1:nprev,]
  }
  initpop[,m+1] = 0
  if(nprev==0){
    rnames = paste0("T0.", 1:n)
  }else if(nprev<n){
    rnames = c(paste0("Pre.", 1:nprev), 
      paste0("T0", 1:(n-nprev)))
  }else{
    rnames = paste0("Pre.", 1:n)
  }
  rownames(initpop) = rnames
  colnames(initpop) = c(paste0("gen", 1:m), "t", "fitval")
  if(type==2)
    initpop = initpop[,1:m]
  return(population=initpop)
}

# Kod 2.16: Başlatma fonksiyonu
# Bağımlık: Kod 2.13, 2.14, 2.15
initialize = function(initfunc, n, m, type, ...){
  func = as.character(match.call()[2])
  if(missing(type)) type=1
  dotargs = list(...)
  dotargs$n = n
  dotargs$m = m
  dotargs$type = type
  initpop=do.call(func, dotargs)
  return(initpop)
}

# Kod 2.17: Normal dağılışla başlatma
initnorm = function(n, m, pmean, psd, type, ...){
  if(missing(m)) m = 8 # Gen sayısı
  if(missing(n)) n = 4*m
  if(missing(pmean)) pmean = rep(0, m)
  if(missing(psd)) psd = rep(1, m)
  if(missing(type)) type = 1
  initpop = matrix(NA, ncol=m+2, nrow=n)
  for(i in 1:n)
    initpop[i,1:m] = rnorm(m, pmean, psd)
  initpop[,m+1] = 0
  rownames(initpop) = paste0("T0", 1:n)
  colnames(initpop) = c(paste0("gen", 1:m), "t", "fitval")
  if(type==2)
    initpop = initpop[,1:m]
  return(population=initpop)
}

# Kod 2.18a: MAXONE uyum fonksiyonu 1
maxone1 = function(x){
  sum1=0
  for(bits in x){
    if(bits) sum1=sum1+1  
  }
  return(sum1)
}

# Kod 2.18b: MAXONE uyum fonksiyonu 2
maxone = function(x) sum(x==1)

# Kod 2.18c: MAXONE uyum fonksiyonu 3
maxone2 = function(x) apply(x, 1, sum)

# Kod 2.19: MINONE uyum fonksiyonu 
# Bağımlılık: Kod 2.18b
minone = function(x) (-1*maxone(x))

# Kod 2.20: Kısıtlı uyum fonksiyonu örneği
penfit = function(x){ 
# Amaç fonksiyon
 fx = function(x) (60*x[1] + 30*x[2])
# Kısıt fonksiyonları
  g1 = function(x) (2*x[1] + x[2])
  g2 = function(x) (10*x[1] + 4*x[2])
  g3 = function(x) (x[1] + x[2])
# Uyum değerinin hesaplanması
  fitval = fx(x)
# Cezaların hesaplanması
  pencf = 1e+50  #Ceza katsayısı
  pen1 = ifelse(g1(x) <= 900, 0, g1(x)*pencf) #1. kısıt cezası
  pen2 = ifelse(g2(x) >= 1200, 0, g2(x)*pencf) #2. kısıt cezası
  pen3 = ifelse(g3(x) >= 300, 0, g3(x)*pencf) #3. kısıt cezası
# Uyumun değerinin cezalandırılması
  fitval = fitval-(pen1+pen2+pen3)
  return(fitval)
}

# Kod 2.21: Popülasyon uyum değerleri hesaplama
evaluate = function(fitfunc, population, objective, ...){
  if(missing(fitfunc)) stop("Uyum fonksiyonu verilmedi")
  if(missing(population)) stop("Popülasyon verilmedi")
  if(missing(objective)) objective="max"
  if(!is.matrix(population)) 
    population = matrix(population, 1, length(population))
  fitvals = rep(NA, nrow(population))
  for(i in 1:nrow(population))
    fitvals[i] = fitfunc(population[i,])
  if(objective=="min") 
    fitvals = -1*fitvals
  return(fitvals)
}

# Kod 2.22: Rastlantısal seçim
selrand = function(fitvals, ns, ...){
  n=length(fitvals)
  if(missing(ns)) ns=n
  matpool = sample(1:n, size=ns, replace=TRUE)
  return(matpool)
}

# Kod 2.23: Uyumla Oranlı Basit İadeli Rastlantısal Seçme (RSWRP)
elrswrp = function(fitvals, ns, ...){
  if (missing(ns)) ns = length(fitvals) #Popülasyon büyüklüğü
  n = ns #Çiftleşme havuzu büyüklüğü
  fitprobs = fitvals/sum(fitvals)  #Uyum olasılıkları
  matpool = sample(1:n, size=ns, prob=fitprobs, replace=TRUE)
  return(matpool)
}

# Kod 2.24: Rulet çarkı seçimi 1
selrws = function(fitvals, ns, ...){
  n = length(fitvals) #Popülasyon büyüklüğü
  if(missing(ns)) ns=n
  p = fitvals/sum(fitvals) 
  q = cumsum(p) #Kümülatif olasılıklar
  matpool = rep(NA, ns) #Çiftleşme havuzu
  i = 1  # Seçilen birey indeksi
  while(i <= ns){
    r = runif(1, 0, 1)
    j = 1
    while(q[j] < r){
      j = j+1
    }
    matpool[i] = j
    i = i+1
   }
   return(matpool)
}

# Kod 2.25: Rulet çarkı seçimi 2
selrws2 = function(fitvals, ns, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  p = fitvals/sum(fitvals)
  q = cumsum(p)  
  matpool = rep(NA, ns) 
  for(i in 1:ns){ 
    r = runif(1, 0, 1)
    matpool[i] = which(q >= r)[1]
  }
  return(matpool)
}

# Kod 2.26: Stokastik kesir kalanıyla seçim
# Bağımlılık: Kod 2.24
selrss = function(fitvals, ns, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  p = fitvals/sum(fitvals)
  expval = p*n
  intpart = floor(expval)
  rempart = expval-intpart
  matpool = c()
  for(i in 1:n)
    matpool = c(matpool, rep(i, intpart[i]))
  nfitvals = expval-intpart * sum(expval)/n
  ndif = n-length(matpool)
  remidx = selrws(nfitvals)[1:ndif]
  matpool = c(matpool, remidx)
  tmatpool = c()
  if(ns <= n){
    matpool = matpool[1:ns]
  }else{
    i=1
    repeat{
      for(j in matpool){
        tmatpool = c(tmatpool,j)
        i=i+1
        if(i>(ns-n)) break
      }
      if(i>(ns-n)) break
    }
    matpool = c(matpool, tmatpool)
  } 
  return(matpool)
}

# Kod 2.27: Stokastik evrensel seçim
selsus = function(fitvals, ns, ...){
  n = length(fitvals)
  if(missing(ns)) ns = n 
  if(ns<1) return(0)
  matpool = c()
  avgfit = sum(fitvals)/ns
  firstpointer = runif(1, 0, avgfit)
  i = 0:(ns-1)
  pointers = firstpointer + i * avgfit 
  for(pointer in pointers){
    idx = 0
    while(sum(fitvals[0:idx]) < pointer)
      idx=idx+1
    matpool = c(matpool, idx)
  }
  matpool = matpool[order(runif(ns))]
  return(matpool)
}

# Kod 2.28 Deterministik seçim
seldet = function(fitvals, ns, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  matpool = c()
  p = fitvals/sum(fitvals) #Uyum oranları
  expval = n*p #Beklenen değerler
  intpart = floor(expval) #Tam kısımlar
  rempart = expval-intpart #Kesirli kısımlar
  for(i in 1:n)
    matpool = c(matpool, rep(i, intpart[i])) # Tam kısımla seçim
  emppos = n-length(matpool) #Kesirli kısımla seçim
  if(emppos!=0){
    for(i in 1:emppos){
      remidx = which.max(rempart)
      matpool = c(matpool, remidx)
      rempart = rempart[-remidx]
    }
  }
  if(ns <= n){
    matpool = matpool[1:ns]
  }else{
    tmatpool = c()
    i=1
    while(i <= (ns-n)){
      for(j in matpool){
        tmatpool[i] = j
        i=i+1
        if(length(c(matpool, tmatpool))==ns) break
      }
    }
    matpool = c(matpool, tmatpool)
  } 
  return(matpool)
}

# Kod 2.29: Pencere ölçeklendirmesi 
selwscale = function(fitvals, ns, fmin, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  fstar = rep(NA, n)
  for(i in 1:n)
    fstar[i] = fitvals[i]-fmin
  p = fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.30: Sigma ölçeklendirmesi
selsscale = function(fitvals, ns, selc, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(selc)) selc=2
  favg = mean(fitvals)
  fsd = sd(fitvals)
  fstar = rep(NA, n)
  for(i in 1:n)
    fstar[i] = ifelse(fsd!=0, 1+(fitvals[i]-favg)/selc*fsd, 1.0)
  fstar[fstar<0]=0
  p = fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.31: Sigma ölçeklendirmesi 2
selsscale2 = function(fitvals, ns, selc, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(selc)) selc=2
  if(any(fitvals<0)) fitvals=fitvals-min(fitvals)
  favg = mean(fitvals)
  fsd = sd(fitvals)
  fstar = fitvals + (favg-selc*fsd)
  fstar[fstar<0]=0
  p = fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.32: Doğrusal uyum ölçeklendirmesi
sellscale = function(fitvals, ns, sells, ...){
  n=length(fitvals)
  if(missing(ns)) ns=n
  if(missing(sells)) sells=1.5 #Ölçeklendirme faktörü
  fmin = min(fitvals)
  fmax = max(fitvals)
  favg = mean(fitvals)
  fstar = rep(NA, n)
  if(sells > (1+(fmax-favg)/(favg-fmin)))
    ms = (fmax-favg)/(favg-fmin)
  else
    ms = sells-1.0
  fstar = 1+ms*(fitvals-favg)/(fmax-favg) #Ölçeklendirilmiş uyum
  p=fstar/sum(fstar)
  matpool=sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.33: Sıra sayısıyla ölçeklendirme
selrscale = function(fitvals, ns, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  fitranks = rank(fitvals, ties.method="min")
  s = max(fitvals)/median(fitvals)
  fstar = s-2*(fitranks-1)*(s-1)/(n-1)
  p = fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.34: Sıra sayısıyla ölçeklendirme 2
selrscale2 = function(fitvals, ns, sels, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(sels)) sels=1.5
  fitranks = rank(fitvals, ties.method="min")
  fstar = (2-sels)/n + (2*fitranks*(sels-1))/(n*(n-1))
  p = fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.35: Güç ölçeklendirmesi
selpscale = function(fitvals, ns, selk, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(selk)) selk=1.005
  fstar = fitvals^selk
  p = fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.36: Üs ölçeklendirmesi
selescale = function(fitvals, ns, selb, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(selb)) selb=0.5
  fstar = exp((selb*fitvals))
  p = fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.37: Turnuva seçimi
seltour = function(fitvals, ns, selt, reptype, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(ns<2) ns=2
  if(missing(selt)) selt=2
  if(selt<2 | selt>n) selt=2
  if(missing(reptype)) reptype=FALSE
  matpool = rep(NA, ns)
  for(i in 1:ns){
    tourgroup = sample(1:n, size=selt, replace=reptype)
    matpool[i] = tourgroup[which.max(fitvals[tourgroup])]
  }
  return(matpool)
}

# Kod 2.38: Turnuva seçimi 2
seltour2 = function(fitvals, ns, selt, reptype, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(ns<2) ns=2
  if(missing(selt)) selt=2
  if(selt<2 | selt>n) selt=2
  if(missing(reptype)) reptype=FALSE
  matpool = c()
  tmppool = 1:n
  for(i in 1:ns){
    if(length(tmppool) < selt) tmppool = 1:n
    tourgroup = sample(tmppool, size=selt, replace=reptype)
    bestidx = which(fitvals==max(fitvals[tourgroup])) 
    matpool[i] = sample(bestidx, size=1)
    tmppool = tmppool[-matpool]
  }
  return(matpool)
}

# Kod 2.39: Boltzmann turnuva seçimi
selboltour = function(fitvals, ns, selt0, selg, selgmax, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(selt0)) selt0=50
  if(selt0<5 | selt0>100) return(0)
  if(missing(selg) | missing(selgmax)) return(0) 
  fstar=rep(0,n)
  fmax=max(fitvals)
  for(i in 1:n){
      alfa=runif(1)
      k=1+selg/selgmax*100
      temp = selt0 * (1-alfa)^k
      fstar[i] = exp(-(fmax-fitvals[i])/temp)
  }
  p=fstar/sum(fstar)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.40: Doğrusal sıra seçimi 1
sellrs = function(fitvals, ns, sels, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(sels)) sels=1.5 #Seçim basıncı orta
  if(sels<1 | sels>2) return(0)
  fitranks = rank(fitvals, ties.method="min")
  fitstars = (2-sels)+2*(sels-1)*((fitranks-1)/(n-1))
  p = fitstars/sum(fitstars)
  q = cumsum(p) #Kümülatif uyum oranları
  matpool = rep(NA, ns) #Çiftleşme havuzu
  i = 1  # Seçilen birey indeksi
  while(i <= ns){
    r = runif(1, 0, 1)
    j = 1
    while(q[j] < r){
      j = j+1
    }
    matpool[i] = j
    i = i+1
   }
   return(matpool)
}

# Kod 2.41: Doğrusal sıra seçimi 2
sellrs2 = function(fitvals, ns, ...){
  n=length(fitvals)
  if(missing(ns)) ns=n
  fitranks = rank(fitvals, ties.method="min")
  p = 2*fitranks/(n*(n+1)) 
  matpool=sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.42: Doğrusal sıra seçimi 3
sellrs3 = function(fitvals, ns, sels, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(sels)) sels=1.5 #Seçim basıncı orta
  if(sels<1 | sels>1.99) sels=1.5
  fitranks = rank(fitvals, ties.method="min")
  p = 1/n*(sels-2*(sels-1)*((fitranks-1)/(n-1)))
  matpool = sample(1:n, size=ns, prob=p)
  return(matpool)
}

# Kod 2.43: Doğrusal olmayan sıra seçimi
selnlrs = function(fitvals, ns, selns, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(selns)) selns=0.5
  fitranks = (n+1)-rank(fitvals, ties.method="min")
  fstar= selns*(1-selns)^(fitranks-1)
  p = pmin(pmax(0, fstar/sum(fstar)), 1)
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.44: Üssel sıra seçimi
selers = function(fitvals, ns, selbc, ...){
  n = length(fitvals)
  if(missing(selbc)) selbc=0.5
  if(missing(ns)) ns=n
  fitranks = rank(fitvals, ties.method="min")
  p = selbc^(n-fitranks) / sum(selbc^(n-fitranks))
  matpool = sample(1:n, size=ns, prob=p, replace=TRUE)
  return(matpool)
}

# Kod 2.45: Kesme seçimi
seltrunc = function(fitvals, ns, selps, ...){
  n = length(fitvals)
  if(missing(ns)) ns=n
  if(missing(selps)) selps=0.5
  if(selps<=0 | selps>1) return(0) 
  ts = ifelse(round(n*selps)<2, 2, round(n*selps))
  matpool = c()
  names(fitvals) = 1:n
  sfitvals = sort(fitvals, decreasing=TRUE)
  if(ts>=ns){
    matpool[1:ns] = as.integer(names(sfitvals)[1:ns])
  }else{
    matpool[1:ts] = as.integer(names(sfitvals)[1:ts])
    for(i in 1: (ns-ts)){
      for(j in matpool){
        if(length(matpool)==ns) break
        matpool = c(matpool, j)
      }
    }
  }
  return(matpool)
}	

# Kod 2.46: Çiftleşme havuzuna ebeveyn seçme
select = function(selfunc, fitvals, ns, selb, selbc,
    selc, selk, sells, selns, selps, sels, selt, 
    selt0,selw, selg, selgmax, fmin, reptype, ...){
  args = as.list(match.call())[-1]
  nargs = length(args)
  if(missing(selfunc)) selfunc=seltour
  if(missing(reptype)) reptype=FALSE
  if(missing(ns)) ns=length(fitvals)
  if(missing(selw)) selw=2
  if(missing(selt)) selt=2
  if(missing(selb)) selb=0.5
  if(missing(selb)) selbc=0.5
  if(missing(selk)) selk=1.005
  if(missing(selc)) selc=0.5
  if(missing(sels)) sels=1.5
  if(missing(sells)) sells=1.5
  if(missing(selns)) selns=0.5
  if(missing(selps)) selps=0.5
  if(missing(selt0)) selt0=50
  selected = do.call(as.character(args[1]), args[2:nargs])
  return(selected)
}

# Kod 2.47: Tek noktalı çaprazlama
px1 = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    v = sample(0:(m-1), size=1) #Kesme noktası
    if(v == 0){
      y1 = x2
      y2 = x1
    }else if(v>0 & v<m){
      y1 = c(x1[1:v], x2[(v+1):m])
      y2 = c(x2[1:v], x1[(v+1):m])
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.48: Çok noktalı çaprazlama
kpx = function(x1, x2, cxon, cxk, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxk)) cxk = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    v = sort(sample(1:(m-1), size=cxk, replace=FALSE)) 
    y1 = y2 = rep(NA, m)
    for(j in 1:v[1]){
      y1[j] = x1[j]
      y2[j] = x2[j]
    }
    direction = 1 
    for(j in 2:cxk){
      for(k in (v[j-1]+1):v[j]){
        if(direction){
          y1[k] = x2[k]
          y2[k] = x1[k]
        }else{
          y1[k] = x1[k]
          y2[k] = x2[k]
        }
      } 
      direction = ifelse(direction,0,1)
    }
    for(j in (v[cxk]+1):m){
      if(direction){
        y1[j] = x2[j]
        y2[j] = x1[j]
      }else{
        y1[j] = x1[j]
        y2[j] = x2[j]
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.49: Karıştırmalı çaprazlama 
sc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  sidx = sample(1:m, size=m) # Karıştırma işlemi
  for(i in seq(from=1, to=cxon, by=2)){
    v = sample(1: (m-1), 1)   # Kesme noktası
    y1 = y2 = rep(NA, m) # Yavrular
    for(j in 1:m){
      if(j < v) {
        y1[sidx[j]] = x1[sidx[j]]
        y2[sidx[j]] = x2[sidx[j]]
      }else{
        y1[sidx[j]] = x2[sidx[j]]
        y2[sidx[j]] = x1[sidx[j]]
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.50 İndirgenmiş ebeveynlik çaprazlaması 
rsc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    pv = c() # Olası kesme noktaları vektörü
    j = 1
    for(k in 1:m){
      if(x1[k] != x2[k]){
         pv[j] = k
         j = j + 1
      }
    }
    y1 = x1; y2 = x2 # Yavrular
    if(length(pv)>0){
      v = sample(pv, size=1) #Kesme noktası
      for(k in 1:m){
        if(k <= v){
          y1[k] = x1[k]
          y2[k] = x2[k]
        }else{
          y1[k] = x2[k]
          y2[k] = x1[k]
        }
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.51: Sezgisel tekdüze çaprazlama
hux = function(x1, x2, cxon, cxps, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxps)) cxps = 0.5
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = x1
    y2 = x2
    ndg = 0  # Farklı genlerin sayısı
    for(j in 1:m){
      if(y1[j] != y2[j]){
        ndg = ndg +1
      }
    }
    sc = 0  # Değiştirme sayacı
    while(sc <= ndg/2){
      for(j in 1:m){
        if(y1[j] != y2[j] & y1[j] != x2[j]){
          v = runif(1, 0, 1)
          if(v > cxps){
            y1[j] = x2[j]
            y2[j] = x1[j]
            sc = sc + 1
          }
        }
      }   
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.52: Tekdüze (Uniform) çaprazlama 1
ux = function(x1, x2, cxon, cxps, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxps)) cxps = 0.5
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = x1
    y2 = x2
    p = runif(m, 0, 1) # Olasılıklar vektörü
    for(j in 1:m){
      if(p[j] > cxps){
        y1[j] = x2[j]
        y2[j] = x1[j]
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.53: Tekdüze (Uniform) çaprazlama 2
ux2 = function(x1, x2, cxon, cxps, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxps)) cxps = 0.5
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = x1
    y2 = x2
    p = runif(m, 0, 1) # Olasılıklar vektörü
    y1[p > cxps] = x2[p > cxps]
    y2[p > cxps] = x1[p > cxps]
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.54: Maske çaprazlaması
mx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = x1
    y2 = x2
    p1 = round(runif(m, 0, 1))
    p2 = round(runif(m, 0, 1))
    for(j in 1:m){
      if(!p1[j] & p2[j]){
        y1[j] = x2[j]
      }else if(p1[j] & !p2[j]){ 
        y2[j] = x1[j]
      }
    }  
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.55: Saygılı Çaprazlama (RRC)
rrc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=1)){
    y = x1 & x2
    v = which(x1!=x2)
    for(j in v)
      y[j] = ifelse(runif(1)>0.5, 1, 0)
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.56: Saygısız Çaprazlama (DISC)
disc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=1)){
    y = rep(NA, m)
    v = sample(1:m, 1)
    for(j in 1:v){
      if(x1[j] != x2[j])
        y[j] = x1[j]
      else
        y[j] = sample(0:1, 1)
    }
    for(j in (v+1):m){
      if(x1[j] != x2[j])
        y[j] = x2[j]
      else
        y[j] = sample(0:1, 1)
    }
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.57: Asimetrik İki Noktalı Çaprazlama (ATC)
atc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    v = c()
    y1 = y2 = rep(NA, m)
    v = sort(sample(1:m, 2, replace=FALSE))
    v[3] = sample(1:m, 1)
    j = v[3]
    for(k in 1:m){
      if(k<v[1])
        y1[k] = x1[j]
      else if(v[1]<=k & k<=v[2]){
        y1[k] = x2[j]
        j = j+1
        if(j>m) j=1
      }else{
        y1[k] = x1[k]
      }
    }
    for(k in 1:m){
      if(k<v[1])
        y2[k] = x2[k]
      else if(v[1]<=k & k<=v[2]){
        y2[k] = x1[k]
      }else{
        y2[k] = x2[k]
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.58: Sayı Korumalı Çaprazlama (CPC)
cpc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  lup = c(); ldown=c(); ll=0
  for(i in seq(from=1, to=cxon, by=2)){
    for(j in 1:m){
      if(x1[j]==1 & x2[j]==0){
       lup = c(lup, j)
       ll = ll+1
      }else if(x1[j]==0 & x2[j]==1){
       ldown = c(ldown, j)
      }
    }
    y1 = x1 ; y2 = x2
    for(j in ll){
      r = runif(1, 0, 1)
      if(r<0.5){
        temp = y1[lup[j]]
        y1[lup[j]] = y2[lup[j]]
        y2[lup[j]] = temp
        temp = y1[ldown[j]]
        y1[ldown[j]] = y2[ldown[j]]
        y2[ldown[j]] = temp
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.59: Bağlantı/Değiştokuş Çaprazlaması (EC,LC)
eclc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=1)){
    y = x2
    v = sort(sample(1:m, size=2, replace=FALSE))
    segx1 = x1[v[1]:v[2]]
    k = sample(1:m, 1)
    for(j in segx1){
      y[k] = j
      k = k+1
      if(k>m) k = 1
    }
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.60: Rastlantısal Ve/Veya Çaprazlaması (RAOC)
raoc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  y1 = y2 = rep(NA, m)
  alfa = runif(1, 0, 1)
  for(i in seq(from=1, to=cxon, by=2)){
    for(j in 1:m){
      r = runif(1, 0, 1)
      if(r>alfa){
        y1[j]= ifelse(x1[j] & x2[j],1,0)
        y2[j]= ifelse(x1[j] | x2[j],1,0)
      }else{
        y1[j]= ifelse(x1[j] | x2[j],1,0)
        y2[j]= ifelse(x1[j] & x2[j],1,0)
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.61: Kesikli çaprazlama
dc = function(x1, x2, cxon, cxps, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxps)) cxps = 0.5
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  y1 = x1; y2 = x2
  for(i in seq(from=1, to=cxon, by=2)){
    v = runif(1, 0, 1)
    for(j in 1:m){
      if(v >= cxps){
        y1[j] = x2[j]
        y2[j] = x1[j]
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.62: Ortalama çaprazlama
ax = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=1)){
    y = (x1 + x2)/2
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.63: Sezgisel aritmetik çaprazlama
hc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  x1fit = fitfunc(x1)
  x2fit = fitfunc(x2)
  for(i in seq(from=1, to=cxon, by=1)){
    r = runif(1, 0, 1)
    if(x2fit>=x1fit)
      y  = r * (x2-x1) + x2
    else
      y  = x1
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.64: Tekli aritmetik çaprazlama 
sax = function(x1, x2, cxon, cxalfa, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxalfa)) cxalfa = runif(1, 0, 1)
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = y2 = rep(NA, m)
    v = sample(1:m, 1)
    y1[1:v] = x1[1:v]
    y2[1:v] = x2[1:v]
    y1[(v+1):m] = (1-cxalfa) * x2[(v+1):m]
    y2[(v+1):m] = (1-cxalfa) * x1[(v+1):m]
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.65: Tam aritmetik çaprazlama 
wax = function(x1, x2, cxon, cxalfa, ...){
  n = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=n)
  for(i in seq(from=1, to=cxon, by=2)){
    if(missing(cxalfa)) cxalfa = runif(1, 0, 1)
    y1 = cxalfa * x1 + (1-cxalfa) * x2
    y2 = (1-cxalfa) * x1 + cxalfa * x2
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.66: Lokal aritmetik çaprazlama
lax = function(x1, x2, cxon, ...){
  n = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=n)
  for(i in seq(from=1, to=cxon, by=2)){
    x = matrix(c(x1, x2), nrow=2, ncol=n, byrow=TRUE)
    y = matrix(NA, nrow=2, ncol=n)
    cxalfa = runif(n) #Alfa her gen için farklı
    y[1,] = cxalfa*x[1,] + (1-cxalfa)*x[2,]
    y[2,] = cxalfa*x[2,] + (1-cxalfa)*x[1,]
    offsprings[i,] = y[1,]
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y[2,]
  }
  return(offsprings)
}

# Kod 2.67: Düz çaprazlama 
bx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon=2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in 1:cxon){
    for(j in 1:m){
      selval = runif(1, min(x1[j], x2[j]), max(x1[j], x2[j]))
      offsprings[i,j] = selval
    }
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.68: Geliştirilmiş kutu çaprazlaması 
ebx = function(x1, x2, lb, ub, cxon, cxalfa, ...){
  m = length(x1)
  if(missing(cxon)) cxon=2
  if(missing(cxalfa)) cxalfa = runif(1, 0, 1)
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in 1:cxon){
    for(j in 1:m){
      minx = min(x1[j],x2[j])
      maxx = max(x1[j],x2[j])
      emin = minx - cxalfa*(maxx-minx)
      emax = maxx + cxalfa*(maxx-minx)
      offsprings[i,j] = runif(1, min(emin,lb[j]), max(emax, ub[j]))
    }
  }
  return(offsprings)
}

# Kod 2.69: Harmanlanmış çaprazlama (BLX-α)
blxa = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  y1 = y2 = rep(NA, m)
  d = abs(x1-x2)
  for(i in seq(from=1, to=cxon, by=2)){
    for(j in 1:m){
      a = runif(1, 0, 1)
      y1[j] = runif(1, min(x1[j], 
              x2[j])-a*d[j], max(x1[j], x2[j])+a*d[j])
      y2[j] = runif(1, min(x1[j], 
              x2[j])-a*d[j], max(x1[j], x2[j])+a*d[j])
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.70: Harmanlanmış Çaprazlama (Alfa-Beta)
blxab = function(x1, x2, cxon, ...){
  n = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=n)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = y2 = rep(NA, n)
    d = abs(x1-x2)
    for(j in 1:n){
      if(x1[j] <= x2[j]){
        alfa = runif(1, 0, 1)
        beta = runif(1, 0, 1)
        y1[j] = runif(1, x1[j]-alfa*d[j], x2[j]+beta*d[j])
        alfa = runif(1, 0, 1)
        beta = runif(1, 0, 1)
        y2[j] = runif(1, x1[j]-alfa*d[j], x2[j]+beta*d[j])
      }else{
        alfa = runif(1, 0, 1)
        beta = runif(1, 0, 1)
        y1[j] = runif(1, x2[j]-beta*d[j], x2[j]+alfa*d[j])
        alfa = runif(1, 0, 1)
        beta = runif(1, 0, 1)
        y2[j] = runif(1, x2[j]-beta*d[j], x2[j]+alfa*d[j])
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.71: Laplace çaprazlaması
lapx = function(x1, x2, cxon, cxa, cxb, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxa)) cxa = 0.0
  if(missing(cxb)) cxb = 0.15
  if(length(cxa) == 1) cxa = rep(cxa, m)
  if(length(cxb) == 1) cxb = rep(cxb, m)
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    x = matrix(c(x1,x2), nrow=2, ncol=m, byrow=TRUE)
    y = matrix(NA, nrow=2, ncol=m)
    pr1 = runif(m)
    pr2 = runif(m)
    beta = cxa + ifelse(pr1 > 0.5, cxb*log(pr2), -cxb*log(pr2))
    bp = beta*abs(x[1,] - x[2,])
    y[1,] = pmin(pmax(x[1,] + bp, range(x[1,])[1]), range(x[1,])[2])
    y[2,] = pmin(pmax(x[2,] + bp, range(x[2,])[1]), range(x[2,])[2])
    offsprings[i,] = y[1,]
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y[2,]
  }
  return(offsprings)
}

# Kod 2.72: Genişletilmiş hat çaprazlaması (ELX)
elx = function(x1, x2, lb, ub, cxon, cxealfa, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxealfa)) cxealfa = 1
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  v1 = -Inf; v2=Inf
  for(i in seq(from=1, to=cxon, by=1)){
    y = rep(NA, m)
    for(j in 1:m){
       if(x1[j]!=x2[j]){
         tl = (lb[j]-x2[j])/(x1[j]-x2[j])
         tu = (ub[j]-x2[j])/(x1[j]-x2[j])
         tmin = min(tl, tu)
         tmax = max(tl, tu)
         v1 = max(v1,tmin)
         v2 = min(v2,tmax)
       }   
    }
    lambda = runif(1, max(v1,-cxealfa), min(v2, 1+cxealfa))
    y = lambda * x1 + (1-lambda)*x2
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.73: Geometrik çaprazlama
geomx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=1)){
    alfa=runif(1)
    y = x1^alfa * x2^(1-alfa)
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.74: Küre çaprazlaması
spherex = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=1)){
    alfa=runif(1)
    y = sqrt(alfa*x1^2 + (1-alfa)*x2^2)
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.75: Kısmi eşleşmeli çaprazlama
pmx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    x = matrix(c(x1,x2), nrow=2, ncol=m, byrow=TRUE)
    y = matrix(NA, nrow=2, ncol=m)
    v = sort(sample(1:m, size=2))
    y[1:2,v[1]:v[2]] = x[2:1,v[1]:v[2]]
    for(j in setdiff(1:m, v)){
      if(!any(x[2,j] == y[1,v])){
        y[1,j] = x[2,j]
      }
      if(!any(x[1,j] == y[2,v])){
        y[2,j] = x[1,j] 
      }
    }
    y[1, is.na(y[1,])] = setdiff(x[2,], y[1,])
    y[2, is.na(y[2,])] = setdiff(x[1,], y[2,])
    offsprings[i,] = y[1,]
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y[2,]
  }
  return(offsprings)
}

# Kod 2.76: Değiştirilmiş kısmi eşleşmeli çaprazlama
mpmx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    v = sort(sample(2: (m-1), size=2)) 
    y1 = y2 = rep(NA, m)
# Orta parça
    y1[v[1]:v[2]] = x1[v[1]:v[2]]
    y2[v[1]:v[2]] = x2[v[1]:v[2]]
# Sol parça
    for(j in 1: (v[1]-1)){
      if(!(x2[j] %in% y1)) y1[j] = x2[j]
      if(!(x1[j] %in% y2)) y2[j] = x1[j]
    }
# Sağ parça
    for(j in (v[2]+1):m){
      if(!(x2[j] %in% y1)) y1[j] = x2[j]
      if(!(x1[j] %in% y2)) y2[j] = x1[j]
    }
# Rastlantısal permütasyon
    y1[is.na(y1)] = sample(setdiff(x2, y1[!is.na(y1)])) 
    y2[is.na(y2)] = sample(setdiff(x1, y2[!is.na(y2)]))
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.77: Tekdüze kısmi eşleşmeli çaprazlama
upmx = function(x1, x2, cxon, ...){
  m = length(x1)
  k = ceiling(m/3)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = x1; y2=x2
    j = 1
    while(j<=k){
      v1 = sample(1:m, 1) 
      gen1 = x2[v1]; gen2 = x1[v1]
      v21 = which(y1==gen1)
      v22 = which(y2==gen2)
      tmp = y1[v1]
      y1[v1] = y1[v21]
      y1[v21] = tmp
      tmp = y2[v1]
      y2[v1] = y1[v22]
      y2[v22] = tmp
      j = j+1
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.78: Sıra çaprazlaması (OX)
ox = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    x = matrix(c(x1,x2), nrow=2, ncol=m, byrow=TRUE)
    y = x
    v = sort(sample(2:(m-1), size=2, replace=TRUE))
    for(j in 1:2){
      part2 = x[j, v[1]:v[2]]
      revx <- x[j,-c(v[1]:v[2])]
      revx <- rev(revx)
      part1 = revx[1: (v[1]-1)]
      part3 = revx[(length(part1)+1):length(revx)] 
      y[i,] = c(part1, part2, part3)
    }
    offsprings[i,] = y[1,]
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y[2,]
  }
  return(offsprings)
}

# Kod 2.79: Sıraya dayalı çaprazlama (OX2)
ox2 = function(x1, x2, cxon, cxoxk, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  if(missing(cxoxk)) cxoxk = max(1, sample(2:round(m/2),1))
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    v = sample(1:m, size=cxoxk, replace=FALSE)
    y1 = x2
    gv = x1[v]
    idx = c()
    for(k in gv)
      idx = c(idx, which(k==y1))
    y1[rev(idx)] = gv
    y2 = x1
    gv = x2[v]
    idx = c()
    for(k in gv)
      idx = c(idx, which(k==y2))
    y2[rev(idx)] = gv
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.80: Maksimal koruyucu çaprazlama
mpx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = y2 = rep(NA, m)
    j = m
    while(j > round(m/2)){
      v = sort(sample(1:m, size=2, replace=FALSE))
      j = length(v[1]:v[2])
    }
    k = j+1
    y1[1:j] = x1[v[1]:v[2]]
    y2[1:j] = x2[v[1]:v[2]]
    y1[k:m] = setdiff(x2, x1[v[1]:v[2]])
    y2[k:m] = setdiff(x1, x2[v[1]:v[2]])
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.81: Kenar rekombinasyonu
erx = function(x1, x2, cxon, ...){
# Kenarları bulma
findedges = function(x){
  m = length(x)
  edges = list()
  for(i in 1:m){
    edge = c()
    if(i == 1)
      edge = c(edge, x[i+1], x[m])
    else if(i < m)
      edge = c(edge, x[i-1], x[i+1])
    else
      edge = c(edge, x[i-1], x[1])
    edges[[x[i]]]=edge
  }
  return(edges)
}
# Kenarları birleştirme
mergeedges = function(x1, x2){
  m = length(x1)
  edges = list()
  for(i in 1:m)
    edges[[i]]=sort(unique(c(x1[[i]], x2[[i]])))
  return(edges)
}
# Çaprazlama 
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  edges1 = findedges(x1)
  edges2 = findedges(x2)
  edges = mergeedges(edges1, edges2)
  for(i in seq(from=1, to=cxon, by=1)){
    tempedges = edges
    y = rep(NA, m)  
    idx = sample(1:m,1)
    y[1] = idx
    k = 2
    while(k<=m){
      for(j in 1:m)
        tempedges[[j]] = setdiff(tempedges[[j]], idx)
      nidx = Inf
      neighbors = tempedges[[idx]]
      for(j in neighbors){
        if(length(tempedges[[j]]) < nidx){
          nidx = length(tempedges[[j]])
          midx = j
        }
      }
      idx = midx
      y[k] = idx
      k = k + 1
    }
    offsprings[i,] = y
    if(i==cxon & cxon%%2==1) break
  }
  return(offsprings)
}

# Kod 2.82: Konuma dayalı çaprazlama (PBX)
pbx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    x = matrix(c(x1,x2), nrow=2, ncol=m, byrow=TRUE)
    y = matrix(rep(NA, m), nrow=2, ncol=m)
    v = unique(sample(1:m, size=m, replace=TRUE))
    y[1,v] = x[2,v]
    y[2,v] = x[1,v]
    for(j in 1:2){
      cidx = which(is.na(y[j,]))
      sdif = setdiff(x[j,], y[j,v])
      y[j,cidx] = sdif
    }
    offsprings[i,] = y[1,]
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y[2,]
  }
  return(offsprings)
}

# Kod 2.83: Konuma dayalı çaprazlama 2(PBX)
pbx2 = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    x = matrix(c(x1,x2), nrow=2, ncol=m, byrow=TRUE)
    y = matrix(rep(NA, m), nrow=2, ncol=m)
    s = round(runif(1, 1, m))
    v = sample(1:m, size=s, replace=FALSE)
    y[1,v] = x[2,v]
    y[2,v] = x[1,v]
    for(j in 1:2){
      cidx = which(is.na(y[j,]))
      sdif = setdiff(x[j,], y[j,v])
      y[j,cidx] = sdif
    }
    offsprings[i,] = y[1,]
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y[2,]
  }
  return(offsprings)
}

# Kod 2.84: Döngü çaprazlaması (CX)
cx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = y2 = rep(NA, m)
    j = 1
    y1[j] = x1[j]
    y2[j] = x2[j]
    repeat{
      k = which(x1==x2[j])
      y1[k] = x2[j]
      y2[k] = x2[k]
      j = k
      if( (x2[j] %in% y1))
        break
    }
    for(j in 1:m){
      if(is.na(y1[j])){
        y1[j] = x2[j]
      }
      if(is.na(y2[j])){
        y2[j] = x1[j]
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.85: İyileştirilmiş Döngü Çaprazlaması (ICX)
icx = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = y2 = rep(NA, m)
    y1[1] = x2[1]
    y2[1] = x1[1]
    j1 = j2 = 1
    k1 = k2 = 2
    repeat{
      if(all(x1==x2)){  #Tüm genler eşitse
        y1=y2=x1
        break
      }
      idx1 = which(x1==y1[j1])
      idx2 = which(x2==y2[j2])
      y1[k1] = x2[idx1]
      y2[k2] = x1[idx2]
      if(x1[1]==y1[k1] & k1<m){   #Tamamlanmamış genler varsa
        remx1 = setdiff(x1, y1[!is.na(y1)])
        remx2 = setdiff(x2, y1[!is.na(y2)])
        y1[(k1+1):m] = remx2
        y2[(k2+1):m] = remx1
        break
      }
      j1 = k1
      j2 = k2
      k1 = k1+1
      k2 = k2+1
      if(k1>m) break
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.86: Sinüs hareketi çaprazlaması (SMC)
smc = function(x1, x2, cxon, ...){
  m = length(x1)
  if(missing(cxon)) cxon = 2
  offsprings = matrix(NA, nrow=cxon, ncol=m)
  for(i in seq(from=1, to=cxon, by=2)){
    y1 = y2 = rep(NA, m)
    y1[1] = x1[1]
    j = 1
    k1 = 2
    k2 = 1
    repeat{
      if(x2[j] %in% y1[!is.na(y1)]){
        y2[k2]=x2[j]
        k2 = k2+1
      }else{
        y1[k1]=x2[j]
        k1 = k1+1
      }
      j = j+1
      if(x1[j] %in% y1[!is.na(y1)]){
        y2[k2]=x1[j]
        k2 = k2+1
      }else{
        y1[k1]=x1[j]
        k1 = k1+1
      }
      if(j==m){
       y2[m] = x2[m]
       break
      }
    }
    offsprings[i,] = y1
    if(i==cxon & cxon%%2==1) break
    offsprings[i+1,] = y2
  }
  return(offsprings)
}

# Kod 2.87: Çaprazlama
cross = function(crossfunc, matpool, cxon, cxpc, gatype, ...){
  func = as.character(match.call()[2])
  if(missing(gatype)) gatype="gga"
  if(missing(cxpc)) cxpc=0.95  #Çaprazlama oranı
  if(missing(cxon)) cxon=2    #Çiftleşme başına yavru sayısı
  dargs = list(...)
  dargs$cxon = cxon
  nm = nrow(matpool) # Çiftleşme havuzundaki birey sayısı
  mm = ncol(matpool) # Gen sayısı
  nc = round(nm*cxpc) 
  nc = ifelse(nc%%2==1, nc-1, nc)
  nc = ifelse(nc>=nm, nm-1, nc) # Çiftleşecek birey sayısı
  nc = ifelse(nc<2, 2, nc) # Çiftleşecek birey sayısı
  if(gatype=="gga")  
    offsprings = matrix(NA, nrow=(nc/2*cxon), ncol=mm) # Yavrular
  else
    offsprings = matrix(NA, nrow=cxon, ncol=mm) # Yavrular
  mc = mm-2
  i = 1
  for(j in seq(from=1, to=nc-1, by=2)){
    dargs$x1 = unname(matpool[j,1:mc])
    dargs$x2 = unname(matpool[j+1,1:mc])
    crossresult = do.call(func, dargs)
    for(k in 1:cxon){
      offsprings[i,1:mc] = crossresult[k,]
      i=i+1	
    }
    if(gatype=="ssga"){  #SSGA ise çık
      return(offsprings)
    }
  }
  if(nc<nm)
    offsprings = rbind(offsprings, matpool[(nc+1):nm,])
  return(offsprings)
}

# Kod 2.88: Bit dönüştürme mutasyonu
bitmut = function(y, ...){
  m = length(y)
  v = sample(1:m, 1)
  y[v] = ifelse(y[v], 0, 1)
  return(list(mutant=y, mutgen=v))
}

# Kod 2.89: Rastlantısal değiştirme mutasyonu
randmut = function(y, lb, ub, ...){
  m = length(y)
  if(missing(lb) | missing(ub)) 
    stop("randmut için lb ve ub gerekli!")
  j = sample(1:m, 1)
  y[j] = runif(1, lb[j], ub[j])
  return(list(mutant=y, mutgen=j))
}

# Kod 2.90: Rastlantısal mutasyon 2
randmut2 = function(y, lb, ub, mutpm, ...){
  m = length(y)
  if(missing(lb) | missing(ub))
    stop("randmut2 için lb ve ub gerekli")
  if(missing(mutpm)) mutpm=0.05
  for(j in 1:m){
    r = runif(1, 0, 1)
    val = rnorm(1, mean=0, sd=(ub[j]-lb[j])/10)
    val = ifelse(r < mutpm, val, 0)
    y[j] = y[j] + val
  }
  return(list(mutant=y))
}

# Kod 2.91: Rastlantısal mutasyon 3 
randmut3 = function(y, lb, ub, mutpm, ...){
  m = length(y)
  if(missing(lb) | missing(ub))
    stop("randmut3 için lb ve ub gerekli")
  if(missing(mutpm)) mutpm=0.05
  for(j in 1:m){
    r = runif(1, 0, 1)
    val = rnorm(1, mean=0, sd=abs(lb[j]-ub[j]))
    val = ifelse(r < mutpm, val, 0)
    y[j] = y[j] + val
  }
  return(list(mutant=y))
}

# Kod 2.92: Rastlantısal mutasyon 4 
randmut4 = function(y, lb, ub,...){
  m = length(y)
  if(missing(lb) | missing(ub))
    stop("randmut4 için lb ve ub gerekli")
  v = sample(1:m, 1)
  r = runif(1, -0.1, 0.1)
  y[v] = y[v] + r * (ub[v]-lb[v])
  return(list(mutant=y, mutgen=v))
}

# Kod 2.93: Tekdüze değer mutasyonu
unimut = function(y, lb, ub, ...){
  if(missing(lb) | missing(ub)) 
    stop("unimut için lb ve ub gerekli")
  v = sample(1:length(y), 1)
  y[v] = runif(1, lb[v], ub[v])
  return(list(mutant=y, mutgen=v))
}

# Kod 2.94: Sınır değeri mutasyonu
boundmut = function(y, lb, ub, ...){
  if(missing(lb) | missing(ub)) 
    stop("boundmut için lb ve ub gerekli")
  v = sample(1:length(y), 1)
  r = runif(1, 0,1)
  y[v] = ifelse(r>0.5, ub[v], lb[v])
  return(list(mutant=y, mutgen=v))
}

# Kod 2.95: Tekdüze olmayan mutasyon
nunimut = function(y, lb, ub, g, gmax, mutb, ...){
  m = length(y)
  if(missing(lb) | missing(ub))
    stop("nunimut için lb ve ub gerekli!")
  if(missing(g))
    stop("nunimut için g gerekli!")
  if(missing(mutb))
    mutb = 0.5
  if(missing(gmax))
    stop("nunimut için gmax gerekli!")
  v = sample(1:m, 1)
  r = runif(1, 0, 1)
  if(r<=0.5)
    y[v] = y[v]+(ub[v]-y[v]) * r * (1-g/gmax)^mutb
  else
    y[v] = y[v]-(y[v]-lb[v]) * r * (1-g/gmax)^mutb
  return(list(mutant=y, mutgen=v))
}

# Kod 2.96: Uyarlamalı tekdüze olmayan mutasyon
nunimut2 = function(y, lb, ub, g, gmax, mutb, ...){
  m = length(y)
  if(missing(lb) | missing(ub))
    stop("nunimut için lb ve ub gerekli!")
  if(missing(g))
    stop("nunimut için g gerekli!")
  if(missing(mutb))
    mutb = 0.5
  if(missing(gmax))
    stop("nunimut için gmax gerekli!")
  v = sample(1:m, 1)
  r = runif(1, 0, 1)
  tau = sample(c(-1,1), 1)
  y[v] = y[v] + tau * (ub[v]-lb[v])*(1-r^(g/gmax)^mutb)
  return(list(mutant=y, mutgen=v))
}

# Kod 2.97: Güç mutasyonu
powmut = function(y, lb, ub, mutpow, ...){
  m = length(y)
  if(missing(lb) | missing(ub)) 
    stop("powmut için lb ve ub gerekli!")
  if(missing(mutpow)) mutpow = 2
  v = sample(1:m, 1)
  r = runif(1)
  p = runif(1)^mutpow
  t = (y[v]-lb[v])/(ub[v]-y[v])
  y[v] = y[v] + ifelse(r>t, -p*(y[v]-lb[v]), +p*(ub[v]-y[v]))
  return(list(mutant=y, mutgen=v))
}

# Kod 2.98: Güç mutasyonu 2
powmut2 = function(y, lb, ub, mutpow, ...){
  m = length(y)
  if(missing(lb) | missing(ub)) 
    stop("powmut için lb ve ub gerekli!")
  if(missing(mutpow)) mutpow = 2
  if(length(mutpow) < m) mutpow = rep(mutpow, m)
  v = sample(1:m, 1)
  r = runif(1)
  p = runif(1)^mutpow[v]
  t = (y[v]-lb[v])/(ub[v]-y[v])
  y[v] = y[v] + ifelse(r>t, -p*(y[v]-lb[v]), +p*(ub[v]-y[v]))
  return(list(mutant=y, mutgen=v))
}

# Kod 2.99: Gauss mutasyonu
gaussmut = function(y, mutsdy, ...){
  m = length(y)
  if(missing(mutsdy)) mutsdy=1
  v = sample(1:m, 1)
  rval = rnorm(1, 0, mutsdy)
  y[v] = y[v]+rval
  return(list(mutant=y, mutgen=v))
}

# Kod 2.100: Gauss mutasyonu 2
gaussmut2 = function(y, ...){
  m = length(y)
  v = sample(1:m, 1)
  y[v] = y[v]+rnorm(1)
  return(list(mutant=y, mutgen=v))
}

# Kod 2.101: Gauss mutasyonu 3
gaussmut3 = function(y, mutmy, mutsdy, ...){
  m = length(y)
  if(missing(mutmy) | missing(mutsdy)) return(0)
  val=NA
  while(is.na(val)){
    v = sample(1:m, 1)
    val = rnorm(1, mean=mutmy[v], sd=mutsdy[v])
  }
  y[v] = val
  return(list(mutant=y, mutgen=v))
}

# Kod 2.102: Sınırlarda aramaya özel mutasyon 1
bsearchmut1 = function(y, mutq, ...){
  m = length(y)
  if(missing(mutq)) mutq=runif(1)
  v = sample(1:m, 2, replace=FALSE)
  y[v[1]] = mutq*y[v[1]]
  y[v[2]] = mutq/y[v[2]]
  return(list(mutant=y, mutgen=v))
}

# Kod 2.103: Sınırlarda aramaya özel mutasyon 2
# Bağımlılık: Kod 2.103
bsearchmut2 = function(y, ...){
  m = length(y)
  v = sample(1:m, 2, replace=FALSE)
  p = runif(1)
  q = sqrt((y[v[1]]/ y[v[2]])^2 * (1-p^2) + 1)
  y[v[1]] = p*y[v[1]]
  y[v[2]] = q*y[v[2]]
  return(list(mutant=y, mutgen=v))
}

# Kod 2.104: Takas mutasyonu
swapmut = function(y, ...){
  n = length(y)
  v = sort(sample(1:n, size=2, replace=FALSE))
  takas = y[v[1]] 
  y[v[1]] = y[v[2]]
  y[v[2]] = takas
  return(list(mutant=y, mutgen=v))
}

# Kod 2.105: Ters çevirme mutasyonu
invmut = function(y, ...){
  n = length(y)
  v = sort(sample(1:n, size=2, replace=FALSE))
  subgenes = y[v[1]:v[2]] 
  subgenes = rev(subgenes)
  y[v[1]:v[2]] = subgenes
  return(list(mutant=y, mutrange=v))
}

# Kod 2.106: Karıştırma mutasyonu
shufmut = function(y, ...){
  n = length(y)
  v = sort(sample(1:n, size=2, replace=FALSE))
  subgenes = y[v[1]:v[2]] 
  subgenes = sample(subgenes)
  y[v[1]:v[2]] = subgenes
  return(list(mutant=y, mutrange=v))
}

# Kod 2.107: Araya girme mutasyonu
insmut = function(y, ...){
  n = length(y)
  v = sample(1:n, 1)
  k = sample(1:(n-1), 1)
  idx = c(setdiff(1:k,v), v, setdiff((k+1):n,v))
  y <- y[idx]
  return(list(mutant=y, mutgen=v, mutpoint=k))
}

# Kod 2.108: Yer değiştirme mutasyonu
dismut = function(y, ...){
  n = length(y)
  v = sort(sample(1:(n-1), size=2, replace=FALSE))
  r = sample(1:n, 1)
  while(r==v[1]) r = sample(1:n, 1)
  ytemp = rep(NA,n)
  yrem = setdiff(1:n, v[1]:v[2])
  idx = r
  for(yval in y[v[1]:v[2]]){
    if(idx>n) idx=1
    ytemp[idx]= yval
    idx =idx + 1
  }
  idx=1
  for(yval in yrem){
    while(!is.na(ytemp[idx]))
      idx=idx+1
    ytemp[idx]=yval
    idx=idx+1
  }
  y=ytemp
  return(list(mutant=y, mutrange=v, r=r))
}

# Kod 2.109: Takas + Ters çevirme mutasyonu
invswapmut = function(y, ...){
  n = length(y)
  v = sort(sample(1:n, size=2, replace=FALSE))
  takas = y[v[1]] 
  y[v[1]] = y[v[2]]
  y[v[2]] = takas
  v = sort(sample(1:n, size=2, replace=FALSE))
  subgenes = y[v[1]:v[2]] 
  subgenes = rev(subgenes)
  y[v[1]:v[2]] = subgenes
  return(list(mutant=y, mutgen=v))
}

# Kod 2.110: Araya girme + Ters çevirme mutasyonu
insswapmut = function(y, ...){
  n = length(y)
  v = sample(1:n, 1)
  k = sample(1:(n-1), 1)
  idx = c(setdiff(1:k,v), v, setdiff((k+1):n,v))
  y = y[idx]
  v = sort(sample(1:n, size=2, replace=FALSE))
  subgenes = y[v[1]:v[2]] 
  subgenes = rev(subgenes)
  y[v[1]:v[2]] = subgenes
  return(list(mutant=y, mutgen=v))
}

# Kod 2.111: Yer değiştirme + Ters çevirme mutasyonu
invdismut = function(y, ...){
  n = length(y)
  v = sort(sample(1:(n-1), size=2, replace=FALSE))
  r = sample(1:n, 1)
  while(r==v[1]) r = sample(1:n, 1)
  ytemp = rep(NA,n)
  yrem = setdiff(1:n, v[1]:v[2])
  yrevpart = rev(y[v[1]:v[2]])
  idx = r
  for(yval in yrevpart){
    if(idx>n) idx=1
    ytemp[idx]= yval
    idx =idx + 1
  }
  idx=1
  for(yval in yrem){
    while(!is.na(ytemp[idx]))
      idx=idx+1
    ytemp[idx]=yval
    idx=idx+1
  }
  y=ytemp
  return(list(mutant=y, mutrange=v, r=r))
}

# Kod 2.112: Mutasyon uygulama fonksiyonu
mutate = function(mutfunc, population, mutpm, gatype, ...){
  func = as.character(match.call()[2])
  dotargs = list(...)
  n = nrow(population)
  m = ncol(population)-2
  if(missing(gatype)) gatype="gga"
  if(missing(mutpm)) mutpm=0.05
  nm = round(n*mutpm)  # Mutasyon uygulanan birey sayısı
  if(gatype=="ssga"){
    if(runif(1,0,1) <= mutpm) nm=1 else nm=0
  }
  if(nm>0){
    midx = sample(1:n, size=nm, replace=FALSE) # Mutant indeksleri
    for(i in midx){
       dotargs$y = unname(population[i,1:m])
       population[i,1:m] = do.call(func, dotargs)$mutant
    }
  }
  return(population)
}

# Kod 2.113: Tümünü silerek yenileme
grdelall = function(parpop, offpop){
  parpop = offpop
  return(parpop)
}

# Kod 2.114: Seçkinci yenileme (elitizm) fonksiyonu
elitism = function(parpop, offpop, reps, ...){
  if(missing(reps)) 
    reps=max(1, round(nrow(parpop)*0.05)) # En iyi birey sayısı
  parpop = parpop[order(parpop[,ncol(parpop)], decreasing=TRUE),]
  offpop = offpop[order(offpop[,ncol(offpop)], decreasing=TRUE),]
  offpop[(nrow(offpop)-reps+1):nrow(offpop),] = parpop[1:reps,]
  parrn = rownames(parpop[1:nrow(parpop),])
  offrn = rownames(offpop[1:(nrow(offpop)-reps),])
  unirn = c(offrn, parrn[1:reps])
  rownames(offpop) = unirn
  return(offpop)
} 

# Kod 2.115: Mu+Lamda yenileme fonksiyonu 1
grmuplambda = function(parpop, offpop, ...){
  n = nrow(parpop)
  unipop = rbind(parpop, offpop)
  unipop = unipop[order(unipop[, ncol(unipop)], decreasing=TRUE),]
  unipop = unipop[1:n,]
  return(unipop)
}

# Kod 2.116: Mu+Lamda yenileme fonksiyonu 2 (en kötü λ silme)
grmuplambda2 = function(parpop, offpop, lambda, ...){
  if(missing(lambda)) lambda = round(nrow(offpop)/2)
  n = nrow(offpop)
  parpop = parpop[order(parpop[, "fitval"], decreasing=TRUE),]
  offpop = offpop[order(offpop[, "fitval"], decreasing=TRUE),]
  parpop[(n-lambda+1):n,] = offpop[1:lambda,]
  return(parpop)
}

# Kod 2.117: Mu+Lamda yenileme fonksiyonu 3 
grmuplambda3 = function(parpop, offpop, lambda, ...){
  if(missing(lambda)) lambda = round(nrow(offpop)/2)
  n = nrow(offpop)
  paridx = sample(1:n, size=lambda, replace=FALSE)
  offpop = offpop[order(offpop[, "fitval"], decreasing=TRUE),]
  parpop[paridx,] = offpop[1:lambda,]
  return(parpop)
}

# Kod 2.118: Mu+Lamda yenileme fonksiyonu 4 
grmuplambda4 = function(parpop, offpop, lambda, ...){
  if(missing(lambda)) lambda = round(nrow(offpop)/2)
  npar = nrow(parpop)
  noff = nrow(offpop)
  paridx = sample(1:npar, size=lambda, replace=FALSE)
  offidx = sample(1:noff, size=lambda, replace=FALSE)
  parpop[paridx,] = offpop[offidx,]
  return(parpop)
}

# Kod 2.119: Mu & Lamda yenileme fonksiyonu
grmuvlambda = function(parpop, offpop, ...){
  n = nrow(parpop)
  offpop = offpop[order(offpop[, "fitval"], decreasing=TRUE),]
  offpop = offpop[1:n,]
  return(offpop)
}

# Kod 2.120: Round Robin yenileme fonksiyonu
grrobin = function(parpop, offpop, repk, ...){
  if(missing(repk)) repk = 10
  n = nrow(parpop)
  m = ncol(parpop)
  unipop = rbind(parpop, offpop)
  wins = rep(0, nrow(unipop))
  unipop = cbind(unipop, wins)
  idx = sample(1:nrow(unipop), replace=FALSE)
  for(i in idx){
    tidx = idx[-which(i==idx)]
    touridx = sample(tidx, size=repk, replace=FALSE)
    for(j in touridx){
      if(unipop[i,m] > unipop[j,m]){ 
        unipop[i,"wins"] = unipop[i,"wins"]+1
      }
    }
  }
  unipop = unipop[order(unipop[,"wins"], decreasing=TRUE),]
  unipop = unipop[1:n,1:m]
  return(unipop)
}

# Kod 2.121: Mu+1 yenileme fonksiyonu
ssrmup1 = function(parpop, offpop, ...){
  n = nrow(parpop)
  m = ncol(parpop)
  parrn = rownames(parpop)
  offrn = rownames(offpop)
  if(nrow(offpop) > 1)
    offpop = offpop[order(offpop[,m], decreasing=TRUE),]
  idx = sample(1:n, size=1)
  if(offpop[1,m] > parpop[idx,m]){
    parpop[idx,] = offpop[1,]
    parrn[idx] = offrn[1]
    rownames(parpop) = parrn
  }
  return(parpop)
}

# Kod 2.122: Genitor yenileme fonksiyonu
ssrgenitor = function(parpop, offpop, ...){
  n = nrow(parpop)
  m = ncol(parpop)
  parrn = rownames(parpop)
  offrn = rownames(offpop)
  parpop = parpop[order(parpop[,m], decreasing=TRUE),]
  if(nrow(offpop) > 1)
    offpop = offpop[order(offpop[,m], decreasing=TRUE),]
  parpop[n,] = offpop[1,]
  parrn[n] = offrn[1]
  rownames(parpop) = parrn
  return(parpop)
}

# Kod 2.123: Aile turnuvası ile yenileme fonksiyonu
ssrfamtour = function(parpop, offpop, reppars, ...){
  if(missing(reppars)) reppars = c(1:2)
  n = nrow(parpop)
  m = ncol(parpop)
  parrn = rownames(parpop)
  offrn = rownames(offpop)
  fampop = rbind(parpop[reppars,], offpop)
  fampop = fampop[order(fampop[,m], decreasing=TRUE),]
  famrn = rownames(fampop)
  parpop[reppars,] = fampop[1:2,]
  parrn[reppars] = famrn[1:2]
  rownames(parpop)=parrn
  return(parpop)
}

# Kod 2.124: Karma yenileme fonksiyonu
ssrx = function(parpop, offpop, reppars, ...){
  n = nrow(parpop)
  m = ncol(parpop)
  if(nrow(offpop) > 1)
    offpop = offpop[order(offpop[,m], decreasing=TRUE),]
  parrn = rownames(parpop)
  offrn = rownames(offpop)
  bidx = which.max(parpop[,m])
  paridx = 1:n
  paridx = paridx[-c(reppars, bidx)]
  pidx = sample(paridx, size=1)
  parpop[pidx,] = offpop[1,]
  parrn[pidx] = offrn[1]
  rownames(parpop) = parrn
  return(parpop)
}

# Kod 2.125: Sonlandırma kontrolü fonksiyonu
terminate = function(tercrit, maxiter, objective, t, genfits,
    fitvals, objval, optdif, rmcnt, rmdif, abdif, mincv,
    sddif, rangedif, simlev, phidif, meandif, bestdif,
    stime, maxtime){
  if(missing(optdif)) optdif=1e-06
  if(missing(rmcnt)) rmcnt=10
  if(missing(rmdif)) rmdif=1e-06
  if(missing(abdif)) abdif=1e-06
  if(missing(phidif)) phidif=1e-04
  if(missing(meandif)) meandif=1e-06
  if(missing(bestdif)) bestdif=1e-03
  if(missing(sddif)) sddif=1e-03
  if(missing(rangedif)) rangedif=1e-03
  if(missing(mincv)) mincv=0.05
  if(missing(simlev)) simlev=0.975
  if(missing(maxtime)) 
    maxtime=ifelse(maxiter>=1e06, 3e-04*maxiter, 60)
  terminate = 0
 # Maksimum iterasyon sayısı
  if(is.element(1, tercrit)){
    if(t > maxiter) terminate = 1
  }
 #Global optimuma (amaç değere) ulaşma
  if(is.element(2,tercrit)){
    if(!missing(objval)){
      if(genfits[t,2]==objval) terminate = 2
    }
  }
 # Global optimuma ulaşma veya çok yaklaşma
  if(is.element(3,tercrit)){
    if(!missing(objval)){
      optdifc = abs(objval-genfits[t,2])
      if(optdifc <= optdif) terminate = 3
    }
  }
 # Son k en iyi uyum ortalaması arası minimum fark
  if(is.element(4, tercrit)){
    if(t > rmcnt){
      tb = t-rmcnt+1 
      rmean = mean(genfits[tb:t,2])
      cbestdif = abs(rmean-genfits[t,2])
    }else{
      cbestdif = Inf
    }
    if(cbestdif <= bestdif) terminate = 4
  }
 # Son iki uyum ortalaması arasındaki minimum fark
  if(is.element(5, tercrit)){
    cmeandif = Inf
    if(t > 1) 
      cmeandif = abs(genfits[t,3]-genfits[(t-1),3])
    if(cmeandif <= meandif) terminate = 5
  }
 # Son k uyum ortalaması arası minimum fark
  if(is.element(6, tercrit)){
    if(t > rmcnt){
      tb = t-rmcnt+1 
      rmean = mean(genfits[tb:t,3])
      crmdif = abs(rmean-genfits[t,3])
    }else{
      crmdif = Inf
    }
    if(crmdif <= rmdif) terminate = 6
  }
 # Ortalama ve en iyi uyum arası minimum fark
  if(is.element(7, tercrit)){
    conv = abs(genfits[t,2]-genfits[t,3])
    if(conv <= abdif) terminate = 7
  }
 # Minimum standart sapma farkı
  if(is.element(8, tercrit)){
    if(t > 1){
      csddif = abs(genfits[t,4] - genfits[t-1,4])
      if(csddif <= sddif) terminate = 8
    }
  }
 # Minimum ve maksimum farkı
  if(is.element(9, tercrit)){
    crangedif = abs(genfits[t,1]-genfits[t,2])
    if(crangedif <= rangedif) terminate = 9
  }
 # Minimum varyasyon katsayısı
  if(is.element(10, tercrit)){
    cv = ifelse(genfits[t,3]==0, 
      Inf, genfits[t,4]/genfits[t,3])
    if(abs(cv) <= mincv) terminate = 10
  }
 # Phi yakınsaması
  if(is.element(11, tercrit)){
    phi = ifelse(genfits[t,2]==0, 
      1, 1-abs(genfits[t,3]/genfits[t,2]))
    if(phi <= phidif) terminate = 11
  }
 # Uyum benzerlik yüzdesi
  if(is.element(12, tercrit)){
    csimlev = (length(fitvals)-length(unique(fitvals))) /
    length(fitvals)
    if(csimlev >= simlev) terminate = 12
  }
 # Zaman aşımı
  if(13 %in% tercrit){
    timedif = as.difftime(Sys.time()-stime, units="mins")
    timedif = as.numeric(timedif, units = "mins")
    if(timedif > maxtime) terminate = 13
  }
  return(terminate)
}

# Kod 2.126: Problem 1 için cezalı uyum fonksiyonu 
fitnessp11 = function(x, r, c, beta, gamma, objective,...){
  if(missing(r)) r=sqrt(.Machine$double.xmax)
  if(missing(c)) c=sqrt(.Machine$double.xmax)
  if(missing(beta)) beta=2
  if(missing(gamma)) gamma=1
  if(missing(objective)) objective="max"
  s = ifelse(objective=="max", -1, 1)
  g1 = function(x) (-x[1]- x[2]+4)
  fobj = function(x) (2*x[1]^2+9*x[2])
  G=c()
  G[1] = max(0, g1(x)^beta)
  fitval = fobj(x) + s*r*sum(G)
  return(fitval)
}

# Kod 2.127: Problem 1 için iç cezalı uyum fonksiyonu 
fitnessp12 = function(x, r, objective,...){
  if(missing(objective)) objective="max"
  if(missing(r)) r=0.1
  if(objective=="max") s=-1 else s=1
  fobj = function(x) (2*x[1]^2+9*x[2])
  g1 = function(x) (-x[1]-x[2]+4)
  G = c()
  repeat{
    if(g1(x)==0)
      G[1] = 0
    else
      G[1] = 1/r*(-1/g1(x)) 
    r=r*10
    if(G[1]<=0) break
  }
  fitval = s*(fobj(x) + sum(G))
  return(fitval)
}

# Kod 2.128: Problem 1 için logaritmik cezalı uyum fonksiyonu 
fitnessp13 = function(x, r, objective,...){
  if(missing(objective)) objective="max"
  if(missing(r)) r=0.1
  if(objective=="max") s=-1 else s=1
  fobj = function(x) (2*x[1]^2+9*x[2])
  g1 = function(x) (-x[1]-x[2]+4)
  G = c()
  repeat{
    if(g1(x)==0)
      G[1] = 0
    else
      G[1] = -1/r*log2(-g1(x)) 
    r=r*10
    if(G[1]<=0) break
  }
  fitval = s*(fobj(x) + sum(G))
  return(fitval)
}

# Kod 2.129: Problem 1 için genişletilmiş cezalı uyum fonksiyonu 
fitnessp14 = function(x, r, epsilon, objective,...){
  if(missing(objective)) objective="max"
  if(missing(r)) r=0.1
  if(missing(epsilon)) epsilon=-0.1
  if(objective=="max") s=-1 else s=1
  fobj = function(x) (2*x[1]^2+9*x[2])
  g1 = function(x) (-x[1]-x[2]+4)
  G = c()
  repeat{
    if(g1(x)==0) {G[1]=0; break}
    if(g1(x)<=epsilon)
      G[1] = 1/r*(-1/g1(x))
    else
      G[1] = 1/r*(1/epsilon^2 * (2*epsilon-g1(x)))
    r=r*10
    if(G[1]<=0) break
  }
  fitval = s*(fobj(x) + sum(G))
  return(fitval)
}

# Kod 2.130: Problem 1 için statik cezalı uyum fonksiyonu 
fitnessp15 = function(x, r, c, beta, gamma, objective,...){
  if(missing(r)) r=sqrt(.Machine$double.xmax)
  if(missing(c)) c=sqrt(.Machine$double.xmax)
  if(missing(beta)) beta=2
  if(missing(gamma)) gamma=1
  if(missing(objective)) objective="max"
  s = ifelse(objective=="max", -1, 1)
  g1 = function(x) (-x[1]- x[2]+4)
  G = c()
  repeat{
    G[1] = max(0, r*g1(x))^beta 
    if(G[1]<=0) break
    r=r/1e06
  }
  fitval = fobj(x) + s*sum(G)
  return(fitval)
}

# Kod 2.131: Problem 2 için cezalı uyum fonksiyonu 
fitnessp2 = function(x, r, c, beta, gamma, objective,...){
  if(missing(r)) r=sqrt(.Machine$double.xmax)
  if(missing(c)) c=sqrt(.Machine$double.xmax)
  if(missing(beta)) beta=2
  if(missing(gamma)) gamma=1
  if(missing(objective)) objective="max"
  s = ifelse(objective=="max", -1, 1)
  g1 = function(x) (-x[1]^2+x[2])
  g2 = function(x) (-x[2])
  h1 = function(x) (x[1]*x[2]^2-1)
  fobj = function(x) (1/2*x[1]^2+x[1]*x[2]^2)
  G=c(); H=c()
  G[1] = max(0, g1(x))^beta
  G[2] = max(0, g2(x))^beta
  H[1] = max(0, abs(h1(x))^gamma)
  fitval = fobj(x) + s*(r*sum(G)+ c*sum(H))
  return(fitval)
}

# Kod 2.132: Problem 3 için cezalı uyum fonksiyonu 
fitnessp3 = function(x, r, c, beta, gamma, objective,...){
  if(missing(r)) r=sqrt(.Machine$double.xmax)
  if(missing(c)) c=sqrt(.Machine$double.xmax)
  if(missing(beta)) beta=2
  if(missing(gamma)) gamma=1
  if(missing(objective)) objective="max"
  s = ifelse(objective=="max", -1, 1)
  g1 = function(x) (-3*x[1]-2*x[2]+6)
  g2 = function(x) (-x[1]+x[2]-3)
  g3 = function(x) (x[1]*x[2]-7)
  g4 = function(x) (2/3*x[1]-x[2]-4/3)
  fobj = function(x) ((x[1]-6)^2+(x[2]-7)^2)
  G=c()
  G[1] = max(0, g1(x)^beta)
  G[2] = max(0, g2(x)^beta)
  G[3] = max(0, g2(x)^beta)
  G[4] = max(0, g2(x)^beta)
  fitval = fobj(x) + s*r*sum(G)
  return(fitval)
}

# Kod 2.133: Statik çaprazlama ve mutasyon oranı
fixpcmut = function(cxpc, mutpm, ...){
  return(list(pc=cxpc, pm=mutpm))
}

# Kod 2.134: ILM/DHC uyarlama fonksiyonu 
ilmdhc = function(g, gmax, ...){
  pmg = g/gmax
  pcg = 1-pmg
  return(list(pc=pcg, pm=pmg))
}

# Kod 2.135: Uyarlanabilir Dinamik Algoritma (Adana 1)
adana1 = function(g, gmax, ...){
  pm = sin((g/gmax))+0.01
  pc = 1-pm
  return(list(pc=pc, pm=pm))
}

# Kod 2.136: Uyarlanabilir Dinamik Algoritma (Adana 2)
adana2 = function(g, gmax, ...){
  pm = sqrt((g/gmax))
  pc = 1-pm
  return(list(pc=pc, pm=pm))
}

# Kod 2.137: Lei & Tingzhi uyarlama fonksiyonu
leitingzhi = function(fitvals, cxpc, cxpc2, 
  mutpm, mutpm2, adapa, adapb, ...){
  if(missing(adapa)) adapa=0.7
  if(missing(adapb)) adapb=0.5
  if(missing(cxpc)) cxpc=0.9
  if(missing(cxpc2)) cxpc2=0.5
  if(missing(mutpm)) mutpm=0.05
  if(missing(mutpm2)) mutpm2=0.2
  fmmr = min(fitvals)/max(fitvals)
  famr = mean(fitvals)/max(fitvals)
  pc = ifelse(famr>adapa & fmmr>adapb, cxpc2*(1-fmmr^-1), cxpc)
  pm = ifelse(famr>adapa & fmmr>adapb, mutpm2*(1-fmmr^-1), mutpm)
  return(list(pc=pc, pm=pm))
}

# Kod 2.138: Dinamik mutasyon ve çaprazlama fonksiyonu (Adana 3)
adana3 = function(fitvals, g, gmax, cxpc, mutpm,
  adapc, adapd, ...){
  if(missing(adapc)) adapc = 0.05
  if(missing(adapd)) adapd = 0.05
  fmin = min(fitvals)
  fmax = max(fitvals)
  favg = mean(fitvals)
  newrate = ifelse(g/gmax>=0.01, g/gmax,0.01)
  pca = ifelse((fmax-favg)/fmax <= adapc, newrate, cxpc)
  pma = ifelse((fmax-favg)/fmax <= adapd, newrate, mutpm)
  return(list(pc=pca, pm=pma))
}

# Kod 2.139: Basit iki değişkenli bir uyum fonksiyonu
fxy = function(x) 2*(x[1]-1)^2 + 5*(x[2]-3)^2 + 10

# Kod 2.140: optim'de değiştirilecek argümanlar listesi 
hgaparams = list(method="Nelder-Mead", poptim=0.05, pressel=0.5,
  control = list(fnscale=1, maxit=100))

# Kod 2.141: GA + optim melezleme fonksiyonu
hgaoptim = function(genpop, fitfunc, hgaparams,
  hgaftype, hgans, ...){
  n = nrow(genpop)
  m = ncol(genpop)
  if(missing(hgans)) hgans = 1
  if(hgans > n) hgans = n
  if(missing(hgaftype)) hgaftype="w" 
  if(hgaftype=="b"){ 
    selpop = genpop[order(genpop[,m], decreasing=TRUE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="w"){
    selpop = genpop[order(genpop[,m], decreasing=FALSE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="r"){
    selpop = genpop[sample(1:n, size=hgans, replace=FALSE),]
  }
  newpop = matrix(NA, nrow=hgans, ncol=m)
  rnselpop = rownames(selpop)
  for(i in rnselpop){
    sol = stats::optim( 
      fn = fitfunc,
      par = genpop[i,1:(m-2)],
      method = hgaparams$method,
      control = hgaparams$control)
    genpop[i,1: (m-2)] = sol$par
    genpop[i,"fitval"] = sol$value
  }
  return(genpop)
}

# Kod 2.142: GA+optimx melezleme fonksiyonu
hgaoptimx = function(genpop, fitfunc, hgaparams,
  hgaftype, hgans, ...){
  n = nrow(genpop)
  m = ncol(genpop)
  if(missing(hgans)) hgans = 1
  if(hgans > n) hgans = n
  if(missing(hgaftype)) hgaftype="w" 
  if(hgaftype=="b"){ 
    selpop = genpop[order(genpop[,m], decreasing=TRUE),]
    selpop = selpop[1:ns,]
  }else if(hgaftype=="w"){
    selpop = genpop[order(genpop[,m], decreasing=FALSE),]
    selpop = selpop[1:ns,]
  }else if(hgaftype=="r"){
    selpop = genpop[sample(1:n, size=hgans, replace=FALSE),]
  }
  newpop = matrix(NA, nrow=hgans, ncol=m)
  rnselpop = rownames(selpop)
  for(i in rnselpop){
    sol = optimx::optimx( 
      fn = fitfunc,
      par = genpop[i,1:(m-2)],
      method = hgaparams$method,
      lower = hgaparams$lower,
      upper = hgaparams$upper,
      control = hgaparams$control)
    sol = c(unname(unlist(sol)))
    genpop[i,1:(m-2)] = sol[1:(m-2)]
    genpop[i,m] = sol[m-2+1]
  }
  return(genpop)
}

# Kod 2.143: optim'de değiştirilecek argümanlar listesi 
lb=rep(-5.12, 2); ub=rep(5.12,2)
hgaparams = list(method="L-BFGS-B", 
  poptim=0.05, pressel=0.5,
  lower=lb, upper=ub,
  control=list(maximize=FALSE, maxit=100))

# Kod 2.144: GA+ROI melezleme fonksiyonu
hgaroi = function(genpop, fitfunc, hgaparams,
  hgaftype, hgans, ...){
  n = nrow(genpop)
  m = ncol(genpop)
  if(missing(hgans)) hgans = 1
  if(hgans > n) hgans = n
  if(missing(hgaftype)) hgaftype="w" 
  if(hgaftype=="b"){ 
    selpop = genpop[order(genpop[,m], decreasing=TRUE),]
    selpop = selpop[1:ns,]
  }else if(hgaftype=="w"){
    selpop = genpop[order(genpop[,m], decreasing=FALSE),]
    selpop = selpop[1:ns,]
  }else if(hgaftype=="r"){
    selpop = genpop[sample(1:n, size=hgans, replace=FALSE),]
  }
  newpop = matrix(NA, nrow=hgans, ncol=m)
  rnselpop = rownames(selpop)
  require(ROI)
  fo = F_objective(F=fitfunc, n=(m-2))
  vb = V_bound(li=1: (m-2), ui=1: (m-2), 
    lb=hgaparams$lower, ub=hgaparams$upper)
  op = OP(objective=fo, bounds = vb, maximum=TRUE)
  for(i in rnselpop){
    sol = ROI_solve(op, 
      solver="optimx", 
      control=hgaparams$control,
      start=genpop[i,1:(m-2)])
    genpop[i,1:(m-2)] = sol$solution
    genpop[i,m] = sol$objval
  }
  return(genpop)
}

# Kod 2.145: Uyum değeri değişimini görüntüleme
monprogress = function(g, genfits, objective, ...){
  if(missing(g)) stop("kuşak numarası, g eksik")
  if(missing(genfits)) stop("kuşak uyum değerleri, genfits eksik")
  if(missing(objective)) stop("objective eksik")
  if(objective=="min") genfits = -1*genfits
  plot(genfits[1:g,3], col=2, lwd=2, type="l",
    ylab="Uyum", xlab="Kuşak")
  title(main=paste("İterasyon", g))
}

# Kod 2.146: İterasyon sonuçlarını görselleştirme fonksiyonu
show = function(monitorfunc, g, genfits, objective, x, ...){
  args = as.list(match.call())[-1]
  nargs = length(args)
  if(missing(monitorfunc)) monitorfunc=monprogress
  if(missing(g)) stop("Parametre g eksik")
  if(missing(genfits)) stop("Parametre genfits eksik")
  if(missing(objective)) stop("Parametre objective eksik")
  if(missing(x)) stop("Parametre x eksik")
  do.call(as.character(args[1]), args[2:nargs])
}

# Kod 3.1: Genel Amaçlı Uyarlamalı GA Optimizasyon Fonksiyonu
adana = function(gatype="gga", objective="max", maxiter=100, 
  initfunc=initbin, fitfunc, selfunc=seltour,
  crossfunc=px1, mutfunc=bitmut, replace=elitism,
  adapfunc=NULL, hgafunc=NULL, monitorfunc=NULL,
  n=100, m=8, lb=rep(0,8), ub=rep(1,8), nmode="real", type=1,
  permset=0:9, prevpop=NULL, 
  selt=2, selbc=0.5, selc=2, selk=1.005, sells=1.5, 
  selns=0.5, selps=0.5, sels=1.5, selt0=50, selw=2, 
  reptype=FALSE, cxpc=0.9, cxpc2=0.8, cxon=2, cxk=2, 
  cxps=0.5, cxa=0, cxb=0.15, cxealfa=1, cxalfa=0.5,
  mutpm=0.05, mutpm2=0.2, mutb=2, mutpow=2, mutq=0.5,
  mutmy=c(), mutsdy=c(), adapa=0.75, adapb=0.5,
  adapc=0.1, adapd=0.1, hgastep=10, hgans=1, hgaftype="w",
  reps=1, repk=10, reppar=c(1:2), lambda=1, tercrit=c(1), 
  abdif=1e-06, bestdif=1e-06, objval=0, optdif=1e-06,
  rmcnt=10, rmdif=1e-06, phidif=1e-06, rangedif=1e-06,
  meandif=1e-06, sddif=1e-06, mincv=0.001, simlev=0.95,
  maxtime=60, keepbest = TRUE, parcontrol=NULL,
  parfile=NULL, verbose=FALSE, ...){
 # Varsa parametreler dosyasını yükle
  if(!is.null(parfile)) source(parfile)
 # Başlangıç popülasyonu
  initpop=initialize(initfunc, n=n, m=m, lb=lb, ub=ub,
    permset=permset, nmode=nmode, prevpop=prevpop, type=1)
  genpop = initpop
  genfits = matrix(NA, nrow=maxiter, ncol=8)
  colnames(genfits) = c("min","max", "mean", "sd",
    "Q1", "Q2", "Q3", "cv%")
  stime = Sys.time()
  g = 1
  tcode = 0
  bestsol = list(chromosome=c(), fitval=-Inf, generation=0)
  while(tcode==0 & g<=maxiter){
    genpop[,m+2] = evaluate(fitfunc, genpop[,1:m], objective)
    genfits[g,1] = min(genpop[,m+2])
    genfits[g,2] = max(genpop[,m+2])
    genfits[g,3] = mean(genpop[,m+2])
    genfits[g,4] = sd(genpop[,m+2])
    genfits[g,5] = quantile(genpop[,m+2])[2]
    genfits[g,6] = quantile(genpop[,m+2])[3]
    genfits[g,7] = quantile(genpop[,m+2])[4]
    genfits[g,8] = round(genfits[g,4]/genfits[g,3]*100,2)
    if(!missing(tercrit))
      tcode = terminate(tercrit=tercrit, maxiter=maxiter,
        objective=objective, t=g, genfits=genfits, 
        fitvals=genpop[,m+2], objval=objval, optdif=optdif,
        rmcnt=rmcnt, rmdif=rmdif, abdif=abdif, mincv=mincv,
        sddif=sddif, rangedif=rangedif, phidif=phidif,
        simlev=simlev, meandif=meandif, bestdif=bestdif,
        stime=stime, maxtime=maxtime)
    fitvals = genpop[,ncol(genpop)]
    if(keepbest){
      if(bestsol$fitval < max(fitvals)){
        bidx = which.max(fitvals)[1]
        bestsol$chromosome = genpop[bidx, 1:m]
        bestsol$fitval = genpop[bidx, m+2]
        bestsol$generation = g
      }
    }
    if(!is.null(adapfunc)){  #Uyarlama
       ocxpc = cxpc; omutpm = mutpm
       adappar = adapfunc(fitvals=fitvals, n=n, g=g, 
         gmax=maxiter, cxpc=cxpc, cxpc2=cxpc2, mutpm=mutpm,
         mutpm2=mutpm2, adapa=adapa, adapb=adapb,
         adapc=adapc, adapd=adapd)
       cxpc = adappar$pc ; mutpm = adappar$pm 
    }
    selidx = select(selfunc, fitvals, selt=selt)
    matpool = genpop[selidx,]
    shfidx = sample(1:nrow(matpool)) 
    matpool = matpool[shfidx,] # Çiftleşme havuzunu karıştır 
    reppars = shfidx[1:2] # SSGA için ilk iki bireyi işaretle 
    offsprings = cross(crossfunc, matpool=matpool,
      cxon=cxon, cxpc=cxpc, lb=lb, ub=ub,
      gatype=gatype, cxk=cxk, cxps=cxps, cxa=cxa, cxb=cxb)
    if(is.vector(offsprings)) 
      offsprings = matrix(offsprings, nrow=1, ncol=m)
    mutoffsprings = mutate(mutfunc=mutfunc, population=offsprings, 
      mutpm=mutpm, gatype=gatype, lb=lb, ub=ub, 
      mutmy=mutmy, mutsdy=mutsdy, mutb=mutb, mutpow=mutpow,
      g=g, gmax=maxiter)
    if(is.vector(mutoffsprings))
      mutoffsprings = matrix(mutoffsprings, nrow=1, ncol=m) 
    rownames(mutoffsprings) = paste0("T",g, ".",
       1:nrow(mutoffsprings))
    mutoffsprings[, m+1] = -1
    mutoffsprings[, m+2] = evaluate(fitfunc, 
       mutoffsprings[,1:m], objective)
    genpop = replace(parpop=genpop, offpop=mutoffsprings, 
       repk=repk, reps=reps, reppars=reppars)
    genpop[, (m+1)] = genpop[, (m+1)]+1
    if(!is.null(hgafunc) & g%%hgastep==0) #Melezleme
      genpop = hgafunc(genpop, fitfunc=fitfunc,
      hgaparams=hgaparams, hgaftype=hgaftype, hgans=hgans)
    if(verbose) cat("Kuşak ", g,":", genfits[g,], "\n")
    if(g>1 & !is.null(monitorfunc)){
      x = bestsol$chromosome
      show(monitorfunc, g=g, genfits=genfits, 
        objective=objective, x=x)
    }
    if(!is.null(adapfunc)){
      cxpc = ocxpc
      mutpm = omutpm
    }
    g = g+1
  }
  genfits = genfits[!is.na(genfits[,1]),]
  if(is.element(tcode, c(4,7)))
    genfits=genfits[1:nrow(genfits)-1,]
  if(objective=="min"){
    if(!is.matrix(genfits)) genfits = matrix(genfits, nr=1,nc=8)
    genfits[,c(1:3,5:8)] = -1*genfits[,c(1:3,5:8)]
    bestsol$fitval = -1*bestsol$fitval
  }
  #cat("Termination criterion:", tcode,"\n")
  return(list(genfits=genfits, initpop=initpop, finalpop=genpop, 
    bestsol=bestsol, objective=objective, tcode=tcode)) 
}

# Kod 3.2a: KSP uyum fonksiyonu (Ölüm cezası ile)
kspfit1 = function(x, ...) {
  tpuan = x %*% kspveri$puan
  tagirlik = x %*% kspveri$agirlik
  fitval = ifelse(tagirlik > kapasite, 0, tpuan) 
  return(fitval)
}

# Kod 3.2b: KSP uyum fonksiyonu 2 (dış ceza ile)
kspfit2 = function(x, ...) {
  tpuan = x %*% kspveri$puan
  tagirlik = x %*% kspveri$agirlik
  G1 = tagirlik-kapasite #Ceza fonksiyonu
  fitval = tpuan-max(0,G1)^2 #Modifiye uyum fonksiyonu
  return(fitval)
}

# Kod 3.3: En iyi çözüm görüntüleme fonksiyonu
bestsol = function(garesult){
  chromosome = garesult$bestsol$chromosome
  fitval = garesult$bestsol$fitval
  generation = garesult$bestsol$generation
  cat("Best solution\n")
  cat(" Chromosome:", chromosome,"\n")
  cat(" Fitness:", fitval,"\n")
  cat(" Generation:", generation,"\n")
}

# Kod 3.4: GA kuşaklara göre uyum istatistikleri grafiği
plotfitness = function(genfits, options){
  if(missing(options)) options=c(3,2)
  if(is.vector(genfits))
    genfits = matrix(genfits, nrow=1, ncol=8)
  x1 = x2 = x3 = x4 = NULL
  lbls = c("min","maks","ort","Q1","med","Q3")
  clrs = c("gray","blue","red", "pink","cyan","green")
  ltys = c(3,2,1,4,5,6)
  x1 = genfits[,options[1]]
  if(length(options)==1)
    ltys[options] = 1
  if(length(options)==2)
    x2=genfits[,options[2]]
  else if(length(options)==3){
    x2=genfits[,options[2]]
    x3=genfits[,options[3]]
  }else if(length(options)==4){
    x2=genfits[,options[2]]
    x3=genfits[,options[3]]
    x4=genfits[,options[4]]
  }
  plot(x1, col=clrs[options[1]], type="l", 
    lwd=5, lty=ltys[options[1]], 
  ylim=c(min(genfits[,1]), max(genfits[,2])),
  xlab="İterasyonlar", ylab="Uyum Değeri",
  main="İterasyonlara Göre Uyum Değerleri")
  lines(x1, col=7, lwd=1, lty=1)
  if(!is.null(x2))
    lines(x2, col=clrs[options[2]], lwd=2, lty=ltys[options[2]])
  if(!is.null(x3))
    lines(x3, col=clrs[options[3]], lwd=2, lty=ltys[options[3]])
  if(!is.null(x4))
    lines(x4, col=clrs[options[4]], lwd=2, lty=ltys[options[4]])
  legend("bottom", inset=.02, 
  lbls[options], col=clrs[options], 
  lty=ltys[options], horiz=TRUE, cex=0.8)
}

# Kod 3.5: adana+optim melezleme fonksiyonu
hgaoptim = function(genpop, fitfunc, hgaparams,
  hgaftype, hgans, ...){
  n = nrow(genpop)
  m = ncol(genpop)
  if(missing(hgans)) hgans = 1
  if(hgans > n) hgans = n
  if(missing(hgaftype)) hgaftype="w" 
  if(hgaftype=="b"){ 
    selpop = genpop[order(genpop[,m], decreasing=TRUE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="w"){
    selpop = genpop[order(genpop[,m], decreasing=FALSE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="r"){
    selpop = genpop[sample(1:n, size=hgans, replace=FALSE),]
  }
  newpop = matrix(NA, nrow=hgans, ncol=m)
  rnselpop = rownames(selpop)
  fobj = function(x){
      (decode(fitfunc(encode(x,
        lb=hgaparams$lower, ub=hgaparams$upper, m=m-2)),
      lb=hgaparams$lower, ub=hgaparams$upper, m=m-2))
  }
  for(i in rnselpop){
    decpar = decode(genpop[i,1: (m-2)], 
      lb=hgaparams$lower, ub=hgaparams$upper, m=m-2)
    sol = stats::optim( 
      fn = fobj,
      par = decpar,
      method = hgaparams$method,
      lower = hgaparams$lower,
      upper = hgaparams$upper,
      control = hgaparams$control)
    encpar = encode(sol$par,
      lb=hgaparams$lower, ub=hgaparams$upper, m=m-2)
    genpop[i,1: (m-2)] = encpar
    genpop[i,"fitval"] = sol$par
  }
  return(genpop)
}

# Kod 3.6: Tur toplamı uyum fonksiyonu 1
# Bağımlılık: Kod 3.7
tspfit = function(route, ...){
  distance = routetotal(route)
  return(distance)
}

# Kod 3.7: Toplam rota uzunluğu hesaplama
routetotal = function(route, ...) {
  route = c(route, route[1])
  embroute = embed(route, 2)[,2:1]
  return(sum(distmat[embroute]))
}

# Kod 3.8: Tur toplamı uyum fonksiyonu 2 
# Bağımlılık: Kod 3.8
tspfit2 = function(route){ 
  route = c(route, route[1])
  route = embed(route, 2)
  route = cbind(route[,2], route[,1])
  return(sum(distmat[route])^-1)
}

# Kod 3.9: Rota grafiği
plotroute = function(route, distmat, ...){
  mds = cmdscale(distmat, add=TRUE)$points
  x = mds[,1]; y = mds[,2]
  i = c(route,route[1])
  tsp = mds[i,]
  opar=par(mar=c(1,1,1,1))
  plot(x, y, asp=4/3, axes=FALSE, xlab="", ylab="", 
    main = "TSP Rota Grafiği")
  points (x, y, pch=19, cex=5, col="orange")
  arrows(tsp[i,1], tsp[i,2], tsp[i+1,1], tsp[i+1,2],
  angle=10, col="lightgray", lwd=12)
  arrows(tsp[i,1], tsp[i,2], tsp[i+1,1], tsp[i+1,2],
    angle=10, col="blue", lwd=2)
  text(x, y+20, pos=3, labels=rownames(distmat), cex = 0.8)
  par(opar)
}

# Kod 3.10: Yeni rota oluşturma 1
changeroute = function(route){
  nr = length(route)  
  route = sample(route, size=nr, replace=FALSE)
  return(route)
}

# Kod 3.11: Yeni rota oluşturma 2
changeroute1 = function(route){  
  i = seq(2, nrow(distmat)-1, 1)
  v = sample(i, size=2, replace=FALSE)
  tmp = route[v[1]]
  route[v[1]] = route[v[2]]
  route[v[2]] = tmp
  return(route)
}

# Kod 3.12: GA+optim melezleme fonksiyonu
hgaoptim2 = function(genpop, fitfunc, hgaparams,
  hgaftype, hgans, ...){
  n = nrow(genpop)
  m = ncol(genpop)
  if(missing(hgans)) hgans = 1
  if(hgans > n) hgans = n
  if(missing(hgaftype)) hgaftype="w" 
  if(hgaftype=="w"){ 
    selpop = genpop[order(genpop[,m], decreasing=TRUE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="w"){
    selpop = genpop[order(genpop[,m], decreasing=FALSE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="r"){
    selpop = genpop[sample(1:n, size=hgans, replace=FALSE),]
  }
  newpop = matrix(NA, nrow=hgans, ncol=m)
  rnselpop = rownames(selpop)
  route = seq(1, m-2)
  for(i in rnselpop){
    sol = optim( 
      par = genpop[i,1: (m-2)],
      fn = routetotal,
      gr = changeroute,
      method = hgaparams$method,
      control = hgaparams$control)
    sol = c(unname(unlist(sol)))
    genpop[i,1: (m-2)] = sol[1: (m-2)]
    genpop[i,m] = sol[m-2+1]
  }
  return(genpop)
}

# Kod 3.13: Rosenbrock uyum fonksiyonu
rosenbrock = function(x){
 m = length(x)
 xi = x[1: (m-1)]
 xj = x[2:m]
 fitval = sum(100*(xi-xj^2)^2+(xj-1)^2)
 return(fitval)
}

# Kod 3.14: Başarım karşılaştırma grafiği fonksiyonu
plotcompare = function(permat, idx=c(3,5)){
  if(is.vector(permat))
    permat = matrix(permat, nrow=1, ncol=8)
  if(!is.vector(idx)) idx = c(idx)
  if(length(idx)==1) 
    opar = par(mfrow=c(1,1))
  if(length(idx)==2) 
    opar = par(mfrow=c(2,1))
  if(length(idx)>2) 
    opar = par(mfrow=c(2,2))
  for(i in idx){
    gtitle = colnames(permat)[i]
    barplot(permat[,i], col=gray.colors(12),
    xlab="Comparison", ylab=gtitle,
    main=paste(gtitle))
  }
par(opar)
}

# Kod 3.15: Rosenbrock test fonksiyonu
rosenbrock = function(x, y, ...){
  if(missing(x)) stop("x eksik")
  if(missing(y)){
    fitval = 0
    for(i in 1: (length(x)-1))
       fitval = fitval + 100*(x[i+1]-x[i]^2)^2+(1-x[i])^2  
  }else{
    fitval = c()
    for(i in 1:length(x)){
      fit = 0
        fit = fit + 100*(y[i]-x[i]^2)^2+(1-x[i])^2 
      fitval[i]=fit 
    }
  }
  return(fitval)
}

# Kod 3.16: Özel renk paleti tanımlama
ml.colors = function(k){
  colvec = c("#00007F", "#0000FF", "#007FFF", "#00FFFF",
  "#7FFF7F", "#FFFF00", "#FF7F00", "#FF0000", "#7F0000")
  palette = grDevices::colorRampPalette(colvec, space="Lab")
  palette(k)
}

# Kod 3.17: Rastrigin fonksiyonu
rastrigin = function(x, y, A){
  if(missing(A)) A=20
  if(missing(x)) stop("x eksik")
  if(missing(y)){
    v = length(x)
    fitval = v*A + sum(x[1:v]^2-A*cos(2*pi*x[1:v]))
  }else{
    fitval = 2*A + x^2-A*cos(2*pi*x) + y^2-A*cos(2*pi*y)
  }
  return(fitval)
}

# Kod 3.18: Himmelblau fonksiyonu 
himmelblau = function(x, y,... ){
  if(missing(x)) stop("x eksik")
  if(missing(y)){
    fitval = (x[1]^2+x[2]-11)^2 +(x[1]+x[2]^2-7)^2
  }else{
    fitval = (x^2+y-11)^2 +(x+y^2-7)^2
  }
  return(fitval)
}

# Kod 3.19: makeSingleObjectiveFunction ile yeni test fonksiyonu 
  if(!require(smoof)){ 
    install.packages("smoof", 
      repo="https://cloud.r-project.org"); require(smoof)}
fn = smoof::makeSingleObjectiveFunction(
   name = "Test",
   fn = function(x) (1/sin(x[1]^2)+cos(x[2])),
     par.set = makeNumericParamSet(
     len=2L, id="t", lower=c(-1.5, -1.5), upper=c(1.5, 1.5)
   )
)

# Kod 3.20: snof ile yeni test fonksiyonu 
  if(!require(smoof)){ 
    install.packages("smoof", 
      repo="https://cloud.r-project.org"); require(smoof)}
fn = snof(
     name = "Test",
     fn = function(x) (1/sin(x[1]^2)+cos(x[2])),
     par.len = 2L, par.id = "t", 
     par.lower = -1.5, par.upper = 1.5
)

# Kod 3.21: Farklı türden parametreleri olan fonksiyon tanımlama
  if(!require(smoof)){ 
    install.packages("smoof", 
      repo="https://cloud.r-project.org"); require(smoof)}
fnkarisik = smoof::makeSingleObjectiveFunction(
  name = "4-parametreli tek amaçlı fonksiyon",
  fn = function(x) {
    if (x$disc1 == "a") {
     (x$x1^2 + x$x2^2) + 10 * as.numeric(x$logic)
    }else{
     x$x1 + x$x2-10 * as.numeric(x$logic)
    }
  },
  has.simple.signature = FALSE,
  par.set = ParamHelpers::makeParamSet(
    makeNumericParam("x1", lower = -5, upper = 5),
    makeNumericParam("x2", lower = -3, upper = 3),
    makeDiscreteParam("disc1", values = c("a", "b")),
    makeLogicalParam("logic")
  )
)

# Kod 3.22: Test fonksiyonları eşyükselti eğrileri monitör fonksiyonu
montestfunction = function(g, x,... ){
  if(!require(smoof)){ 
    install.packages("smoof", 
      repo="https://cloud.r-project.org"); require(smoof)}
  if(!require(ggplot2)){
    install.packages("ggplot2", 
      repo="https://cloud.r-project.org"); require(ggplot2)}
  plot(fn, col=4, lwd=2)
  points(x[1], x[2], pch=19, col=2, cex=2)
  text(x[1], x[2], 
    label=paste("UD:",round(fn(x),3)), cex=0.8, pos=3)
  text(x[1], x[2], label=paste("Kuşak:", g), cex=0.8, pos=1)
}

# Kod 3.23: GA+optimx melezleme fonksiyonu
hgaoptimx = function(genpop, fitfunc, hgaparams,
  hgaftype, hgans, ...){
  n = nrow(genpop)
  m = ncol(genpop)
  if(missing(hgans)) hgans = 1
  if(hgans > n) hgans = n
  if(missing(hgaftype)) hgaftype="w" 
  if(hgaftype=="b"){ 
    selpop = genpop[order(genpop[,m], decreasing=TRUE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="w"){
    selpop = genpop[order(genpop[,m], decreasing=FALSE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="r"){
    selpop = genpop[sample(1:n, size=hgans, replace=FALSE),]
  }
  newpop = matrix(NA, nrow=hgans, ncol=m)
  rnselpop = rownames(selpop)
  for(i in rnselpop){
    sol = optimx::optimx( 
      fn = fitfunc,
      par = genpop[i,1: (m-2)],
      method = hgaparams$method,
      lower = hgaparams$lower,
      upper = hgaparams$upper,
      control = hgaparams$control)
    sol = c(unname(unlist(sol)))
    genpop[i,1: (m-2)] = sol[1 : (m-2)]
    genpop[i,m] = sol[m-2+1]
  }
  return(genpop)
}

# Kod 3.24: GA'yı üstsezgisel algoritmalarla melezleme fonksiyonu
hgamho = function(genpop, fitfunc, hgaparams,
  hgaftype, hgans, ... ){
  if(!require(metaheuristicOpt)){
    install.packages("metaheuristicOpt", 
    repo="https://cloud.r-project.org"); library(metaheuristicOpt)
  }
  n = nrow(genpop)
  m = ncol(genpop)
  if(missing(hgans)) hgans = 1
  if(hgans > n) hgans = n
  if(missing(hgaftype)) hgaftype="w" 
  if(hgaftype=="b"){ 
    selpop = genpop[order(genpop[,m], decreasing=TRUE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="w"){
    selpop = genpop[order(genpop[,m], decreasing=FALSE),]
    selpop = selpop[1:hgans,]
  }else if(hgaftype=="r"){
    selpop = genpop[sample(1:n, size=hgans, replace=FALSE),]
  }
  rnselpop = rownames(selpop)
  nv = length(hgaparams$lower)
  bounds = matrix(c(hgaparams$lower[1],hgaparams$upper[1]), nr=2)
   for(i in rnselpop){
# SI optimizasyonu yap
  resultMHO = metaOpt(FUN=hgaparams$fitfunc,
    algorithm=hgaparams$algo,
    optimType=hgaparams$objective, 
    numVar=nv, rangeVar=bounds,
    control=hgaparams$control)
    genpop[i,1: (m-2)] = resultMHO$result
    genpop[i,m] = resultMHO$optimumValue
  }
  return(genpop)
}

# Kod 4.1: Reklam sayısı optimizasyonu uyum fonksiyonu 1
reklamfit1 = function(x, penr, penc,  penbeta, pengamma,
  objective,...){
  if(missing(penr)) penr=sqrt(.Machine$double.xmax)
  if(missing(penc)) penc=sqrt(.Machine$double.xmax)
  if(missing(penbeta)) penbeta=2
  if(missing(pengamma)) pengamma=1
  if(missing(objective)) objective="max"
  s = ifelse(objective=="min", 1, -1)
  fobj = function(x) (1200000 * x[1] + 2000000 * x[2])
  g1 = function(x) (-1100*x[1]-1600*x[2]+50000)
  g2 = function(x) (x[1]+x[2]-15)
  g3 = function(x) (x[1]-3)
  g4 = function(x) (x[2]-3)
  G=c()
  G[1] = max(0, g1(x))^penbeta
  G[2] = max(0, g2(x))^penbeta
  G[3] = max(0, g3(x))^penbeta
  fitval = fobj(x) + s*sum(penr*G)
  return(fitval)
}

# Kod 4.2: Reklam sayısı optimizasyonu uyum fonksiyonu 2
reklamfit2 = function(x, ...){
  fobj = function(x) (1200000*x[1] + 2000000*x[2])
  g1 = function(x) (1100*x[1]+1600*x[2])
  g2 = function(x) (x[1]+x[2])
  g3 = function(x) (x[1])
  g4 = function(x) (x[2])
  fitval = fobj(x)
  if(g1(x)>50000) fitval = -1*sqrt(.Machine$double.xmax)
  if(g2(x)<15) fitval = -1*sqrt(.Machine$double.xmax)
  if(g3(x)<3) fitval = -1*sqrt(.Machine$double.xmax)
  if(g4(x)<3) fitval = -1*sqrt(.Machine$double.xmax)
  return(fitval)
}

# Kod 4.3: Reklam sayısı optimizasyonu uyum fonksiyonu 3
reklamfit3 = function(x) (-1*reklamfit1(x))

# Kod 4.4: Reklam sayısı optimizasyonu uyum fonksiyonu 3
reklamfit4 = function(x) (-1*reklamfit2(x))

# Kod 5.1: Model uyum iyiliği ölçütleri
modelmetrics = function(yactual, ypredict){
  res = yactual-ypredict
  rss = sum(res^2)
  tss = sum((yactual-mean(yactual))^2)
  r2 = 1-rss/tss
  mae = mean(abs(res))
  mse = mean(res^2)
  rmse = sqrt(mse)
  metrics = data.frame(
   MAE = mae, MSE=mse, RMSE=rmse, R2=r2
  )
  return(metrics)
}

# Kod 5.2: Ridge regresyon grafiğine değişken adı ekleme
addvarnames = function(model, osx=1){
  lamda = length(model$lambda)
  c1 = log(model$lambda[lamda]) + osx
  c2 = model$beta[,lamda]
  varnames = names(c2)
  text(c1, c2, labels=varnames)
}

# Kod 5.3: Doğrusal regresyon modelleri için uyum fonksiyonu
slumpfit = function(x){
  predictors = xnames[x==1]
  nselected = length(predictors)
  response = yname
  if(is.null(criterion)) criterion = 6
  if(criterion < 1 | criterion > 6) criterion = 6
  critval = numeric(6)
  if(nselected > 0)
    model = lm(slumptrain[c(response, predictors)])
  else 
    model = lm(as.formula(paste(response,"~1")),
     data=slumptrain)
  critval[1] = summary(model)$r.squared # R-kare
  critval[2] = summary(model)$adj.r.squared #Düz. R-kare
  critval[3] = sum(resid(model)^2)  #RSS
  critval[4] = sqrt(deviance(model)/df.residual(model)) #RSE
  critval[5] = critval[3]/length(model$residuals) #MSE
  critval[6] = sqrt(mean(model$residuals^2)) #RMSE
  return(critval[criterion])
}

# Kod 5.4: Lojistik regresyon için uyum fonksiyonu
pimafit = function(x){
  if(!require(ModelMetrics)){
    install.packages("ModelMetrics",
    repo="https://cloud.r-project.org");
    library(ModelMetrics)
  }
  predictors = xnames[x==1]
  response = yname
  nselected = length(predictors)
  if(is.null(criterion)) criterion = 1
  if(criterion < 1 | criterion > 7) criterion = 1
  critval = c()
  if(nselected > 0){
    model = glm(pimatrain[c(response, predictors)], 
     family=binomial(link='logit'))
  }else{
     model = glm(y~1, data=pimatrain, 
      family=binomial(link='logit')) 
  } 
  critval[1] = ModelMetrics::auc(model)
  critval[2] = ModelMetrics::ce(model)
  critval[3] = ModelMetrics::mse(model)
  critval[4] = ModelMetrics::rmse(model)
  critval[5]= sapply(list(model), AIC)
  critval[6]= sapply(list(model), BIC)
  critval[7]= ModelMetrics::auc(model)/nselected
  return(critval[criterion])
}

# Kod 5.5: Yem formülasyon uyum fonksiyonu (Ölüm cezalı)
cfpfit1 = function(x,...){
  fobj = function(x) (t(feeds[,4]) %*% x)
  g1 = function(x) (t(feeds[,1]) %*% x)
  g2 = function(x) (t(feeds[,2]) %*% x)
  g3 = function(x) (t(feeds[,3]) %*% x)
  g4 = function(x) (x[2])
  g5 = function(x) (x[4])
  g6 = function(x) (sum(x))
  fitval = fobj(x)
  if (g1(x)<30) fitval = sqrt(.Machine$double.xmax)
  if (g2(x)<250) fitval = sqrt(.Machine$double.xmax)
  if (g3(x)<0.5 | g3(x)>1.5) fitval = sqrt(.Machine$double.xmax)
  if (g4(x)<8) fitval = sqrt(.Machine$double.xmax)
  if (g5(x)<20) fitval = sqrt(.Machine$double.xmax)
  if (g6(x)<100 | g6(x)>490) fitval= sqrt(.Machine$double.xmax)
  return(fitval)
}

# Kod 5.6: Yem formülasyon uyum fonksiyonu (Dış ceza)
cfpfit2 = function(x,...){
  r = sqrt(.Machine$double.xmax)
  c = sqrt(.Machine$double.xmax)
  fobj = function(x) (sum(feeds[,4]*x))
  g1 = function(x) (t(-feeds[,1]) %*% x + 30)
  g2 = function(x) (t(-feeds[,2]) %*% x + 250)
  g3 = function(x) (t(feeds[,3]) %*% x - 1.5)
  g4 = function(x) (t(-feeds[,3]) %*% x + 0.5)
  g5 = function(x) (-x[2]+8)
  g6 = function(x) (-x[4]+20)
  h1 = function(x) (sum(-x)+100)
  G=c(); H=c()
  G[1] = max(0, g1(x))^2
  G[2] = max(0, g2(x))^2
  G[3] = max(0, g3(x))^2
  G[4] = max(0, g4(x))^2
  G[5] = max(0, g5(x))^2
  G[6] = max(0, g6(x))^2
  H[1] = max(0, h1(x))^1
  fitval = fobj(x) + 1*(r*sum(G) + c*sum(H))
  return(fitval)
}

# Kod 5.7: Uyum fonksiyonu
# Bağımlılık: Kod 2.9, 2.10
sspfit = function(x,...){
  r = sqrt(.Machine$integer.max) ; s=1
  k = length(M)
  xint = decode4int(x, M=M)
  fobj = function(x) (sum(x))
  g1 = function (x) sum(-x[-7])+3
  g2 = function (x) sum(-x[-1])+5
  g3 = function (x) sum(-x[-2])+10
  g4 = function (x) sum(-x[-3])+12
  g5 = function (x) sum(-x[-4])+8
  g6 = function (x) sum(-x[-5])+3
  g7 = function (x) sum(-x[-6])+2
  G = c()
  G[1] = max(0, g1(xint))^2
  G[2] = max(0, g2(xint))^2
  G[3] = max(0, g3(xint))^2
  G[4] = max(0, g4(xint))^2
  G[5] = max(0, g5(xint))^2
  G[6] = max(0, g6(xint))^2
  G[7] = max(0, g7(xint))^2
  fitval = fobj(xint) + s*(r*sum(G))
  return(fitval)
}

# Kod 5.8: Görüntü RGB kanallarını veri çerçevesine okuma
getRGB = function(goruntu){
  boyut = dim(goruntu)
  df = data.frame(
    x = rep(1:boyut[2], each=boyut[1]),
    y = rep(boyut[1]:1, boyut[2]),
    R = as.vector(goruntu[,,1]),
    G = as.vector(goruntu[,,2]),
    B = as.vector(goruntu[,,3])
  )
  return(df)
}

# Kod 5.9: Uyum fonksiyonu - Ortalama silüet genişliği (ASW)
aswfit = function(x, penalty=TRUE, ...){
  if(!require(cluster)){ 
    install.packages("cluster", 
    repo="https://cloud.r-project.org")
  }
  if(!require(Rfast)){
    install.packages("Rfast", 
    repo="https://cloud.r-project.org")
  }
  fitval = NA
  xmat = matrix(x, nrow=clusterGA$k, ncol=length(x)/clusterGA$k)
  clusters = Rfast::dista(clusterGA$ds, xmat, 
    "euclidean", square=TRUE)
  clusters = apply(clusters, 1, which.min)
  if(length(unique(clusters)) < clusterGA$k){
    fitval = -1  
  }else{
    asw = cluster::silhouette(clusters, clusterGA$dist)  
    fitval = tryCatch(
      summary(asw)$avg.width, 
      error=function(e) {return(0)}
    )
    # Ceza fonksiyonu
    if(penalty){
      penalty = 0
      total = apply(xmat, 1, sum)
      ofval = which(total > 100)
      nconsts = length(ofval)
      if(nc > 0) {
        penalty = ncosts*max(abs(total[ofval]-100)/total[ofval])
      }
      fitval = fitval-penalty
    }
  }
  return(fitval)
}

# Kod 6.1a: MOOP için Binh'in test fonksiyonu 4
binhtest4 = function(x){
  f1 = function(x) (x[1]^2-x[2])
  f2 = function(x) (-0.5*x[1]-x[2]-1)
  fitvals = c(f1(x), f2(x))
  return(fitvals)
}

# Kod 6.1b: Binh test fonksiyonunun 4 kısıtları
binhtest4const = function(x){
  g1 = function(x) (x[1]/6+x[2]-6.5)
  g2 = function(x) (0.5*x[1]+x[2]-7.5)
  g3 = function(x) (5*x[1]+x[2]-30)
  const = c(g1(x), g2(x), g3(x))
  return(const)
}

# Kod 6.1c: MOOP için Binh'in test fonksiyonu 4 için cezalı uyum f.
binhtest4_2 = function(x){
  f1 = function(x) (x[1]^2-x[2])
  f2 = function(x) (-0.5*x[1]-x[2]-1)
  fitvals = c(f1(x), f2(x))
  g1 = function(x) (x[1]/6+x[2]-6.5)
  g2 = function(x) (0.5*x[1]+x[2]-7.5)
  g3 = function(x) (5*x[1]+x[2]-30)
  if(g1(x) > 0) rep(sqrt(.Machine$double.xmax), 2)
  if(g2(x) > 0) rep(sqrt(.Machine$double.xmax), 2)
  if(g3(x) > 0) rep(sqrt(.Machine$double.xmax), 2)
  return(fitvals)
}

# Kod 6.1d: MOOP için Binh'in test fonksiyonu 4 için cezalı uyum f.
binhtest4_3 = function(i){
  f1 = x[i,1]^2-x[i,2]
  f2 = -0.5*x[i,1]-x[i,2]-1
  fitvals = c(f1, f2)
  g1 = (x[i,1]/6+x[i,2]-6.5)
  g2 = (0.5*x[i,1]+x[i,2]-7.5)
  g3 = (5*x[i,1]+x[i,2]-30)
  if(g1 > 0) fitvals = rep(sqrt(.Machine$double.xmax), 2)
  if(g2 > 0) fitvals = rep(sqrt(.Machine$double.xmax), 2)
  if(g3 > 0) fitvals = rep(sqrt(.Machine$double.xmax), 2)
  return(fitvals)
}

# Kod 6.2: Çekicilik hesaplama
desrank = function(objvals, desidx){
  desvals = apply(objvals, 1, desidx)
  idx = order(desvals, decreasing=TRUE)
  desmatrix = data.frame(objvals, desvals)
  colnames(desmatrix)[ncol(desmatrix)] = "Çekicilik"
  desmatrix = desmatrix[idx, ]
  return(desmatrix)
}

# ALIŞTIRMA KODLARI LDA ANALİZİ.......................................................................
# Kod A.1: LDA uyum fonksiyonu 1 (yanlış sınıflama oranı, FDR)
pimaldafit1 = function(x){
  if(!require(MASS)){
     install.packages("MASS", 
     repo="https://cloud.r-project.org");
     library(MASS)
  }
  maxz = function(z) which(z == max(z))
  fitval = 1
  if(sum(x) > 1){
    lda.model = MASS::lda(idvmat[,x==1], dvmat, CV=TRUE)
    posteriors = lda.model$posterior
    fitval = sum(dvmat != dimnames(posteriors)[[2]]
    [apply(posteriors, 1, maxz)])/nrow(dvmat)
  }
  return(fitval)
}

# Kod A.2: LDA tabanlı uyum fonksiyonu 2 (Doğruluk, ACC)
# Bağımlılık: Kod 3.1, 5.5; Örnek 5.25
pimaldafit2 = function(x){
  if(!require(MASS)){
     install.packages("MASS", 
     repo="https://cloud.r-project.org");
     library(MASS)
  }
  fitval = 0
  if(sum(x) > 1){
    lda.model = MASS::lda(idvmat[,x==1], dvmat, CV=TRUE)
    fitval = mean(lda.model$class==dvmat)
  }
  return(fitval)
}

# Kod A.3: QDA uyum fonksiyonu 1 (yanlış sınıflama oranı, FDR)
pimaqdafit1 = function(x){
  if(!require(MASS)){
     install.packages("MASS", 
     repo="https://cloud.r-project.org");
     library(MASS)
  }
  maxz = function(z) which(z == max(z))
  fitval = 1
  if(sum(x) > 1){
    qda.model = MASS::qda(idvmat[,x==1], dvmat, CV=TRUE)
    posteriors = qda.model$posterior
    fitval = sum(dvmat != dimnames(posteriors)[[2]]
    [apply(posteriors, 1, maxz)])/nrow(dvmat)
  }
  return(fitval)
}

# Kod A.4: QDA uyum fonksiyonu 2 (Doğruluk, ACC)
pimaqdafit2 = function(x){
  if(!require(MASS)){
     install.packages("MASS", 
     repo="https://cloud.r-project.org");
     library(MASS)
  }
  fitval = 0
  if(sum(x) > 1){
    qda.model = MASS::qda(idvmat[,x==1], dvmat, CV=TRUE)
    fitval = mean(qda.model$class==dvmat)
  }
  return(fitval)
}

# Kod A.5: MDA uyum fonksiyonu 1 (Doğruluk maksimizasyonu)
pimamdafit1 = function(x){
  if(!require(mda)){
     install.packages("mda", 
     repo="https://cloud.r-project.org");
     library(mda)
  }
  fitval = 0
  if(sum(x)>1){
    mda.model = mda::mda(dvmat~idvmat[,x==1])
    preddvar = predict(mda.model, idvmat[,x==1])
    fitval = mean(preddvar == dvmat)
  }
  return(fitval)
}

# Kod A.6: MDA uyum fonksiyonu 2 (Sapma minimizasyonu)
pimamdafit2 = function(x){
  fitval = 0
  if(!require(mda)){
     install.packages("mda", repo="https://cloud.r-project.org");
     library(mda)
  }
  if(sum(x)>1){
    mda.model = mda::mda(dvmat~idvmat[,x==1])
    fitval = mda.model$deviance
  }
  return(fitval)
}

# Kod A.7: FDA uyum fonksiyonu (Doğruluk maksimizasyonu)
pimafdafit = function(x){
  if(!require(mda)){
     install.packages("mda", 
      repo="https://cloud.r-project.org");
     library(mda)
  }
  fitval = 0
  if(sum(x)>1){
    fda.model = mda::fda(dvmat~idvmat[,x==1])
    preddvar = predict(fda.model, idvmat[,x==1])
    fitval = mean(preddvar == dvmat)
  }
  return(fitval)
}

# Kod A.8: RDA uyum fonksiyonu 1 (Doğruluk değeri maksimizasyonu)
pimardafit = function(x){
  if(!require(klaR)){
     install.packages("klaR", 
     repo="https://cloud.r-project.org");
     library(klaR)
  }
  fitval = 1
  if(sum(x)>1){
    rda.model = klaR::rda(dvmat~idvmat[,x==1])
    fitval = rda.model$error.rate[2]
  }
  return(fitval)
}
