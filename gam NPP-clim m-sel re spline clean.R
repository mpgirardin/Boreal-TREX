
vif.lme <- function (fit) {
  #https://stackoverflow.com/questions/26633483/collinearity-after-accounting-for-random-mixed-effects
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }

m1.collinearity <- function (max.vif){
  
  xVar <- names(xClim)
  colinearity.m1<- data.frame() 
  arret <- 1
  while(arret > 0){
    
    temp <- paste(xVar, collapse=" + ")
    form <- formula(paste("NPP ~ ", temp))
    model.lme <- lme(form  , random=~1|Bloc, data = y.clim)
    vif_lme<-data.frame(vif = vif.lme(model.lme))
    
    setorder(vif_lme, -vif)
    vif.max <- vif_lme[1,"vif"]
    if (vif.max > max.vif) {
      var.remove <- row.names(vif_lme)[1]
      xVar <- setdiff(xVar, var.remove)
      tmp1<-data.frame(var = var.remove, vif = vif.max)
      colinearity.m1 <- rbind(colinearity.m1, tmp1)
    }
    else arret <- 0
  }
  
  if (dim(colinearity.m1)[2] > 0) colinearity.m1$rm <- 1
  vif_lme$var<- row.names(vif_lme)
  vif_lme$rm<- 0
  colinearity.m1<- rbind(colinearity.m1,vif_lme)
  return(colinearity.m1)
}

m3.collinearity <- function (max.vif){
  xVar <- names(xClim)
  colinearity.m3<- data.frame() 
  arret <- 1
  while(arret > 0){
    
    temp <- paste(xVar, collapse=" + ")
    
    form <- formula(paste("tV ~ ", temp))
    
    model.lm <- lm(form  , data = y.clim)
    
    vif_lm<-as.data.frame(ols_vif_tol(model.lm))
    
    setorder(vif_lm, -VIF)
    vif.max <- vif_lm[1,"VIF"]
    if (vif.max > max.vif) {
      
      var.remove <- (vif_lm[1,"Variables"])
      xVar <- setdiff(xVar, var.remove)
      tmp1<-data.frame(var = var.remove, VIF = vif.max)
      colinearity.m3 <- rbind(colinearity.m3, tmp1)
    }
    else arret <- 0
  }
  
  if (dim(colinearity.m3)[2] > 0) colinearity.m3$rm <- 1
  vif_lm$Tolerance<- NULL
  vif_lm$rm<- 0
  names(vif_lm)[1]<-"var"
  colinearity.m3<- rbind(colinearity.m3,vif_lm)
  return(colinearity.m3)
  
}


gam.model.sel.pValue.Grp <-function (y, N.grp = 1, max.Var = 6, R2.tol = 0.01, bloc,resp.label, sel.by ){
  
  
  dd<-y
  # forT <- paste0("Y.",names(y))
  sf_sel_df <- data.frame()
  # print(forT)
  
  # Obj_old <- -Inf; 
  Obj <- -99999
  
  # obj	<- Obj
  nX.sel <- 0
  all.s <- data.frame()
  # selected.s <- data.table()
  
  for (iG in (1:N.grp)){
    i	<- 0
    Obj_old <- -Inf
    # print(iG)
    X<- as.data.frame( dt.use[,eval(parse(text= paste0("xVar.g", iG))), with = FALSE])
    X.iG <- X
    nX <- dim(X)[2]
    ii.range <- 1:nX
    # sf_list	<- 1:nX
    R2.sel <- data.frame()
    while( ( Obj - Obj_old >=R2.tol ) &(nX.sel < max.Var)){
      
      
      Obj_old	<- Obj
      # Obj_list<- NULL
      #print (paste0("i = ", i))
      
      for( ii in ii.range ){
        if (ii > 0){
          #print (paste0("ii = ", ii))
          isf <-as.data.frame(X[,ii,drop = FALSE])
          # names(isf) <- names(X)[ii]
          #print (paste0(ii, "   ", names(isf)))
          
          d	<- data.table( dd, isf )
          names(d)[1] <- "NPP"
          temp <- paste(names(d[, - 1]), collapse=",bs='cr')+s(")
          #print(ii)
          #print(names(d)[1])
          
          form <- formula(paste("NPP ~  s(",temp,",bs='cr')"))
        }
        print(form)
        if (bloc == 0){ 
          d <-data.table(d,dt.pb[,2])
          names(d)[1] <- "NPP"
          
          
          result <- tryCatch(
            {
              form.re <-formula(paste("NPP ~  s(",temp,",bs='cr')", "+ s(Bloc, bs = 're')"))
              #sfmod0	<- gamm(form, method = 'REML', random = list(Bloc=~1),data = d)$gam 
              sfmod0	<- gam(form.re, method = 'REML',data = d)
            },
            error = function(e) e
          )
          
        }else{
          result <- tryCatch(
            {
              sfmod0	<- gam(form, method = 'REML', data = d)},
            error = function(e) e
          )
          
        }
        
        if (is.null(result$message)) {
          Obja	<-  summary(sfmod0)$r.sq
          s.table <- as.data.frame(summary(sfmod0)$s.table)
          names(s.table)[endsWith("p-value",names(s.table) )]<-"pValue"
          s.table$rn <- row.names(s.table)
          setDT(s.table)[, var:=stringr::str_sub(rn, 3,-2)]
          s.table<-s.table[var %in% eval(parse(text= paste0("xVar.g", iG)))]
          s.table[,R2 := Obja]
          s.table[,aic :=AIC(sfmod0)]
          
          s.table[, c("bloc", "iG", "i0", "ii") := list(bloc, iG, i, ii)]
          s.table[,selected := 0]
        } 
        # else  Obj	<- -99999
        
        all.s <- rbind(all.s, s.table)
      }
      
      # all.s.i0 <- setorder(all.s[i0 == i], -R2)
      
      s.i0 <- all.s[i0 == i,]
      i0.05 <- s.i0[ pValue <= 0.05,]
      #remove the existing
      
      if (dim(i0.05)[1] > 0) {
        tmp.i0 <- i0.05[, .N, by = .(iG, i0, ii, R2, aic)][N == nX.sel+1]
        
        if (dim(tmp.i0)[1] == 0){
          tmp <- i0.05[!(var %in% names(sf_sel_df))]
          if (dim(tmp)[1] > 0) {
            if (sel.by == "R2"){
              setorder(tmp, R2 )
              var.sel <- tmp[, .SD[.N]]}
            if (sel.by == "AIC"){
              setorder(tmp, aic )
              var.sel <- tmp[, .SD[1]]}
            var.sel[,selected.i0 := 1][,selected:= NULL]
            sf_sel_df <- as.data.frame(X.iG[,var.sel$var, drop = FALSE])
            
            Obj_old <- -Inf
            all.s<- var.sel[, c("i0", "ii","var", "selected.i0")][all.s, on = c("i0", "ii","var"), nomatch = NA]
            all.s[selected.i0 == 1, selected := selected.i0]
            all.s[,selected.i0  := NULL]
            X<-X.iG[, names(sf_sel_df), drop = FALSE]
            dd	<- as.data.frame( y, drop = FALSE )
            d <-data.frame(dd, X)
            nX.sel <- dim(sf_sel_df)[2]-1
            var.sel[,nX := nX.sel+1]
            var.sel[,re.sel := 1]
            var.sel <- all.s[, -"selected" ][var.sel[, c("i0", "ii","selected.i0","nX",  "re.sel")], on =c("i0", "ii") ]
            R2.sel <- rbind(R2.sel, var.sel) 
            ii.range <- 0:0
            temp <- paste(names(sf_sel_df), collapse=",bs='cr')+s(")
            
            #form <- formula(paste("y ~s(",xChr, ",bs='cr') + s(",temp,",bs='cr')"))
            form <- formula(paste("NPP ~  s(",temp,",bs='cr')"))
            # dd <-data.frame(y, sf_sel_df)
            #print(names(dd))
            names(d)[1] <- "NPP"
            if (nX.sel == 0) Obj <- -99999
          } else Obj_old <- Obj
          
          
        }else{
          
          if (sel.by == "R2"){
            setorder(tmp.i0, R2)
            var.sel <- tmp.i0[, .SD[.N]]
            
          }
          
          if (sel.by == "AIC"){
            setorder(tmp.i0, aic)
            var.sel <- tmp.i0[, .SD[1]]
            
          }
          
          Obj <- var.sel$R2
          # chk R2
          if (Obj - Obj_old >=R2.tol){
            var.sel <- var.sel[,-c('R2','aic', 'N')]
            var.sel[,selected.i0 := 1]
            
            
            
            var.sel <- i0.05[,-'selected'][var.sel, on = .(iG, i0, ii)]
            
            sf_sel_df <- as.data.frame(X.iG[,var.sel$var, drop = FALSE])
            nX.sel <- dim(sf_sel_df)[2] 
            if (nX.sel == 1) Obj_old <- -Inf else Obj_old <-  Obj
            var.sel[,nX := nX.sel]
            var.sel[,re.sel := 0]
            setDF(var.sel)
            R2.sel <- rbind(R2.sel, var.sel)
            all.s<- setDT(var.sel)[, c("i0", "ii", "var", "selected.i0")][all.s, on = c("i0", "ii", "var"), nomatch = NA]
            all.s[selected.i0 == 1, selected := selected.i0]
            all.s[,selected.i0  := NULL]
            X<-X.iG[,!(names(X.iG) %in% names(sf_sel_df))]
            dd	<- as.data.frame( cbind(y, sf_sel_df) )
            
            
            ii.range <- 1:dim(X)[2]
          }
          
        }
      }  
      # 
      names(X)
      
      
      
      
      
      
      i		<- i + 1
      # print(paste0("i = ", i))
      # print(paste0("nX.sel = ", nX.sel))
      # print(paste0("Obj = ", Obj))
      # print(paste0("Obj_old = ", Obj_old))
      # print(names(X))
    }#end of while (selection)
    
  }#end of iG 
  #final model
  # sel <- all.s[selected == 1,]
  if (dim(sf_sel_df)[2] > 0) { 
    #temp <- paste(names(subset(dd, select=-c(1,2))), collapse=",bs='cr')+s(")
    temp <- paste(names(sf_sel_df), collapse=",bs='cr')+s(")
    
    #form <- formula(paste("y ~s(",xChr, ",bs='cr') + s(",temp,",bs='cr')"))
    form <- formula(paste("NPP ~  s(",temp,",bs='cr')"))
    dd <-data.frame(y, sf_sel_df)
    #print(names(dd))
    names(dd)[1] <- "NPP"
    if (bloc == 0){
      # print("bloc  0")
      
      dd <-data.table(dd,dt.pb[,2, drop = FALSE])
      
      
      form.re <-formula(paste("NPP ~  s(",temp,",bs='cr')", "+ s(Bloc, bs = 're')"))
      #sfmod0	<- gamm(form, method = 'REML', random = list(Bloc=~1),data = d)$gam 
      sfmod	<- gam(form.re, method = 'REML',data = dd)
      
    }else{ 
      # print("bloc non 0")
      #print(form)
      sfmod	<- gam(form, method = 'REML', data = dd)}
    
  }else {
    dd <- y
    names(dd)[1] <- "NPP"
    if (bloc == 0){
      print("no variable selected bloc  0")
      
      dd <-data.table(dd,dt.pb[,2, drop = FALSE])
      sfmod	<- gam(NPP~ s(Bloc, bs = 're'), method = 'REML', data = dd)
    } else
      
      #no selection finish the process
      sfmod	<- gam(NPP ~ 1, method = 'REML', data = dd)
  }
  # par(mfrow=c(2,3))
  # print(form)
  if (dim(sf_sel_df)[1] > 0 ){
    P.vars <- names(sf_sel_df)
    
    map(P.vars, function(x){
      p <- plotGAM(sfmod, smooth.cov = x, rawOrFitted = "fitted") +
        geom_point(data = dd, aes_string(y = "NPP", x = x), color="#0073C2FF", alpha = 0.2) +
        #  geom_rug(data = train.college, aes_string(y = "Outstate", x = x ), alpha = 0.2)
        ylab(paste0(names(y))) +
        # ggtitle()
        ggtitle(paste0(names(y), " vs. s(", x,")"))
      
      g <- ggplotGrob(p)
    }) %>%
    {grid.arrange(grobs = (.), ncol = 2)}
    
    
    myplot <- recordPlot()
    graphics.off()
    # replayPlot(myplot)
  }else myplot <-NULL
  #  print(is.null(sfmod$lme))
  #2+ linear predictors
  
  
  pred	<- predict( sfmod )
  resid	<- residuals ( sfmod)
  
  aic <- AIC(sfmod)
  # obs = pred+resid
  sfmod.summary <-summary(sfmod)
  
  r.sq = sfmod.summary$r.sq
  R.sq = data.frame(r.sq,aic)
  names(R.sq) <- c('adjRsq', 'aic')
  
  p_table <- data.frame(sfmod.summary$p.table)
  p_table <- within(p_table, {lci <- Estimate - qnorm(0.975) * Std..Error
  uci <- Estimate + qnorm(0.975) * Std..Error})
  
  p.r <-data.frame(dd[,1],pred,resid)
  
  names(p.r)[1:3] <-c(resp.label,'pred','resid')
  p.r <-data.frame(dt.pb, p.r)
  
  #  p_table
  
  # Nsf <- dim(sf_sel_df)[2] 
  if (dim(sf_sel_df)[2] > 0)
    p.r <-data.frame(p.r, sf_sel_df)
  # p.r$var <- row.names(p.r)
  names.na <- setdiff(xchr.complete, as.list(names(sf_sel_df)) )
  p.r.na <-as.data.frame(matrix(NA_real_, nrow = dim(p.r)[1], ncol = length(names.na)))
  names(p.r.na) <- names.na
  p.r <-data.frame(p.r, p.r.na)
  #setdiff(xchr.complete,c("Geno2_p" , "Geno3_p" ))
  s_table <-data.frame(sfmod.summary$s.table)
  # setDF(all.s)
  # setDF(R2.sel)
  if (dim(R2.sel)[1]  > 0){
    print("BBB")
    print(names(y))
    
    R2.sel$resp <- names(y)}
  if (dim(s_table)[1]  > 0){
    s_table$resp <- names(y)}
  
  
  p_table$resp <- names(y)
  all.s$resp <- names(y)
  p.r$resp <- names(y)
  
  #  s_table
  # hist(y)
  list(p.table = p_table, R.sq = R2.sel, s.table = s_table, Pred = p.r, sequential.sel = all.s, plots = myplot)
}

fun_var.sel <- function(y,effect.i, effect.s, max.Var = 2, R2.tol = 0.01, bloc, iG =1, sel.by = "AIC" ){
  
  dd<-y
  print(names(y))
  # forT <- paste0("Y.",names(y))
  sf_sel_df <- data.frame()
  # print(forT)
  
  # Obj_old <- -Inf; 
  Obj <- -99999
  
  # obj	<- Obj
  nX.sel <- 0
  all.s <- data.frame()
  
  i	<- 0
  Obj_old <- -Inf
  # print(iG)
   X<- as.data.frame( dt.use[,eval(parse(text= paste0("xVar.g", iG))), with = FALSE])
  # X<- as.data.frame( dt.use[,xVar, with = FALSE])
  X.iG <- X
  nX <- dim(X)[2]
  ii.range <- 1:nX
  # sf_list	<- 1:nX
  R2.sel <- data.frame()
  while( ( Obj - Obj_old >=R2.tol ) &(nX.sel < max.Var)){
    
    
    Obj_old	<- Obj
    # Obj_list<- NULL
    #print (paste0("i = ", i))
    
    for( ii in ii.range ){
      if (ii > 0){
        #print (paste0("ii = ", ii))
        isf <-as.data.frame(X[,ii,drop = FALSE])
        # names(isf) <- names(X)[ii]
        #print (paste0(ii, "   ", names(isf)))
        
        d	<- data.table( dd, isf )
        names(d)[1] <- "NPP"
        temp <- paste(names(d[, - 1]), collapse=paste0(",bs='cr'", effect.s, ")+s("))
        #print(ii)
        #print(names(d)[1])
        
        form <- formula(paste("NPP ~ ", effect.i, " s(",temp,",bs='cr'", effect.s, ")"))
      }
      # print(form)
      if (bloc == 0){ 
        d <-data.table(d,dt.pb)
        names(d)[1] <- "NPP"
        
        
        result <- tryCatch(
          {
            form.re <-formula(paste("NPP ~ ", effect.i, " s(",temp,",bs='cr'", effect.s, ")", "+ s(Bloc, bs = 're')"))
            #sfmod0	<- gamm(form, method = 'REML', random = list(Bloc=~1),data = d)$gam 
            sfmod0	<- gam(form.re, method = 'REML',data = d)
          },
          error = function(e) e
        )
        
      }else{
        result <- tryCatch(
          {
            sfmod0	<- gam(form, method = 'REML', data = d)},
          error = function(e) e
        )
        
      }
      
      if (is.null(result$message)) {
        Obja	<-  summary(sfmod0)$r.sq
        p.table <- as.data.frame(summary(sfmod0)$p.table)
        names(p.table)[endsWith("Pr(>|t|)",names(p.table) )]<-"pValue"
        p.table$rn <- row.names(p.table)
        p.table <- setDT(p.table)[, var:=stringr::str_sub(rn, 1,8)][var =="genoType" ]
        if (dim(p.table)[1] > 1)
          pV.i <-metap::sumlog(p.table$pValue)$p
        else pV.i <- NA_real_
        
        s.table <- as.data.frame(summary(sfmod0)$s.table)
        names(s.table)[endsWith("p-value",names(s.table) )]<-"pValue"
        s.table$rn <- row.names(s.table)
        setDT(s.table)[,rn.var:=tstrsplit(rn, ":", fixed=TRUE, keep=1L)][, var:=stringr::str_sub(rn.var, 3,-2)]
        s.table<-s.table[var %in% eval(parse(text= paste0("xVar.g", iG)))]
        # pV.s <-metap::sumlog(s.table$pValue)$p
        if (dim(s.table)[1] > 2){
          pV.s<-s.table[,list(pValue.slope =metap::sumlog(s.table$pValue)$p ), by = "var"]
        }else {
          pV.s<-s.table[,c("var","pValue")]
          names(pV.s)[2]<- "pValue.slope"
        }
        
        stats.tmp <- data.table(pV.s
                                , pValue.intercept = pV.i,R2 = Obja, aic=AIC(sfmod0))
        
        
        stats.tmp[, c("bloc", "iG", "i0", "ii") := list(bloc, iG, i, ii)]
        stats.tmp[,selected := 0]
      } 
      # else  Obj	<- -99999
      if(!(is.null(stats.tmp)))
      all.s <- rbind(all.s, stats.tmp)
    }
    
    # all.s.i0 <- setorder(all.s[i0 == i], -R2)
    
    s.i0 <- all.s[i0 == i,]
    s.i0[ is.na(pValue.intercept) & pValue.slope <= 0.05, pre.sel := 1]
    s.i0[ !is.na(pValue.intercept) & pValue.intercept <= 0.05 & pValue.slope <= 0.05, pre.sel := 1]
    i0.05 <- s.i0[ pre.sel == 1,]
    i0.05[ ,pre.sel := NULL]
    #remove the existing
    
    if (dim(i0.05)[1] > 0) {
      tmp.i0 <- i0.05[, .N, by = .(iG, i0, ii, R2, aic)][N == nX.sel+1]
      
      if (dim(tmp.i0)[1] == 0){
        tmp <- i0.05[!(var %in% names(sf_sel_df))]
        if (dim(tmp)[1] > 0) {
          if (sel.by == "R2"){
            setorder(tmp, R2 )
            var.sel <- tmp[, .SD[.N]]}
          if (sel.by == "AIC"){
            setorder(tmp, aic )
            var.sel <- tmp[, .SD[1]]}
          var.sel[,selected.i0 := 1][,selected:= NULL]
          sf_sel_df <- as.data.frame(X.iG[,var.sel$var, drop = FALSE])
          
          Obj_old <- -Inf
          all.s<- var.sel[, c("i0", "ii","var", "selected.i0")][all.s, on = c("i0", "ii","var"), nomatch = NA]
          all.s[selected.i0 == 1, selected := selected.i0]
          all.s[,selected.i0  := NULL]
          X<-X.iG[, names(sf_sel_df), drop = FALSE]
          dd	<- as.data.frame( y, drop = FALSE )
          d <-data.frame(dd, X)
          nX.sel <- dim(sf_sel_df)[2]-1
          var.sel[,nX := nX.sel+1]
          var.sel[,re.sel := 1]
          var.sel <- all.s[, -"selected" ][var.sel[, c("i0", "ii","selected.i0","nX",  "re.sel")], on =c("i0", "ii") ]
          R2.sel <- rbind(R2.sel, var.sel) 
          ii.range <- 0:0
          temp <- paste(names(sf_sel_df), collapse=",bs='cr')+s(")
          
          #form <- formula(paste("y ~s(",xChr, ",bs='cr') + s(",temp,",bs='cr')"))
          form <- formula(paste("NPP ~  s(",temp,",bs='cr')"))
          # dd <-data.frame(y, sf_sel_df)
          #print(names(dd))
          names(d)[1] <- "NPP"
          if (nX.sel == 0) Obj <- -99999
        } else Obj_old <- Obj
        
        
      }else{
        
        if (sel.by == "R2"){
          setorder(tmp.i0, R2)
          var.sel <- tmp.i0[, .SD[.N]]
          
        }
        
        if (sel.by == "AIC"){
          setorder(tmp.i0, aic)
          var.sel <- tmp.i0[, .SD[1]]
          
        }
        
        Obj <- var.sel$R2
        # chk R2
        if (Obj - Obj_old >=R2.tol){
          var.sel <- var.sel[,-c('R2','aic', 'N')]
          var.sel[,selected.i0 := 1]
          
          
          
          var.sel <- i0.05[,-'selected'][var.sel, on = .(iG, i0, ii)]
          
          sf_sel_df <- as.data.frame(X.iG[,var.sel$var, drop = FALSE])
          nX.sel <- dim(sf_sel_df)[2] 
          if (nX.sel == 1) Obj_old <- -Inf else Obj_old <-  Obj
          var.sel[,nX := nX.sel]
          var.sel[,re.sel := 0]
          setDF(var.sel)
          R2.sel <- rbind(R2.sel, var.sel)
          all.s<- setDT(var.sel)[, c("i0", "ii", "var", "selected.i0")][all.s, on = c("i0", "ii", "var"), nomatch = NA]
          all.s[selected.i0 == 1, selected := selected.i0]
          all.s[,selected.i0  := NULL]
          X<-X.iG[,!(names(X.iG) %in% names(sf_sel_df))]
          dd	<- as.data.frame( cbind(y, sf_sel_df) )
          
          
          ii.range <- 1:dim(X)[2]
        }
        
      }
    }  
    # 
    names(X)
    
    
    
    
    
    
    i		<- i + 1
    # print(paste0("i = ", i))
    # print(paste0("nX.sel = ", nX.sel))
    # print(paste0("Obj = ", Obj))
    # print(paste0("Obj_old = ", Obj_old))
    # print(names(X))
  }#end of while (selection)
  if (dim(R2.sel)[1]>0){ 
    setDT(R2.sel)
    setkey(R2.sel, iG, nX)
    R2.sel <- R2.sel[,.SD[.N], by = .(iG, nX)]
    R2.sel[,c("effect.i","effect.s", "resp") := list(effect.i,effect.s, names(y))]
  }
  
  if (dim(sf_sel_df)[2] > 0) { 
    #temp <- paste(names(subset(dd, select=-c(1,2))), collapse=",bs='cr')+s(")
    
    temp <- paste(names(sf_sel_df), collapse=paste0( effect.s, ")+s("))
    
    #print(ii)
    #print(names(d)[1])
    
    form <- formula(paste("NPP ~ ", effect.i, " s(",temp, effect.s, ")"))
    
    dd <-data.frame(y, sf_sel_df)
    
    dd <-data.table(dd,dt.pb)
    #print(names(dd))
    names(dd)[1] <- "NPP"
    if (bloc == 0){
      # print("bloc  0")
      
      form.re <- formula(paste("NPP ~ ", effect.i, " s(",temp, effect.s, ")", "+ s(Bloc, bs = 're')"))
      
      sfmod	<- gam(form.re, method = 'REML',data = dd)
      
    }else{ 
      # print("bloc non 0")
      #print(form)
      sfmod	<- gam(form, method = 'REML', data = dd)}
    
  }else {
    dd <- y
    names(dd)[1] <- "NPP"
    if (bloc == 0){
      print("no variable selected bloc  0")
      
      dd <-data.table(dd,dt.pb[,2, drop = FALSE])
      sfmod	<- gam(NPP~ s(Bloc, bs = 're'), method = 'REML', data = dd)
    } else
      
      #no selection finish the process
      sfmod	<- gam(NPP ~ 1, method = 'REML', data = dd)
  }
  
  return(list(R2.sel = R2.sel, all.s = all.s, dd = dd, sfmod = sfmod, sf_sel_df = sf_sel_df))
  
}
# 
#R2 in vertical
extract.gam.m13.R2v <-function ( model.s, max.Var = 4){
  #one year only
  if (length(model.s) == 5){
    R2<- model.s $R.sq
    Pred <- model.s$Pred
    s.table <- model.s$s.table
    all.s <- model.s$sequential.sel
    
  } else {
    p.table <-do.call("rbind", lapply(model.s, '[[', "p.table"))
    s.table <-do.call("rbind", lapply(model.s, '[[', "s.table"))
    R2 <- do.call("rbind", lapply(model.s, '[[', "R.sq"))
    all.s <- do.call("rbind", lapply(model.s, '[[', "sequential.sel"))
    #obs vs. pred
    # R2.inc <-do.call("rbind", lapply(model.s, '[[', "r2.inc"))
    Pred <-do.call("rbind.fill", lapply(model.s, '[[', "Pred"))
    plot.list <-lapply(model.s, '[[', "plots")
    
  }
  setDF(Pred)
  Pred<-Pred[,(which(colSums(!is.na(Pred)) > 0))]
  
  #sTable <- out$s.table
  s.table$rn <- row.names(s.table)
  # all.s$rn <- row.names(all.s)
  
  list(R2 = R2,  s.table = s.table, Pred = Pred, all.s = all.s, plot.list = plot.list)
}

outcsv <- function(model.label , pbop = "pb"){
  
  write.csv(R2, file = file.path(folder.dest, paste0(model.label,".", respL,".", pbop," R2-Vsel.csv")), row.names = FALSE)
  write.csv(Pred, file = file.path(folder.dest, paste0(model.label,".", respL,".",pbop," Pred.csv")), row.names = FALSE)
  write.csv(sTable, file = file.path(folder.dest, paste0(model.label,".", respL,".",pbop," sTable.csv")), row.names = FALSE)
  write.csv(seq.sel, file = file.path(folder.dest, paste0(model.label,".", respL,".",pbop," sequential selection.csv")), row.names = FALSE)
}



#for model2
gam.model2 <-function (y, form.s){
  y0<-data.table(y,dt.py)
  names(y0)[1] <- "V"
  dt.m2.b <-  dcast(y0, Year + Bloc ~ var, value.var="V")
  dt.m2.b <- dt.m2.b[clim.local, on = .(Year), nomatch = 0]
  setorder(dt.m2.b, Bloc, Year)  
  names(dt.m2.b)
  
  
  result <- tryCatch(
    {
      NPP.clim.b0a <- gamm(form.s, 
                           
                           random = list(Bloc=~1), 
                           correlation=corAR1(), select = TRUE,
                           data =  dt.m2.b)
      
      
      #r2.m2a.all <- rbind(r2.m2a.all, r.sq)
      
      
      #setDF(t.lme)
    },
    error = function(e) e
  )
  if (is.null(result$message)) {
    r.sq <- data.frame( r2 = summary(NPP.clim.b0a$gam)$r.sq)
    t.lme <- as.data.frame(summary(NPP.clim.b0a$lme)$tTable)
    t.lme$var <- row.names(t.lme) }else  {
      t.lme <- data.frame()
      r.sq <-NULL}
  list(r2 = r.sq, tValue = t.lme)
  
}

extract.gam.m2 <-function (model.s)
{
  
  R2.m2 <-  do.call("rbind", lapply(model.s, '[[', "r2"))
  R2.m2 $Prov<- row.names(R2.m2)
  tValue.m2 <-  do.call("rbind", lapply(model.s, '[[', "tValue"))
  names(tValue.m2)[4] <-"tvalue"
  tValue.m2 $rn<- row.names(tValue.m2)
  setDT(tValue.m2)[, c("Prov", "xVar") := tstrsplit(rn, ".X", fixed=TRUE)] 
  tValue.m2[, rn:=NULL]
  setcolorder(tValue.m2, c("Prov", "xVar"))
  
  
  tmp <-tValue.m2[!(stringr::str_detect(tValue.m2$var, "(Year)")), ]
  
  tmp <-tmp[!(stringr::str_detect(tmp$var, "(Intercept)")), ]
  
  tValue.m2.wide <- dcast(tmp,Prov  ~ xVar, value.var = "tvalue")
  # tValue.m2.wide <- transform(tValue.m2.wide, Prov = as.numeric(Prov))
  tValue.m2.wide<-tValue.m2.wide[,  c("Prov", "prevSummerTemp" , "prevFallTemp" , "WinterTemp" , "SpringTemp" , "SummerTemp" , "FallTemp",
                                      "prevSummerSMI" , "prevFallSMI" , "Snowfall" , "SpringSMI" , "SummerSMI" , "FallSMI")]
  #setcolorder(tValue.m2.wide, c("Prov", "prevSummerTemp" , "prevFallTemp" , "WinterTemp" , "SpringTemp" , "SummerTemp" , "FallTemp",
  #"prevSummerSMI" , "prevFallSMI" , "Snowfall" , "SpringSMI" , "SummerSMI" , "FallSMI"))
  #setDF(tValue.m2)
  list(r2 = R2.m2, tTable = tValue.m2, tValue = tValue.m2.wide)
  
}
  

fun_var.sel.s.default <- function(y,effect.i, effect.s, max.Var = 2, R2.tol = 0.01, bloc, iG =1, sel.by = "AIC",standalone = 0 ){
  #standalone version for variable selection
  #it can be used within by setting standalone = 0
 # if (standalone == 1) {
  dd<-y
  print(names(y))
  # }
  # forT <- paste0("Y.",names(y))
  sf_sel_df <- data.frame()
  # print(forT)
  
  # Obj_old <- -Inf; 
  Obj <- -99999
  
  # obj	<- Obj
  nX.sel <- 0
  all.s <- data.frame()
  
  i	<- 0
  Obj_old <- -Inf
  # print(iG)
  X<- as.data.frame( dt.use[,eval(parse(text= paste0("xVar.g", iG))), with = FALSE])
  # X<- as.data.frame( dt.use[,xVar, with = FALSE])
  X.iG <- X
  nX <- dim(X)[2]
  ii.range <- 1:nX
  # sf_list	<- 1:nX
  R2.sel <- data.frame()
  
  while( ( Obj - Obj_old >=R2.tol ) &(nX.sel < max.Var)){
    
    
    Obj_old	<- Obj
    # Obj_list<- NULL
    #print (paste0("i = ", i))
    
    for( ii in ii.range ){
      if (ii > 0){
        #print (paste0("ii = ", ii))
        isf <-as.data.frame(X[,ii,drop = FALSE])
        # names(isf) <- names(X)[ii]
        #print (paste0(ii, "   ", names(isf)))
        
        d	<- data.table( dd, isf )
        names(d)[1] <- "NPP"
        temp <- paste(names(d[, - 1]), collapse=paste0( effect.s, ")+s("))
        #print(ii)
        #print(names(d)[1])
        
        form <- formula(paste("NPP ~ ", effect.i, " s(",temp, effect.s, ")"))
      }
      # print(form)
      d <-data.table(d,dt.pb)
      names(d)[1] <- "NPP"
      if (bloc == 0){ 
        
        
        
        result <- tryCatch(
          {
            form.re <-formula(paste("NPP ~ ", effect.i, " s(",temp, effect.s, ")", "+ s(Bloc, bs = 're')"))
            #sfmod0	<- gamm(form, method = 'REML', random = list(Bloc=~1),data = d)$gam 
            sfmod0	<- gam(form.re, method = 'REML',data = d)
          },
          error = function(e) e
        )
        
      }else{
        result <- tryCatch(
          {
            sfmod0	<- gam(form, method = 'REML', data = d)},
          error = function(e) e
        )
        
      }
      
      if (is.null(result$message)) {
        Obja	<-  summary(sfmod0)$r.sq
        p.table <- as.data.frame(summary(sfmod0)$p.table)
        names(p.table)[endsWith("Pr(>|t|)",names(p.table) )]<-"pValue"
        p.table$rn <- row.names(p.table)
        p.table <- setDT(p.table)[, var:=stringr::str_sub(rn, 1,8)][var =="genoType" ]
        if (dim(p.table)[1] > 1)
          pV.i <-metap::sumlog(p.table$pValue)$p
        else pV.i <- NA_real_
        
        s.table <- as.data.frame(summary(sfmod0)$s.table)
        names(s.table)[endsWith("p-value",names(s.table) )]<-"pValue"
        s.table$rn <- row.names(s.table)
        setDT(s.table)[,rn.var:=tstrsplit(rn, ":", fixed=TRUE, keep=1L)][, var:=stringr::str_sub(rn.var, 3,-2)][,var:=tstrsplit(var, ",", fixed=TRUE, keep=1L)]
        s.table<-s.table[var %in% eval(parse(text= paste0("xVar.g", iG)))]
        # pV.s <-metap::sumlog(s.table$pValue)$p
        if (dim(s.table)[1] > 2){
          pV.s<-s.table[,list(pValue.slope =metap::sumlog(s.table$pValue)$p ), by = "var"]
        }else {
          pV.s<-s.table[,c("var","pValue")]
          names(pV.s)[2]<- "pValue.slope"
        }
        
        stats.tmp <- data.table(pV.s
                                , pValue.intercept = pV.i,R2 = Obja, aic=AIC(sfmod0))
        
        
        stats.tmp[, c("bloc", "iG", "i0", "ii") := list(bloc, iG, i, ii)]
        stats.tmp[,selected := 0]
      } 
      # else  Obj	<- -99999
      if (!(is.null(stats.tmp)))
      all.s <- rbind(all.s, stats.tmp)
    }
    
    # all.s.i0 <- setorder(all.s[i0 == i], -R2)
    
    s.i0 <- all.s[i0 == i,]
    s.i0[ is.na(pValue.intercept) & pValue.slope <= 0.05, pre.sel := 1]
    s.i0[ !is.na(pValue.intercept) & pValue.intercept <= 0.05 & pValue.slope <= 0.05, pre.sel := 1]
    i0.05 <- s.i0[ pre.sel == 1,]
    i0.05[ ,pre.sel := NULL]
    #remove the existing
    
    if (dim(i0.05)[1] > 0) {
      tmp.i0 <- i0.05[, .N, by = .(iG, i0, ii, R2, aic)][N == nX.sel+1]
      
      if (dim(tmp.i0)[1] == 0){
        tmp <- i0.05[!(var %in% names(sf_sel_df))]
        if (dim(tmp)[1] > 0) {
          if (sel.by == "R2"){
            setorder(tmp, R2 )
            var.sel <- tmp[, .SD[.N]]}
          if (sel.by == "AIC"){
            setorder(tmp, aic )
            var.sel <- tmp[, .SD[1]]}
          var.sel[,selected.i0 := 1][,selected:= NULL]
          sf_sel_df <- as.data.frame(X.iG[,var.sel$var, drop = FALSE])
          
          Obj_old <- -Inf
          all.s<- var.sel[, c("i0", "ii","var", "selected.i0")][all.s, on = c("i0", "ii","var"), nomatch = NA]
          all.s[selected.i0 == 1, selected := selected.i0]
          all.s[,selected.i0  := NULL]
          X<-X.iG[, names(sf_sel_df), drop = FALSE]
          dd	<- as.data.frame( y, drop = FALSE )
          d <-data.frame(dd, X)
          nX.sel <- dim(sf_sel_df)[2]-1
          var.sel[,nX := nX.sel+1]
          var.sel[,re.sel := 1]
          var.sel <- all.s[, -"selected" ][var.sel[, c("i0", "ii","selected.i0","nX",  "re.sel")], on =c("i0", "ii") ]
          R2.sel <- rbind(R2.sel, var.sel) 
          ii.range <- 0:0
          temp <- paste(names(sf_sel_df), collapse=paste0( effect.s, ")+s("))
          
          #print(ii)
          #print(names(d)[1])
          
          form <- formula(paste("NPP ~ ", effect.i, " s(",temp, effect.s, ")"))
                 # dd <-data.frame(y, sf_sel_df)
          #print(names(dd))
          names(d)[1] <- "NPP"
          if (nX.sel == 0) Obj <- -99999
        } else Obj_old <- Obj
        
        
      }else{
        
        if (sel.by == "R2"){
          setorder(tmp.i0, R2)
          var.sel <- tmp.i0[, .SD[.N]]
          
        }
        
        if (sel.by == "AIC"){
          setorder(tmp.i0, aic)
          var.sel <- tmp.i0[, .SD[1]]
          
        }
        
        Obj <- var.sel$R2
        # chk R2
        if (Obj - Obj_old >=R2.tol){
          var.sel <- var.sel[,-c('R2','aic', 'N')]
          var.sel[,selected.i0 := 1]
          
          
          
          var.sel <- i0.05[,-'selected'][var.sel, on = .(iG, i0, ii)]
          
          sf_sel_df <- as.data.frame(X.iG[,var.sel$var, drop = FALSE])
          nX.sel <- dim(sf_sel_df)[2] 
          if (nX.sel == 1) Obj_old <- -Inf else Obj_old <-  Obj
          var.sel[,nX := nX.sel]
          var.sel[,re.sel := 0]
          setDF(var.sel)
          R2.sel <- rbind(R2.sel, var.sel)
          all.s<- setDT(var.sel)[, c("i0", "ii", "var", "selected.i0")][all.s, on = c("i0", "ii", "var"), nomatch = NA]
          all.s[selected.i0 == 1, selected := selected.i0]
          all.s[,selected.i0  := NULL]
          X<-X.iG[,!(names(X.iG) %in% names(sf_sel_df))]
          dd	<- as.data.frame( cbind(y, sf_sel_df) )
          
          
          ii.range <- 1:dim(X)[2]
        }
        
      }
    }  
    # 
    names(X)
    
    i		<- i + 1
 
  }#end of while (selection)
  if (dim(R2.sel)[1]>0){ 
    setDT(R2.sel)
    setkey(R2.sel, iG, nX)
    R2.sel <- R2.sel[,.SD[.N], by = .(iG, nX)]
    R2.sel[,c("effect.i","effect.s", "resp") := list(effect.i,effect.s, names(y))]
  }
  
  if (dim(sf_sel_df)[2] > 0) { 
    #temp <- paste(names(subset(dd, select=-c(1,2))), collapse=",bs='cr')+s(")
 
    temp <- paste(names(sf_sel_df), collapse=paste0( effect.s, ")+s("))
    
    #print(ii)
    #print(names(d)[1])
    
    form <- formula(paste("NPP ~ ", effect.i, " s(",temp, effect.s, ")"))
 
    dd <-data.frame(y, sf_sel_df)
    
    dd <-data.table(dd,dt.pb)
    #print(names(dd))
    names(dd)[1] <- "NPP"
    if (bloc == 0){
      # print("bloc  0")
    
      form.re <- formula(paste("NPP ~ ", effect.i, " s(",temp, effect.s, ")", "+ s(Bloc, bs = 're')"))
    
      sfmod	<- gam(form.re, method = 'REML',data = dd)
      
    }else{ 
      # print("bloc non 0")
      #print(form)
      sfmod	<- gam(form, method = 'REML', data = dd)}
    
  }else {
    dd <- y
    names(dd)[1] <- "NPP"
    if (bloc == 0){
      print("no variable selected bloc  0")
      
      dd <-data.table(dd,dt.pb[,2, drop = FALSE])
      sfmod	<- gam(NPP~ s(Bloc, bs = 're'), method = 'REML', data = dd)
    } else
      
      #no selection finish the process
      sfmod	<- gam(NPP ~ 1, method = 'REML', data = dd)
  }
  
  return(list(R2.sel = R2.sel, all.s = all.s, dd = dd, sfmod = sfmod, sf_sel_df = sf_sel_df))
}


plotGAM.CI <- function(gamFit, smooth.cov ,
                       groupCovs = NULL) {
  
  
  
  if (missing(gamFit)) {
    stop("gamFit is missing")}
  if (missing(smooth.cov)) {
    stop("smooth.cov is missing")}
  if (class(smooth.cov) != "character") {
    stop("smooth.cov must be a character")}
  if (class(groupCovs) != "character" & (!is.null(groupCovs))) {
    stop("groupCovs must be a character")}
  # if (!(rawOrFitted == F | rawOrFitted == "raw" | rawOrFitted == "fitted")) {
  #   stop("Wrong input for rawOrFitted")
  # }
  
  
  
  fit <- fitted <- group <- se.fit <- NULL
  
  if (base::class(gamFit)[1] != "gam") {
    stop("gamFit is not a object of type gam")
  }
  
  #for groupCovs only
  
  gam <- gamFit
  temp.data <- gam$model
  
  old.formula <- gam$formula
  old.covariates <- base::as.character(old.formula)[3]
  old.covariates <- base::gsub(" ", "", old.covariates)
  main.smooth <- base::paste0("s(", smooth.cov)
  interaction.smooth <- base::paste0("s(", smooth.cov, ",by=")
  
  if(regexpr(groupCovs, as.character(old.covariates))[1] == 1) {
    old.covariates <- paste0("+", old.covariates)
  }
  
  # main.group <- base::paste0("+", groupCovs)
  # 
  # if (!base::grepl(pattern = main.group, old.covariates, fixed=T)) {
  #   base::stop("Error groupCovs is not included as a main effect in original model, include in original fit")
  # }
  # 
  temp.fm <- strsplit(old.covariates, split='+', fixed = T)[[1]]
  for (i in 1:length(temp.fm)) {
    if (grepl(main.smooth, temp.fm[i], fixed=T)) {
      temp.term <- strsplit(temp.fm[i], split=',', fixed = T)[[1]]
      if (length(temp.term) > 1) {
        for (j in 2:length(temp.term)) {
          if (grepl('=',temp.term[j])) {
            if(grepl('by',temp.term[j])) {
              order <- 1:length(temp.term)
              order <- setdiff(order, j)
              if (all(order == 1)) {
                order <- c(1,2)
                temp.term <- temp.term[order]
              }
              else {
                order <- c(order, NA)
                order <- order[c(1, length(order), 2:(length(order) - 1))]
                order[2] <- j
                temp.term <- temp.term[order]
                temp.term <- gsub(pattern=")", replacement = "", x = temp.term)
                temp.term[length(temp.term)]<- paste0(temp.term[length(temp.term)], ")")
              }
            }
          }
          else {
            base::stop("spline containing smooth.cov has more than one variable, refit with only one")
          }
        } #end of for j
      }
      temp.fm[i] <- toString(temp.term)
    }
  } # end of for i
  
  old.covariates <- paste0(temp.fm, collapse="+")
  old.covariates <- base::gsub(" ", "", old.covariates)
  
  
  
  foo <-  base::seq(min(temp.data[smooth.cov]), max(temp.data[smooth.cov]), length.out=200)
  plot.df = base::data.frame( x = rep(foo, nlevels(as.factor(temp.data[groupCovs][,1]))))
  base::names(plot.df) <- smooth.cov
  
  plot.df$temp <- rep(levels(as.factor(temp.data[groupCovs][,1])), each = 200)
  base::names(plot.df)[2] <- groupCovs
  
  for (i in base::names(temp.data)[-1]) {
    if (!(i == smooth.cov | i == groupCovs)) {
      if (base::any(base::class(temp.data[i][,1])[1] == c("numeric", "integer","boolean"))) {
        plot.df[, base::dim(plot.df)[2] + 1] <- base::mean(temp.data[i][,1])
        base::names(plot.df)[base::dim(plot.df)[2]] <- i
      } else if (base::any(base::class(temp.data[i][,1])[1] == c("character", "factor","ordered"))) {
        plot.df[, base::dim(plot.df)[2] + 1] <- temp.data[i][1,1]
        base::names(plot.df)[base::dim(plot.df)[2]] <- i
      } else {
        stop(base::paste("Unrecognized data type for",i,"please refit cluster model with different datatype"))
      }
    }
  }
  
  plot.df = base::cbind(plot.df, base::as.data.frame(mgcv::predict.gam(gam, plot.df, se.fit = TRUE)))
  temp.data$fitted <- gam$fitted.values
  # plot.df$group <- plot.df[,2]
  
  
  plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
           + ggplot2::geom_line(ggplot2::aes(y=fit))
           + ggplot2::geom_ribbon(data=plot.df, ggplot2::aes(ymax = fit+1.96*se.fit, ymin = fit-1.96*se.fit,  linetype=NA), alpha = .2)
           + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], "  vs. s(",smooth.cov,")"))
           + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
           + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))         
           
           + ggplot2::theme_bw())
  
  
  plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x=purrr::as_vector(temp.data[smooth.cov]), color="#0073C2FF", alpha = 0.2, y=temp.data[,1]))
  
  
  plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x=purrr::as_vector(temp.data[smooth.cov]), y=fitted))
  
  plot <- plot + facet_wrap(~genoType)
  plot <- plot + theme(legend.position = "none") 
  
  # return(list(plot.df = plot.df, temp.data = temp.data)  )
  return(plot)
}


gam.model.sel.pValue.Grp.clean <-function (y, N.grp = 1, bloc,resp.label, effect.i , effect.s ){
  
  
  
  # df<-y
  # # forT <- paste0("Y.",names(y))
  # sf_sel_df <- data.frame()
  # # print(forT)
  # 
  # # Obj_old <- -Inf; 
  # Obj <- -99999
  # 
  # # obj	<- Obj
  # nX.sel <- 0
  # all.s <- data.frame()
  # # selected.s <- data.table()
  # 
  #    i	<- 0
  #   Obj_old <- -Inf
  #   # print(iG)
  #   X<- as.data.frame( dt.use[,xVar, with = FALSE])
  #   X.iG <- X
  #   nX <- dim(X)[2]
  #   ii.range <- 1:nX
  #   # sf_list	<- 1:nX
  #   R2.sel <- data.frame()
  sel <- fun_var.sel(y,effect.i = effect.i, effect.s = effect.s, max.Var = 2, R2.tol = 0.01, bloc = 0)
  sf_sel_df<-sel$sf_sel_df
  sfmod <-sel$sfmod
  dd <- sel$dd
  R2.sel<-sel$R2.sel
  all.s<-sel$all.s
  # rm(all.s, R2.sel, sfmod, sf_sel_df,sel)
  #graphing
  if (dim(sf_sel_df)[1] > 0 ){
    P.vars <- names(sf_sel_df)
    
    map(P.vars, function(x){
      
      p <- plotGAM.CI(sfmod, smooth.cov = x, groupCovs = "genoType") +
        #  geom_rug(data = train.college, aes_string(y = "Outstate", x = x ), alpha = 0.2)
        ylab(paste0(names(y))) +
        # # ggtitle()
        ggtitle(paste0(names(y), " vs. s(", x,")"))
      
      g <- ggplotGrob(p)
    }) %>%
    {grid.arrange(grobs = (.), ncol = 2)}
    
    
    myplot <- recordPlot()
    graphics.off()
    # replayPlot(myplot)
  }else myplot <-NULL
  #  print(is.null(sfmod$lme))
  #2+ linear predictors
  
  
  pred	<- predict( sfmod )
  resid	<- residuals ( sfmod)
  
  aic <- AIC(sfmod)
  # obs = pred+resid
  sfmod.summary <-summary(sfmod)
  
  r.sq = sfmod.summary$r.sq
  R.sq = data.frame(r.sq,aic)
  names(R.sq) <- c('adjRsq', 'aic')
  
  p_table <- data.frame(sfmod.summary$p.table)
  p_table <- within(p_table, {lci <- Estimate - qnorm(0.975) * Std..Error
  uci <- Estimate + qnorm(0.975) * Std..Error})
  
  p.r <-data.frame(dd[,1],pred,resid)
  
  names(p.r)[1:3] <-c(resp.label,'pred','resid')
  p.r <-data.frame(dt.pb, p.r)
  
  #  p_table
  
  # Nsf <- dim(sf_sel_df)[2] 
  if (dim(sf_sel_df)[2] > 0)
    p.r <-data.frame(p.r, sf_sel_df)
  # p.r$var <- row.names(p.r)
  names.na <- setdiff(xchr.complete, as.list(names(sf_sel_df)) )
  p.r.na <-as.data.frame(matrix(NA_real_, nrow = dim(p.r)[1], ncol = length(names.na)))
  names(p.r.na) <- names.na
  p.r <-data.frame(p.r, p.r.na)
  #setdiff(xchr.complete,c("Geno2_p" , "Geno3_p" ))
  s_table <-data.frame(sfmod.summary$s.table)
  # setDF(all.s)
  # setDF(R2.sel)
  if (dim(R2.sel)[1]  > 0){
    print("BBB")
    print(names(y))
    
    R2.sel$resp <- names(y)}
  if (dim(s_table)[1]  > 0){
    s_table$resp <- names(y)}
  
  
  p_table$resp <- names(y)
  all.s$resp <- names(y)
  p.r$resp <- names(y)
  
  #  s_table
  # hist(y)
  list(p.table = p_table, R.sq = R2.sel, s.table = s_table, Pred = p.r, sequential.sel = all.s, plots = myplot)
}
